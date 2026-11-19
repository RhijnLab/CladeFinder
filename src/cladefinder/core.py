"""VCF/metadata handling and transparent SNP-distance classification."""
import csv
import gzip
import hashlib
import heapq
import json
import re
from collections import Counter
from pathlib import Path

BASES = set('ACGT')

def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path, 'rt')

def metadata(path):
    with open(path, newline='', encoding='utf-8-sig') as fh:
        rows = list(csv.DictReader(fh))
    if not rows or not {'ID', 'Cluster'}.issubset(rows[0]):
        raise ValueError('Metadata must contain ID and Cluster columns')
    out = {}
    for r in rows:
        ident, cl = r['ID'].strip(), r['Cluster'].strip()
        if ident and cl:
            if ident in out and out[ident] != cl:
                raise ValueError(f'Conflicting clusters for {ident}')
            out[ident] = cl
    return out

def nexus_tips(path):
    """Extract leaf labels from annotated BEAST-style Nexus tree; no topology inference."""
    text = Path(path).read_text(encoding='utf-8')
    match = re.search(r'\btree\s+[^=]+\s*=\s*(?:\[&R\]\s*)?(\(.*?;)', text, re.I | re.S)
    if not match:
        raise ValueError('No Newick tree found in Nexus file')
    tree = match.group(1)
    # Remove bracket annotations (nested braces inside annotations are OK).
    tree = re.sub(r'\[[^\[\]]*\]', '', tree)
    leaves = re.findall(r'(?<=[(,])\s*([^():,;\s]+)\s*(?=[:),])', tree)
    if len(leaves) != len(set(leaves)) or not leaves:
        raise ValueError('Nexus tips missing or duplicated; inspect tree formatting')
    return leaves

def resolve(label, ids):
    if label in ids:
        return label
    # BEAST tree labels can append a sampling year, e.g. C103_2015.
    stripped = re.sub(r'_\d{4}$', '', label)
    if stripped in ids:
        return stripped
    raise ValueError(f'Cannot match isolate {label!r} to metadata')

def genotype_base(gt, ref, alt):
    if gt in ('.', './.', '.|.', ''):
        return None
    alleles = re.split(r'[|/]', gt)
    if not alleles or any(x == '.' for x in alleles):
        return None
    if len(set(alleles)) != 1:
        return None  # heterozygous calls excluded for haploid-distance model
    try:
        idx = int(alleles[0])
    except ValueError:
        return None
    candidates = [ref] + alt.split(',')
    if idx >= len(candidates) or len(candidates[idx]) != 1:
        return None
    return candidates[idx].upper() if candidates[idx].upper() in BASES else None

def selected_sites(vcf, tree_tips, metadata_ids, max_sites, min_called_fraction=0.98):
    """Stream uncompressed VCF once; deterministic bounded-memory hash sampling."""
    if max_sites < 1:
        raise ValueError('--max-sites must be positive')
    selected = []
    header = None
    columns = None
    total_eligible = 0
    with open_text(vcf) as fh:
        for line in fh:
            if line.startswith('##'):
                continue
            if line.startswith('#CHROM'):
                header = line.rstrip('\n').split('\t')
                sample_to_col = {}
                for i, name in enumerate(header[9:], 9):
                    try:
                        canonical = resolve(name, metadata_ids)
                    except ValueError:
                        continue
                    if canonical in sample_to_col:
                        raise ValueError(f'Duplicate VCF samples for {canonical}')
                    sample_to_col[canonical] = i
                columns = []
                for tip in tree_tips:
                    canonical = resolve(tip, metadata_ids)
                    if canonical not in sample_to_col:
                        raise ValueError(f'Tree isolate {tip} ({canonical}) absent from VCF')
                    columns.append(sample_to_col[canonical])
                continue
            if line.startswith('#'):
                continue
            if columns is None:
                raise ValueError('VCF #CHROM header not found')
            parts = line.rstrip('\n').split('\t')
            if len(parts) <= max(columns):
                raise ValueError('Malformed VCF record with missing sample columns')
            chrom, pos, _, ref, alt = parts[:5]
            if len(ref) != 1 or ref.upper() not in BASES or not alt or alt == '.':
                continue
            alts = alt.split(',')
            if any(len(a) != 1 or a.upper() not in BASES for a in alts):
                continue
            fmt = parts[8].split(':')
            if 'GT' not in fmt:
                continue
            gt_idx = fmt.index('GT')
            calls = []
            for i in columns:
                fields = parts[i].split(':')
                gt = fields[gt_idx] if len(fields) > gt_idx else '.'
                calls.append(genotype_base(gt, ref, alt) or 'N')
            if sum(b != 'N' for b in calls) / len(calls) < min_called_fraction:
                continue
            if len(set(b for b in calls if b != 'N')) < 2:
                continue
            total_eligible += 1
            key = f'{chrom}:{pos}:{ref}:{alt}'.encode()
            score = int.from_bytes(hashlib.blake2b(key, digest_size=8).digest(), 'big')
            item = (-score, chrom, int(pos), ref.upper(), alt.upper(), ''.join(calls))
            if len(selected) < max_sites:
                heapq.heappush(selected, item)
            elif item > selected[0]:
                heapq.heapreplace(selected, item)
    if not selected:
        raise ValueError('No suitable SNPs found. Check VCF samples, genotypes and filters')
    sites = sorted([(c, p, r, a, calls) for _, c, p, r, a, calls in selected], key=lambda x:(x[0],x[1]))
    if len({(c,p) for c,p,*_ in sites}) != len(sites):
        raise ValueError('Multiple VCF records at a selected position; normalize VCF first')
    return sites, total_eligible

def distances(query, references, clusters, min_compared=100, min_margin=0.005):
    """Missing-aware nearest-reference p-distance; no phylogenetic placement.

    A nearest cluster is descriptive even when the assignment is ambiguous.
    A competing-cluster margin is evaluated independently of the nearest
    within-cluster reference. Thresholds are heuristics, not calibrated support.
    """
    scored = []
    for tip, seq in references.items():
        if len(query) != len(seq):
            raise ValueError('Sequence lengths do not match')
        comparable = 0
        mismatches = 0
        for a, b in zip(query, seq):
            if a in BASES and b in BASES:
                comparable += 1
                mismatches += (a != b)
        if comparable >= min_compared:
            scored.append({'isolate': tip, 'cluster': clusters[tip],
                           'mismatches': mismatches, 'compared_sites': comparable,
                           'p_distance': mismatches / comparable})
    scored.sort(key=lambda r: (r['p_distance'], r['mismatches'], r['isolate']))
    if not scored:
        return {'status': 'insufficient_data', 'predicted_cluster': None,
                'nearest_cluster': None, 'nearest': None, 'second_cluster': None,
                'second_cluster_distance': None, 'margin': None,
                'tied_clusters': [], 'top_matches': [], 'cluster_best': [],
                'eligible_references': 0, 'ambiguity_reason': 'no_reference_meets_min_compared'}
    best = scored[0]
    competing = next((r for r in scored if r['cluster'] != best['cluster']), None)
    margin = competing['p_distance'] - best['p_distance'] if competing else None
    # Floating point tolerance is only for detecting genuine distance ties.
    tied = sorted({r['cluster'] for r in scored
                   if abs(r['p_distance'] - best['p_distance']) <= 1e-12})
    if competing is None:
        status, reason = 'ambiguous', 'no_competing_cluster_reference'
    elif len(tied) > 1:
        status, reason = 'ambiguous', 'equal_nearest_distance_across_clusters'
    elif margin < min_margin:
        status, reason = 'ambiguous', 'margin_below_threshold'
    else:
        status, reason = 'assigned', None
    by_cluster = {}
    for r in scored:
        by_cluster.setdefault(r['cluster'], r)
    return {'status': status,
            'predicted_cluster': best['cluster'] if status == 'assigned' else None,
            'nearest_cluster': best['cluster'], 'nearest': best,
            'second_cluster': competing['cluster'] if competing else None,
            'second_cluster_distance': competing['p_distance'] if competing else None,
            'margin': margin, 'tied_clusters': tied,
            'top_matches': scored[:10],
            'cluster_best': sorted(by_cluster.values(), key=lambda r: (r['p_distance'], r['cluster'])),
            'eligible_references': len(scored), 'ambiguity_reason': reason}

def load_database(path):
    with open(Path(path)/'database.json') as fh:
        return json.load(fh)
