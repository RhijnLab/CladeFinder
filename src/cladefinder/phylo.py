"""Nexus topology diagnostics and EPA-ng phylogenetic placement.

Tree support is separate from SNP nearest-neighbour evidence. BEAST time-scale
branch lengths are NEVER passed directly to EPA-ng: optimize on alignment first.
"""
import csv
import io
import json
import re
import subprocess
import tempfile
import os
import shutil
from collections import Counter
from pathlib import Path
from Bio import Phylo
from .core import load_database


def load_tree(path):
    tree = Phylo.read(str(path), 'nexus')
    tips = [x.name for x in tree.get_terminals()]
    if not tips or len(tips) != len(set(tips)):
        raise ValueError('Tree has missing or duplicate tip labels')
    return tree


def tree_report(database, output):
    db = load_database(database)
    tree = load_tree(Path(database) / 'reference.nex')
    labels = db['clusters']
    if set(labels) != {x.name for x in tree.get_terminals()}:
        raise ValueError('Nexus tree tips and database isolates differ')
    total = Counter(labels.values())
    rows = []
    for cluster in sorted(total):
        tips = [x for x in tree.get_terminals() if labels[x.name] == cluster]
        ancestor = tree.common_ancestor(tips)
        descendants = [x.name for x in ancestor.get_terminals()]
        foreign = Counter(labels[x] for x in descendants if labels[x] != cluster)
        # Number of maximal pure subtrees for this cluster.
        pure_roots = []
        def walk(node, parent_pure=False):
            subtree = [x.name for x in node.get_terminals()]
            pure = all(labels[x] == cluster for x in subtree)
            if pure and not parent_pure:
                pure_roots.append(node)
                return
            for child in node.clades:
                walk(child, pure)
        walk(tree.root)
        rows.append({'cluster':cluster,'tip_count':total[cluster],
                     'mrca_descendant_count':len(descendants),
                     'foreign_tips_in_mrca':sum(foreign.values()),
                     'monophyletic':sum(foreign.values()) == 0,
                     'pure_subtree_count':len(pure_roots),
                     'foreign_clusters':','.join(f'{k}:{v}' for k,v in sorted(foreign.items()))})
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    with (out/'tree_clusters.tsv').open('w',newline='') as fh:
        w=csv.DictWriter(fh,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    summary={'tip_count':len(labels),'cluster_report':rows,
             'warning':'Nexus branch lengths are not automatically substitutions/site. Tree topology alone cannot place a new sequence.'}
    (out/'tree_report.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary


def _write_verified(destination, source):
    """Copy a finished local artifact onto a potentially unreliable Docker Desktop mount.

    Never use NamedTemporaryFile/rename on the bind mount. Docker Desktop file
    sharing has produced transient ENOENT for such files in this workflow.
    """
    destination = Path(destination)
    source = Path(source)
    if not source.is_file() or source.stat().st_size == 0:
        raise RuntimeError(f"Missing or empty generated artifact: {source}")
    with source.open('rb') as reader, destination.open('wb') as writer:
        shutil.copyfileobj(reader, writer, length=1024 * 1024)
        writer.flush()
        os.fsync(writer.fileno())
    if destination.stat().st_size != source.stat().st_size:
        raise RuntimeError(f"Incomplete copy to {destination}")


def _stage_reference(database, stage):
    """Generate inputs inside container-local scratch, never on mounted database."""
    db = load_database(database)
    refs = db['references']
    if not refs:
        raise ValueError('No reference sequences in database.json')
    lengths = {len(s) for s in refs.values()}
    if len(lengths) != 1 or 0 in lengths:
        raise ValueError(f'Invalid reference alignment lengths: {lengths}')
    expected = db.get('site_count')
    if expected is not None and next(iter(lengths)) != expected:
        raise ValueError(f'Alignment length does not match site_count {expected}')
    tree = load_tree(Path(database) / 'reference.nex')
    if set(refs) != {tip.name for tip in tree.get_terminals()}:
        raise ValueError('Reference alignment IDs do not match Nexus tree tips')
    alignment = stage / 'reference_snps.fasta'
    with alignment.open('w') as fh:
        for tip, seq in refs.items():
            fh.write(f'>{tip}\n{seq}\n')
    for node in tree.find_clades():
        node.branch_length = None
        node.comment = None
    topology = stage / 'nexus_topology.nwk'
    Phylo.write(tree, str(topology), 'newick')
    if alignment.stat().st_size == 0 or topology.stat().st_size == 0:
        raise RuntimeError('Empty staged alignment or topology')
    return alignment, topology


def write_alignment(database):
    """Safely regenerate alignment and topology from database.json."""
    folder = Path(database)
    with tempfile.TemporaryDirectory(prefix='cladefinder_stage_') as tmp:
        alignment, topology = _stage_reference(folder, Path(tmp))
        _write_verified(folder / 'reference_snps.fasta', alignment)
        _write_verified(folder / 'nexus_topology.nwk', topology)
    return folder / 'reference_snps.fasta'


def prepare_placement(database, threads=4):
    """Fit the fixed Nexus topology in local scratch and publish verified outputs.

    A failure in RAxML never touches existing placement files. Publication
    happens only after both fitted outputs exist and are readable.
    """
    folder = Path(database).resolve()
    if not (folder / 'database.json').is_file():
        raise FileNotFoundError(f'Database missing: {folder / "database.json"}')
    if threads < 1:
        raise ValueError('threads must be >= 1')
    with tempfile.TemporaryDirectory(prefix='cladefinder_fit_') as tmp:
        stage = Path(tmp)
        alignment, topology = _stage_reference(folder, stage)
        prefix = stage / 'fit'
        cmd = ['raxml-ng', '--evaluate', '--msa', str(alignment),
               '--tree', str(topology), '--model', 'GTR+G+ASC_LEWIS',
               '--prefix', str(prefix), '--threads', str(threads)]
        subprocess.run(cmd, check=True)
        best = stage / 'fit.raxml.bestTree'
        model = stage / 'fit.raxml.bestModel'
        if not best.is_file() or not model.is_file() or not best.stat().st_size or not model.stat().st_size:
            raise RuntimeError('RAxML-NG did not produce a nonempty bestTree and bestModel')
        # Validate Newick and tip labels before publishing.
        fitted = Phylo.read(str(best), 'newick')
        reference_tips = {x.name for x in load_tree(folder / 'reference.nex').get_terminals()}
        if {x.name for x in fitted.get_terminals()} != reference_tips:
            raise RuntimeError('Fitted tree tip labels differ from Nexus tree')
        manifest = {
            'method': 'EPA-ng fixed-reference placement (experimental)',
            'alignment': 'reference_snps.fasta',
            'tree': 'placement_tree.nwk',
            'fitted_model': 'placement_model.txt',
            'raxml_model': 'GTR+G+ASC_LEWIS',
            'epa_model': 'GTR+G',
            'validated': False,
            'caveat': 'RAxML ascertainment correction and EPA-ng model differ; SNP-only likelihood weights are not calibrated cluster probabilities.'
        }
        manifest_stage = stage / 'placement_manifest.json'
        manifest_stage.write_text(json.dumps(manifest, indent=2) + '\n')
        # Publish manifest LAST; it is the completion marker.
        for source, name in [(alignment, 'reference_snps.fasta'),
                             (topology, 'nexus_topology.nwk'),
                             (best, 'placement_tree.nwk'),
                             (model, 'placement_model.txt'),
                             (manifest_stage, 'placement_manifest.json')]:
            _write_verified(folder / name, source)
    return str(folder / 'placement_tree.nwk')


def filter_variable_sites(reference_sequences, query_sequence):
    """Remove columns invariant among *training* references only.

    Apply precisely the same retained indices to the withheld query. Columns
    with ambiguous reference states are excluded to avoid treating N as an
    alternative allele. The withheld query is never used to select sites.
    """
    if not reference_sequences:
        raise ValueError('No training reference sequences')
    lengths = {len(s) for s in reference_sequences.values()}
    if len(lengths) != 1 or len(query_sequence) != next(iter(lengths)):
        raise ValueError('Reference/query SNP alignment lengths disagree')
    length = next(iter(lengths))
    if length == 0:
        raise ValueError('Empty SNP alignment')
    names = list(reference_sequences)
    sequences = [reference_sequences[name].upper() for name in names]
    allowed = set('ACGT')
    keep = []
    for i in range(length):
        states = {seq[i] for seq in sequences}
        if states <= allowed and len(states) >= 2:
            keep.append(i)
    if not keep:
        raise ValueError('No variable A/C/G/T sites remain after withholding the isolate')
    filtered = {name: ''.join(sequences[j][i] for i in keep)
                for j, name in enumerate(names)}
    query = query_sequence.upper()
    filtered_query = ''.join(query[i] for i in keep)
    return filtered, filtered_query, {
        'original_sites': length,
        'retained_sites': len(keep),
        'removed_invariant_or_ambiguous_sites': length-len(keep),
    }


def _edge_clusters(jplace_tree, labels):
    # jplace edge IDs are encoded as {123} on each branch.
    annotated=re.sub(r'\{(\d+)\}',r'[edge_\1]',jplace_tree)
    tree=Phylo.read(io.StringIO(annotated),'newick')
    edges={}
    all_names=set(labels)
    for node in tree.find_clades():
        m=re.search(r'edge_(\d+)',node.comment or '')
        if not m: continue
        desc={tip.name for tip in node.get_terminals()}
        if not desc <= all_names:
            raise ValueError('EPA-ng tree includes unknown reference tips')
        side1=Counter(labels[x] for x in desc)
        side2=Counter(labels[x] for x in all_names-desc)
        # A pure side indicates a cluster-defining split. Mixed edges are unresolved.
        pure=[]
        for side in (side1,side2):
            if len(side)==1 and side:
                pure.append(next(iter(side)))
        edges[int(m.group(1))]=next(iter(side1)) if len(side1)==1 else (next(iter(side2)) if len(side2)==1 else None)
    return edges


def placement_from_jplace(jplace_path, clusters, min_lwr=0.8):
    data=json.loads(Path(jplace_path).read_text())
    fields=data['fields']
    if 'edge_num' not in fields or 'like_weight_ratio' not in fields:
        raise ValueError('jplace missing edge_num or like_weight_ratio')
    idx_edge=fields.index('edge_num');idx_lwr=fields.index('like_weight_ratio')
    edges=_edge_clusters(data['tree'],clusters)
    if len(data['placements']) != 1:
        raise ValueError('Expected exactly one query placement')
    p=data['placements'][0]
    support=Counter(); unresolved=0.0
    for entry in p['p']:
        cluster=edges.get(int(entry[idx_edge])); weight=float(entry[idx_lwr])
        if cluster is None: unresolved+=weight
        else: support[cluster]+=weight
    ranked=sorted(support.items(),key=lambda x:(-x[1],x[0]))
    best=ranked[0] if ranked else None
    assigned=bool(best and best[1]>=min_lwr and unresolved<=1-min_lwr+1e-9)
    return {'status':'assigned' if assigned else 'ambiguous',
            'predicted_cluster':best[0] if assigned else None,
            'nearest_tree_cluster':best[0] if best else None,
            'cluster_likelihood_weights':dict(support),
            'unresolved_edge_weight':unresolved,
            'minimum_weight':min_lwr,
            'note':'Weights are conditional EPA-ng placement likelihood weights, not calibrated cluster probabilities.'}


def place_sequence(database, sequence, output, threads=4, min_lwr=0.8):
    db=load_database(database); folder=Path(database)
    if len(sequence)!=db['site_count']:
        raise ValueError('Query sequence length differs from database sites')
    for name in ('reference_snps.fasta','placement_tree.nwk'):
        if not (folder/name).exists():
            raise ValueError('Placement files missing: run prepare-placement --database first')
    out=Path(output);out.mkdir(parents=True,exist_ok=True)
    query=out/'query_snps.fasta';query.write_text('>query\n'+sequence+'\n')
    cmd=['epa-ng','--ref-msa',str(folder/'reference_snps.fasta'),
         '--tree',str(folder/'placement_tree.nwk'),'--query',str(query),
         '--model','GTR+G','--threads',str(threads),'--outdir',str(out)]
    subprocess.run(cmd,check=True)
    jplace=out/'epa_result.jplace'
    if not jplace.exists(): raise RuntimeError('EPA-ng did not create epa_result.jplace')
    return placement_from_jplace(jplace,db['clusters'],min_lwr)
