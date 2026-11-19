import argparse
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path
from .core import metadata, nexus_tips, resolve, selected_sites, distances, load_database, genotype_base

def run(cmd, stdout=None):
    print('+ ' + ' '.join(map(str, cmd)), file=sys.stderr)
    subprocess.run(list(map(str,cmd)), check=True, stdout=stdout)

def build(args):
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    labels = metadata(args.clusters)
    tips = nexus_tips(args.tree)
    mapped = {tip: labels[resolve(tip, labels)] for tip in tips}
    sites, eligible = selected_sites(args.vcf, tips, labels, args.max_sites, args.min_called_fraction)
    if len(set(mapped.values())) < 2:
        raise ValueError('Need at least two distinct clusters')
    # Validate contig names against FASTA, avoiding sequence load.
    fasta_contigs = set()
    with open(args.reference) as fh:
        for line in fh:
            if line.startswith('>'):
                fasta_contigs.add(line[1:].split()[0])
    missing = sorted(set(s[0] for s in sites)-fasta_contigs)
    if missing:
        raise ValueError(f'VCF contig names absent from FASTA: {missing[:10]}')
    shutil.copyfile(args.reference, out/'reference.fna')
    shutil.copyfile(args.tree, out/'reference.nex')
    with open(out/'sites.bed','w') as fh:
        for chrom,pos,ref,alt,calls in sites:
            fh.write(f'{chrom}\t{pos-1}\t{pos}\n')
    with open(out/'sites.tsv','w') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['CHROM','POS','REF','ALT'])
        for chrom,pos,ref,alt,calls in sites:
            writer.writerow([chrom,pos,ref,alt])
    refs = {tip: ''.join(site[4][i] for site in sites) for i,tip in enumerate(tips)}
    database = {'schema_version':1, 'method':'nearest-reference SNP p-distance (not phylogenetic placement)',
                'clusters':mapped,'references':refs,'site_count':len(sites), 'eligible_site_count':eligible,
                'tree_tip_count':len(tips),'max_sites':args.max_sites}
    (out/'database.json').write_text(json.dumps(database, separators=(',',':')))
    print(json.dumps({'database':str(out),'sites':len(sites),'eligible_sites':eligible,'tree_tips':len(tips),'clusters':sorted(set(mapped.values()))},indent=2))

def parse_calls(text, site_lookup, min_depth):
    calls = {}
    for line in text.splitlines():
        parts = line.split('\t')
        if len(parts) < 6:
            continue
        chrom,pos,ref,alt,gt,dp = parts[:6]
        key = (chrom,int(pos))
        if key not in site_lookup:
            continue
        try:
            depth = int(dp)
        except ValueError:
            continue
        if depth < min_depth:
            continue
        base = genotype_base(gt,ref,alt)
        if base:
            calls[key] = base
    return calls

def classify(args):
    db = load_database(args.database)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    reference = Path(args.database)/'reference.fna'
    bed = Path(args.database)/'sites.bed'
    bam = out/'sample.sorted.bam'
    # No shell interpolation: paths with spaces work safely.
    run(['bwa','index',reference]) if not Path(str(reference)+'.bwt').exists() else None
    bwa = subprocess.Popen(['bwa','mem','-t',str(args.threads),str(reference),args.r1]+([args.r2] if args.r2 else []),stdout=subprocess.PIPE)
    sort = subprocess.run(['samtools','sort','-@',str(args.threads),'-o',str(bam),'-'],stdin=bwa.stdout)
    bwa.stdout.close()
    bwa_rc = bwa.wait()
    if bwa_rc or sort.returncode:
        raise RuntimeError(f'Mapping failed (bwa={bwa_rc}, samtools={sort.returncode})')
    run(['samtools','index',bam])
    # Emit all targeted sites, including invariant reference calls, with depth.
    mp = subprocess.Popen(['bcftools','mpileup','-Ou','-f',str(reference),'-R',str(bed),'-a','FORMAT/DP',str(bam)],stdout=subprocess.PIPE)
    call = subprocess.Popen(['bcftools','call','-m','-A','-Ou'],stdin=mp.stdout,stdout=subprocess.PIPE)
    mp.stdout.close()
    query = subprocess.run(['bcftools','query','-f','%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%DP]\n'],stdin=call.stdout,capture_output=True,text=True)
    call.stdout.close()
    call_rc = call.wait()
    mp_rc = mp.wait()
    if query.returncode or call_rc or mp_rc:
        raise RuntimeError(f'Variant genotyping failed: {query.stderr}')
    with open(Path(args.database)/'sites.tsv') as fh:
        rows = list(csv.DictReader(fh,delimiter='\t'))
    lookup = {(r['CHROM'],int(r['POS'])) for r in rows}
    called = parse_calls(query.stdout,lookup,args.min_depth)
    seq = ''.join(called.get((r['CHROM'],int(r['POS'])),'N') for r in rows)
    result = distances(seq,db['references'],db['clusters'],args.min_compared,args.min_margin)
    if args.phylogeny != 'off':
        from .phylo import place_sequence
        try:
            result['phylogenetic_placement'] = place_sequence(args.database, seq, out/'placement', args.threads, args.min_lwr)
        except (ValueError, RuntimeError, subprocess.CalledProcessError) as exc:
            if args.phylogeny == 'required':
                raise
            result['phylogenetic_placement'] = {'status':'unavailable','reason':str(exc)}
    result.update({'sample':args.sample or Path(args.r1).name,'method':db['method'],
                   'called_sites':len(called),'total_sites':len(rows),'call_fraction':len(called)/len(rows),
                   'min_depth':args.min_depth,'min_compared':args.min_compared,'min_margin':args.min_margin})
    (out/'classification.json').write_text(json.dumps(result,indent=2)+'\n')
    with open(out/'classification.tsv','w') as fh:
        fields = ['sample','status','predicted_cluster','called_sites','total_sites','call_fraction','nearest_isolate','nearest_cluster','second_cluster','second_cluster_distance','mismatches','compared_sites','p_distance','margin','ambiguity_reason','tied_clusters','eligible_references']
        writer = csv.DictWriter(fh,fieldnames=fields,delimiter='\t'); writer.writeheader()
        nearest = result['nearest'] or {}
        writer.writerow({'sample':result['sample'],'status':result['status'],'predicted_cluster':result['predicted_cluster'] or '',
                         'called_sites':len(called),'total_sites':len(rows),'call_fraction':result['call_fraction'],
                         'nearest_isolate':nearest.get('isolate',''),'nearest_cluster':nearest.get('cluster',''),
                         'mismatches':nearest.get('mismatches',''),'compared_sites':nearest.get('compared_sites',''),
                         'p_distance':nearest.get('p_distance',''),'margin':result['margin'] if result['margin'] is not None else '',
                         'second_cluster':result['second_cluster'] or '',
                         'second_cluster_distance':result['second_cluster_distance'] if result['second_cluster_distance'] is not None else '',
                         'ambiguity_reason':result['ambiguity_reason'] or '',
                         'tied_clusters':','.join(result['tied_clusters']),
                         'eligible_references':result['eligible_references']})
    if not args.keep_bam:
        bam.unlink(missing_ok=True); Path(str(bam)+'.bai').unlink(missing_ok=True)
    print(json.dumps(result,indent=2))

def validate(args):
    from collections import Counter, defaultdict
    db = load_database(args.database)
    refs, clusters = db['references'], db['clusters']
    results = []
    for i, (tip, seq) in enumerate(refs.items(), 1):
        remaining = {k: v for k, v in refs.items() if k != tip}
        res = distances(seq, remaining, clusters, args.min_compared, args.min_margin)
        nearest = res['nearest'] or {}
        results.append({'isolate': tip, 'true_cluster': clusters[tip],
                        'predicted_cluster': res['predicted_cluster'] or '',
                        'nearest_cluster': res['nearest_cluster'] or '',
                        'nearest_isolate': nearest.get('isolate', ''),
                        'nearest_mismatches': nearest.get('mismatches', ''),
                        'nearest_compared_sites': nearest.get('compared_sites', ''),
                        'nearest_p_distance': nearest.get('p_distance', ''),
                        'second_cluster': res['second_cluster'] or '',
                        'second_cluster_distance': res['second_cluster_distance'] if res['second_cluster_distance'] is not None else '',
                        'margin': res['margin'] if res['margin'] is not None else '',
                        'tied_clusters': ','.join(res['tied_clusters']),
                        'eligible_references': res['eligible_references'],
                        'ambiguity_reason': res['ambiguity_reason'] or '',
                        'status': res['status'],
                        'correct': res['predicted_cluster'] == clusters[tip]})
        if i % 50 == 0:
            print(f'Validated {i}/{len(refs)} isolates', file=sys.stderr)
    out = Path(args.output)
    out.mkdir(parents=True, exist_ok=True)
    with open(out/'leave_one_out.tsv', 'w') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(results[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(results)
    cluster_counts = defaultdict(Counter)
    for row in results:
        cluster_counts[row['true_cluster']]['total'] += 1
        cluster_counts[row['true_cluster']]['correct'] += bool(row['correct'])
        cluster_counts[row['true_cluster']]['assigned'] += row['status'] == 'assigned'
        cluster_counts[row['true_cluster']][row['status']] += 1
    with open(out/'cluster_summary.tsv', 'w') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['cluster', 'total', 'assigned', 'correct', 'ambiguous', 'insufficient_data', 'correct_fraction', 'assignment_fraction'])
        for cluster in sorted(cluster_counts):
            x = cluster_counts[cluster]
            writer.writerow([cluster, x['total'], x['assigned'], x['correct'], x['ambiguous'], x['insufficient_data'],
                             x['correct']/x['total'], x['assigned']/x['total']])
    with open(out/'confusion_matrix.tsv', 'w') as fh:
        writer = csv.writer(fh, delimiter='\t')
        writer.writerow(['true_cluster', 'predicted_cluster', 'count'])
        for (true, pred), count in sorted(Counter((r['true_cluster'], r['predicted_cluster'] or 'unassigned') for r in results).items()):
            writer.writerow([true, pred, count])
    assigned = sum(r['status'] == 'assigned' for r in results)
    correct = sum(r['correct'] for r in results)
    summary = {'total': len(results), 'correct': correct, 'assigned': assigned,
               'accuracy_including_abstentions': correct/len(results),
               'accuracy_among_assigned': correct/assigned if assigned else None,
               'assignment_fraction': assigned/len(results),
               'unassigned': len(results)-assigned,
               'ambiguity_reasons': dict(Counter(r['ambiguity_reason'] for r in results if r['ambiguity_reason'])),
               'min_compared': args.min_compared, 'min_margin': args.min_margin,
               'note': 'Internal leave-one-out check only. SNP sites were selected using all isolates. Not independent or phylogenetic-placement validation.'}
    (out/'validation.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))

def main():
    p=argparse.ArgumentParser(prog='cladefinder',description='FASTQ to reference-cluster SNP-distance classifier')
    sub=p.add_subparsers(dest='command',required=True)
    b=sub.add_parser('build',help='Build a database locally from VCF, Nexus tree, metadata and FASTA')
    for opt in ('vcf','reference','tree','clusters','output'): b.add_argument('--'+opt,required=True)
    b.add_argument('--max-sites',type=int,default=20000)
    b.add_argument('--min-called-fraction',type=float,default=0.98)
    b.set_defaults(func=build)
    c=sub.add_parser('classify',help='Map FASTQ and classify')
    for opt in ('r1','database','output'):c.add_argument('--'+opt,required=True)
    c.add_argument('--r2');c.add_argument('--sample');c.add_argument('--threads',type=int,default=4)
    c.add_argument('--min-depth',type=int,default=5);c.add_argument('--min-compared',type=int,default=100)
    c.add_argument('--min-margin',type=float,default=0.005);c.add_argument('--keep-bam',action='store_true')
    c.add_argument('--phylogeny',choices=['off','auto','required'],default='auto')
    c.add_argument('--min-lwr',type=float,default=0.8)
    c.set_defaults(func=classify)
    v=sub.add_parser('validate',help='Leave-one-out nearest-reference validation')
    v.add_argument('--database',required=True);v.add_argument('--output',required=True)
    v.add_argument('--min-compared',type=int,default=100);v.add_argument('--min-margin',type=float,default=0.005)
    v.set_defaults(func=validate)
    t=sub.add_parser('tree-report',help='Assess cluster monophyly on the supplied Nexus topology')
    t.add_argument('--database',required=True);t.add_argument('--output',required=True)
    def report_command(a):
        from .phylo import tree_report
        print(json.dumps(tree_report(a.database,a.output),indent=2))
    t.set_defaults(func=report_command)
    prep=sub.add_parser('prepare-placement',help='Reconstruct SNP alignment and optimize Nexus topology branch lengths')
    prep.add_argument('--database',required=True);prep.add_argument('--threads',type=int,default=4)
    def prep_command(a):
        from .phylo import prepare_placement
        print(prepare_placement(a.database,a.threads))
    prep.set_defaults(func=prep_command)
    vp=sub.add_parser('validate-placement',help='Leave-one-out EPA-ng placement validation (computationally intensive)')
    vp.add_argument('--database',required=True);vp.add_argument('--output',required=True)
    vp.add_argument('--limit',type=int,default=20);vp.add_argument('--threads',type=int,default=4)
    vp.add_argument('--min-lwr',type=float,default=0.8)
    def validate_placement(a):
        from .phylo import placement_from_jplace, load_tree, filter_variable_sites
        from Bio import Phylo
        import tempfile
        if a.limit < 1: raise ValueError('--limit must be positive')
        db=load_database(a.database);out=Path(a.output);out.mkdir(parents=True,exist_ok=True)
        # TRUE leave-one-out: prune withheld tip and re-optimize branches for each fold.
        results=[]
        for i,(tip,seq) in enumerate(db['references'].items()):
            if i>=a.limit:break
            with tempfile.TemporaryDirectory() as tmp:
                d=Path(tmp);refs={k:v for k,v in db['references'].items() if k!=tip}
                refs, filtered_query, filter_stats = filter_variable_sites(refs, seq)
                with (d/'ref.fa').open('w') as f:
                    for name,s in refs.items():f.write(f'>{name}\n{s}\n')
                (d/'query.fa').write_text(f'>query\n{filtered_query}\n')
                tree=load_tree(Path(a.database)/'reference.nex');tree.prune(tip)
                for node in tree.find_clades():node.branch_length=None;node.comment=None
                Phylo.write(tree,str(d/'topology.nwk'),'newick')
                subprocess.run(['raxml-ng','--evaluate','--msa',str(d/'ref.fa'),'--tree',str(d/'topology.nwk'),
                    '--model','GTR+G+ASC_LEWIS','--prefix',str(d/'fit'),'--threads',str(a.threads)],check=True)
                (d/'epa').mkdir(parents=True, exist_ok=True)
                subprocess.run(['epa-ng','--ref-msa',str(d/'ref.fa'),'--tree',str(d/'fit.raxml.bestTree'),
                    '--query',str(d/'query.fa'),'--model','GTR+G','--threads',str(a.threads),
                    '--outdir',str(d/'epa')],check=True)
                result=placement_from_jplace(d/'epa'/'epa_result.jplace',
                    {k:v for k,v in db['clusters'].items() if k!=tip},a.min_lwr)
                weights=result.get('cluster_likelihood_weights', {})
                row={'isolate':tip,'true_cluster':db['clusters'][tip],
                    'tree_cluster':result['predicted_cluster'] or '',
                    'status':result['status'],'correct':result['predicted_cluster']==db['clusters'][tip],
                    **filter_stats}
                for cluster in ('1','2','3','4','5'):
                    row['score_cluster_'+cluster]=float(weights.get(cluster, 0.0))
                row['unresolved_weight']=float(result.get('unresolved_edge_weight', 0.0))
                results.append(row)
            print(f'Completed {i+1}/{min(a.limit,len(db["references"]))} folds',file=sys.stderr)
        with (out/'placement_leave_one_out.tsv').open('w',newline='') as fh:
            w=csv.DictWriter(fh,fieldnames=list(results[0]),delimiter='\t');w.writeheader();w.writerows(results)
        summary={'tested':len(results),'correct':sum(r['correct'] for r in results),
            'assigned':sum(r['status']=='assigned' for r in results),
            'note':'True tip-withheld placement validation on the selected SNP panel. Variable-site filtering is fitted separately within each training fold and applied to its query. Not independent FASTQ validation; EPA-ng model lacks the RAxML ascertainment correction.'}
        (out/'placement_validation.json').write_text(json.dumps(summary,indent=2)+'\n')
        print(json.dumps(summary,indent=2))
    vp.set_defaults(func=validate_placement)
    args=p.parse_args()
    try: args.func(args)
    except (ValueError,OSError,RuntimeError,subprocess.CalledProcessError) as exc:
        p.exit(2,f'ERROR: {exc}\n')
if __name__=='__main__':main()
