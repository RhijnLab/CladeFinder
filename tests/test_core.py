import csv
from pathlib import Path
from cladefinder.core import genotype_base, nexus_tips, metadata, selected_sites, distances

def test_gt():
    assert genotype_base('1/1','A','G')=='G'
    assert genotype_base('0','A','G')=='A'
    assert genotype_base('0/1','A','G') is None
    assert genotype_base('./.','A','G') is None

def test_nexus_and_vcf(tmp_path):
    tree=tmp_path/'tree.nex'
    tree.write_text('#NEXUS\nbegin trees;\ntree TREE1 = [&R] (A_2015[&height=1]:1,(B_2016:1,C_2017:1):1);\nend;\n')
    assert set(nexus_tips(tree))=={'A_2015','B_2016','C_2017'}
    meta=tmp_path/'meta.csv';meta.write_text('ID,Cluster\nA,1\nB,2\nC,2\n')
    labels=metadata(meta)
    vcf=tmp_path/'tiny.vcf'
    vcf.write_text('##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n'+
                   'chr1\t3\t.\tA\tG\t.\tPASS\t.\tGT\t0\t1\t1\n'+
                   'chr1\t6\t.\tC\tT\t.\tPASS\t.\tGT\t0\t1\t1\n')
    sites,eligible=selected_sites(vcf,nexus_tips(tree),labels,10)
    assert eligible==2 and len(sites)==2
    refs={'A_2015':'AA','B_2016':'GG','C_2017':'GG'}
    result=distances('GG',refs,{'A_2015':'1','B_2016':'2','C_2017':'2'},min_compared=2)
    assert result['predicted_cluster']=='2'

def test_abstention():
    result=distances('NN',{'A':'AA'},{'A':'1'},min_compared=1)
    assert result['status']=='insufficient_data'
    result=distances('AA',{'A':'AA','B':'AA'},{'A':'1','B':'2'},min_compared=2)
    assert result['status']=='ambiguous'

def test_diagnostic_margin_and_ties():
    refs = {'A':'AAAA', 'B':'AAAT', 'C':'TTTT'}
    clusters = {'A':'1','B':'1','C':'2'}
    r = distances('AAAA', refs, clusters, min_compared=4, min_margin=.005)
    assert r['status'] == 'assigned' and r['nearest_cluster'] == '1'
    assert r['second_cluster'] == '2' and r['nearest']['compared_sites'] == 4
    tied = distances('AAAA', {'A':'AAAA','C':'AAAA'}, {'A':'1','C':'2'}, min_compared=4)
    assert tied['ambiguity_reason'] == 'equal_nearest_distance_across_clusters'
    assert tied['tied_clusters'] == ['1','2']
    assert tied['predicted_cluster'] is None

def test_margin_below_threshold():
    r = distances('AAAA', {'A':'AAAA','B':'AAAT'}, {'A':'1','B':'2'}, min_compared=4, min_margin=.3)
    assert r['status'] == 'ambiguous'
    assert r['nearest_cluster'] == '1'
    assert r['ambiguity_reason'] == 'margin_below_threshold'
