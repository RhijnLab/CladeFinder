import json
from pathlib import Path
from cladefinder.phylo import _edge_clusters, placement_from_jplace, tree_report


def test_edge_annotation_and_purity():
    labels={'A':'1','B':'1','C':'2','D':'2'}
    tree='((A:0.1{1},B:0.1{2}):0.1{3},(C:0.1{4},D:0.1{5}):0.1{6}):0.1{7};'
    edges=_edge_clusters(tree,labels)
    assert edges[3]=='1'
    assert edges[6]=='2'


def test_jplace_weighted_assignment(tmp_path):
    labels={'A':'1','B':'1','C':'2','D':'2'}
    tree='((A:0.1{1},B:0.1{2}):0.1{3},(C:0.1{4},D:0.1{5}):0.1{6}):0.1{7};'
    data={'tree':tree,'fields':['edge_num','like_weight_ratio'],
          'placements':[{'n':['query'],'p':[[3,0.95],[6,0.05]]}]}
    f=tmp_path/'x.jplace';f.write_text(json.dumps(data))
    r=placement_from_jplace(f,labels)
    assert r['predicted_cluster']=='1'
    data['placements'][0]['p']=[[3,0.6],[6,0.4]]
    f.write_text(json.dumps(data))
    assert placement_from_jplace(f,labels)['status']=='ambiguous'
