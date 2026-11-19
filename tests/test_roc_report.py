import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from roc_report import roc

def test_roc_perfect():
    rows=[{'true_cluster':'1','score_cluster_1':'0.9'}, {'true_cluster':'2','score_cluster_1':'0.1'}]
    assert roc(rows,'1')['auc']==1.0

def test_roc_tied():
    rows=[{'true_cluster':'1','score_cluster_1':'0.5'}, {'true_cluster':'2','score_cluster_1':'0.5'}]
    assert roc(rows,'1')['auc']==0.5

def test_roc_no_negative():
    assert roc([{'true_cluster':'1','score_cluster_1':'0.8'}],'1') is None
