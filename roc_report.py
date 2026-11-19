#!/usr/bin/env python3
"""Make one-vs-rest ROC SVG and AUROC JSON from withheld-isolate scores. No dependencies."""
import argparse, csv, json
from pathlib import Path

def roc(rows, cluster):
    pairs = sorted(((float(r['score_cluster_'+cluster]), int(r['true_cluster']==cluster)) for r in rows), reverse=True)
    positives=sum(y for _,y in pairs); negatives=len(pairs)-positives
    if not positives or not negatives: return None
    points=[(0.,0.)]; tp=fp=0; i=0
    while i<len(pairs):
        score=pairs[i][0]
        while i<len(pairs) and pairs[i][0]==score:
            if pairs[i][1]: tp+=1
            else: fp+=1
            i+=1
        points.append((fp/negatives,tp/positives))
    auc=sum((x2-x1)*(y1+y2)/2 for (x1,y1),(x2,y2) in zip(points,points[1:]))
    return {'auc':auc,'positives':positives,'negatives':negatives,'points':points}

def main():
    ap=argparse.ArgumentParser();ap.add_argument('tsv');ap.add_argument('--output',default='roc-results');args=ap.parse_args()
    with open(args.tsv,newline='') as f: rows=list(csv.DictReader(f,delimiter='\t'))
    if not rows: raise SystemExit('No validation rows')
    for c in '12345':
        if 'score_cluster_'+c not in rows[0]: raise SystemExit('Missing score_cluster_'+c+'; rebuild CladeFinder and rerun validation')
    out=Path(args.output);out.mkdir(parents=True,exist_ok=True)
    result={c:roc(rows,c) for c in '12345'}
    (out/'roc_summary.json').write_text(json.dumps({'tested':len(rows),'clusters':result},indent=2)+'\n')
    colors=['#1976d2','#e67e22','#2e7d32','#9c27b0','#d32f2f']
    svg=['<svg xmlns="http://www.w3.org/2000/svg" width="760" height="600" viewBox="0 0 760 600">','<rect width="760" height="600" fill="white"/>','<text x="75" y="30" font-family="Arial" font-size="20">Leave-one-out placement ROC (one-vs-rest)</text>','<path d="M75 60 V520 H535" stroke="black" stroke-width="2" fill="none"/>','<path d="M75 520 L535 60" stroke="#999" stroke-dasharray="6 6" fill="none"/>','<text x="240" y="560" font-family="Arial">False positive rate</text>','<text x="15" y="300" transform="rotate(-90 15 300)" font-family="Arial">True positive rate</text>']
    for j,c in enumerate('12345'):
        r=result[c]
        if r is None: continue
        pts=' '.join(f'{75+460*x:.2f},{520-460*y:.2f}' for x,y in r['points'])
        svg.append(f'<polyline points="{pts}" stroke="{colors[j]}" stroke-width="2.5" fill="none"/>')
        svg.append(f'<text x="555" y="{90+j*30}" font-family="Arial" font-size="14" fill="{colors[j]}">Cluster {c}: AUC {r["auc"]:.3f}</text>')
    svg.append('</svg>');(out/'roc_curves.svg').write_text('\n'.join(svg))
    print('Saved',out/'roc_curves.svg','and',out/'roc_summary.json')
if __name__=='__main__':main()
