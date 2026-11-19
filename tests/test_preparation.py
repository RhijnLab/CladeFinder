import json
from pathlib import Path
import pytest
from cladefinder import phylo


def fixture_db(tmp_path):
    db = tmp_path / 'db'
    db.mkdir()
    (db / 'database.json').write_text(json.dumps({
        'references': {'A': 'ACGT', 'B': 'ACGA', 'C': 'TCGA'},
        'clusters': {'A': '1', 'B': '1', 'C': '2'}, 'site_count': 4
    }))
    (db / 'reference.nex').write_text('#NEXUS\nBegin trees;\nTree t = ((A:1,B:1):1,C:1);\nEnd;\n')
    return db


def fake_raxml(cmd, check):
    prefix = Path(cmd[cmd.index('--prefix') + 1])
    Path(str(prefix) + '.raxml.bestTree').write_text('((A:0.1,B:0.2):0.1,C:0.2);\n')
    Path(str(prefix) + '.raxml.bestModel').write_text('GTR+G\n')


def test_prepare_twice_and_nonempty(tmp_path, monkeypatch):
    db = fixture_db(tmp_path)
    monkeypatch.setattr(phylo.subprocess, 'run', fake_raxml)
    for _ in range(2):
        phylo.prepare_placement(db, threads=2)
        for name in ('reference_snps.fasta', 'nexus_topology.nwk',
                     'placement_tree.nwk', 'placement_model.txt', 'placement_manifest.json'):
            assert (db / name).stat().st_size > 0


def test_raxml_failure_preserves_outputs(tmp_path, monkeypatch):
    db = fixture_db(tmp_path)
    monkeypatch.setattr(phylo.subprocess, 'run', fake_raxml)
    phylo.prepare_placement(db)
    before = {name: (db / name).read_bytes() for name in ('reference_snps.fasta', 'placement_tree.nwk', 'placement_manifest.json')}
    def broken(cmd, check):
        raise RuntimeError('simulated RAxML failure')
    monkeypatch.setattr(phylo.subprocess, 'run', broken)
    with pytest.raises(RuntimeError, match='simulated'):
        phylo.prepare_placement(db)
    assert before == {name: (db / name).read_bytes() for name in before}
