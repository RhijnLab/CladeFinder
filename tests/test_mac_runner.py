import importlib.machinery
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

RUNNER = Path(__file__).resolve().parents[1] / 'cladefinder-mac'
loader = importlib.machinery.SourceFileLoader('macrunner', str(RUNNER))
spec = importlib.util.spec_from_loader(loader.name, loader)
module = importlib.util.module_from_spec(spec)
loader.exec_module(module)

class PublishTests(unittest.TestCase):
    def test_publishes_complete_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            stage = root/'stage'; dest = root/'db'
            stage.mkdir(); dest.mkdir()
            for name in module.ARTIFACTS:
                (stage/name).write_text(json.dumps({'tree':'placement_tree.nwk'}) if name.endswith('.json') else 'test')
            module.publish_artifacts(stage, dest)
            self.assertEqual((dest/'placement_tree.nwk').read_text(), 'test')

    def test_missing_artifact_preserves_old_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); stage=root/'stage'; dest=root/'db'
            stage.mkdir(); dest.mkdir()
            (dest/'placement_tree.nwk').write_text('old')
            (stage/'placement_manifest.json').write_text('{"tree":"placement_tree.nwk"}')
            with self.assertRaises(RuntimeError):
                module.publish_artifacts(stage, dest)
            self.assertEqual((dest/'placement_tree.nwk').read_text(), 'old')

if __name__ == '__main__': unittest.main()
