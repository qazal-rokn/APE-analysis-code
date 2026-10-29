import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from tools.ape_unified_orchestrator import STAGES, execute

class TestUnifiedOrchestrator(unittest.TestCase):
    def test_preflight_and_reuse_without_mutation(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d)
            p=root/"immutable";p.write_bytes(b"frozen")
            h=hashlib.sha256(p.read_bytes()).hexdigest()
            config=root/"plan.json"
            config.write_text(json.dumps({"schema_version":1,"run_dir":str(root/"run"),
                "stages":{stage:{"mode":"reuse","outputs":[{"path":str(p),"sha256":h}]} for stage in STAGES}}))
            self.assertEqual(execute(config,preflight=True)["status"],"PREFLIGHT_OK_NO_EXECUTION")
            self.assertEqual(execute(config)["status"],"STAGES_COMPLETED_PROVISIONAL")
            self.assertEqual(p.read_bytes(),b"frozen")
            with self.assertRaises(FileExistsError):
                execute(config)
    def test_bad_hash_fails(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);p=root/"source";p.write_text("unchanged")
            config=root/"plan.json"
            config.write_text(json.dumps({"schema_version":1,"run_dir":str(root/"run"),
                "stages":{stage:{"mode":"reuse","outputs":[{"path":str(p),"sha256":"0"*64}]} for stage in STAGES}}))
            with self.assertRaises(RuntimeError):
                execute(config,preflight=True)
if __name__ == "__main__": unittest.main()
