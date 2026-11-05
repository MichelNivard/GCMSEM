"""Checks for reference storage and independent GWAS matching."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import zipfile

import numpy as np
from scipy.io import loadmat

PATH = Path(__file__).resolve().parents[2]/'tutorials'/'prepare_original_sld4m.py'
spec = importlib.util.spec_from_file_location('original_sld4m', PATH)
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)


class OriginalSLD4MTest(unittest.TestCase):
    def test_matches_first_eligible_occurrence_in_reference_order(self):
        with tempfile.TemporaryDirectory() as temp:
            out = Path(temp)
            # Neither array is sorted. Duplicate GWAS records have distinct values.
            ref = np.array([50,20,70,10,30])
            ids = np.array([20,30,20,50,70,10,999])
            chi = np.array([2.,3.,200.,5.,7.,1.,9.])
            valid = np.array([True,True,True,True,False,True,True])
            module.match(ref,ids,chi,valid,out,'test')
            result = loadmat(out/'test.input.mat')
            np.testing.assert_array_equal(result['idx_ref'].ravel(),[1,2,4,5])
            np.testing.assert_array_equal(result['chisq'].ravel(),[5.,2.,1.,3.])
            self.assertEqual(json.loads((out/'test.input.json').read_text())['matched_snps'],4)

    def test_stored_zip_member_reads_and_seeks(self):
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp)/'example.zip'; payload = bytes(range(256))*3
            with zipfile.ZipFile(path,'w',compression=zipfile.ZIP_STORED) as z:
                z.writestr('unrelated','prefix'); z.writestr('reference.mat',payload)
            with zipfile.ZipFile(path) as z: info = z.getinfo('reference.mat')
            with module.StoredZipMember(path,info) as f:
                self.assertEqual(f.read(9),payload[:9]); f.seek(-5,2)
                self.assertEqual(f.read(),payload[-5:]); f.seek(77)
                b = bytearray(17); self.assertEqual(f.readinto(b),17)
                self.assertEqual(bytes(b),payload[77:94])


if __name__ == '__main__': unittest.main()
