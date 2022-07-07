import unittest
from smn_analysis import compute_total_reads

class SMN_test(unittest.TestCase):
    def test_compute_total_reads(self):
        total_reads = compute_total_reads()