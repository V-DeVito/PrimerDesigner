import unittest

from fastapi.testclient import TestClient

from api.main import app
from primerdesignr.design import design_pcr_primers
from primerdesignr.thermo import analyze_pair, analyze_primer


TEMPLATE = (
    "TGACCTAAACGATGGTAGTAGGTTGGGAGCTTTCGAGAGGTCCGCCTTGGAAACGCGTTACTCGTGCGTGAATGTAGTGC"
    "AAGAGGAGGGTCAGTCGTGCTGAGACTGGGACTCTAAAATTACTGAAGCTCTTCCCATCCTCTCTAAGGTTTCGTGAGAC"
    "ACCACTTGGCGCTGGCCAGTACTGATCCCCTATTAGACCTATTGCCAAAACGTAGAAGGAACCTGCTCCAAAGCCCATGC"
    "CATGTTTGCGGATAGAACCATCGGTTAAGCGCCCTGGGTCATTGGTGAACTGGGTGAAAGATTCAATCTGCTCTGGTGCC"
    "CGGCCAGTCGCCATACGGAGCTCATTAGTATCGATTATAGAAGAGTTCTAGATCTTGCACTTGCGCGTCACAGAGTATAA"
    "TTACTCCGTTACGGCATGGCGATGAAGCGATACTATAGTGAAATGAAACTTGTGGCCATCCAGGTCGTTACCGGCCGTGG"
    "AACGGTGATATGTCAGACTTATAAGGTATATTGAGGGTTTAATAGCATTCCCGGCCAATTGTGTGTGCGATTGTTAAATG"
    "GCAGCTTCCGGACTCACTAGGAACTATTGAACCTTATTGTCGGTGGGGATGTCCAATTTTAGAATTGAGCGTCGATGCAA"
    "GGATCACAGTTTTAAAAGCAAGTTAAATCTAGGGTAAATAGCGGTGCTCTCATGCGTGTGGGTCGGGTGATGTTCAGAAA"
    "ATTTGCCTCCGAGAATACACTCATGGGTAGGACTCGCACTACCTAAGATTCCGCGCGCAGCACCGTTCGAGATTCTGCCC"
)

REPETITIVE_TEMPLATE = (
    "ATGCGTACGTAGCTAGCTACGATCGATCGTACGTAGCTAGCTAGCGATCGATCGTACGTAGCTAGCTAGC"
    "ATCGATCGATGCTAGCTAGCTAGCATCGATCGTACGTAGCTAGCTAGCATCGATCGATCGATCGATCGTAC"
    "GTAGCTAGCTAGCATCGATCGATCGTACGATCGATCGTAGCTAGCTAGCTAGCTAGCATCGATCGTACGAT"
    "CGATCGATCGATCGATCGTAGCTAGCTAGCATCGATCGATCGTAGCTAGCATCGATCGATCGTAGCTAGCT"
    "AGCATCGATCGATCGTAGCTAGCTAGCTAGCATCGATCG"
)


class CoreWorkflowTests(unittest.TestCase):
    def test_analyze_primer_reports_core_metrics(self):
        report = analyze_primer("GTCTTCACATCGGTTTGAAAGGAGG")
        self.assertEqual(report.length, 25)
        self.assertGreater(report.tm.tm, 50)
        self.assertLess(report.tm.tm, 70)
        self.assertIsInstance(report.warnings, list)

    def test_pair_analysis_uses_shared_conditions(self):
        low_salt = analyze_pair(
            "GTCTTCACATCGGTTTGAAAGGAGG",
            "AACCCGCTCCGATTAAAGCTACTTT",
            mv_conc=10,
            dv_conc=0,
            dna_conc=250,
        )
        high_salt = analyze_pair(
            "GTCTTCACATCGGTTTGAAAGGAGG",
            "AACCCGCTCCGATTAAAGCTACTTT",
            mv_conc=100,
            dv_conc=2,
            dna_conc=250,
        )
        self.assertNotEqual(low_salt.forward.tm.tm, high_salt.forward.tm.tm)

    def test_design_returns_ranked_candidates(self):
        result = design_pcr_primers(TEMPLATE, product_min=120, product_max=320, primer_count=3)
        self.assertEqual(result.template_length, len(TEMPLATE))
        self.assertEqual(len(result.candidates), 1)
        self.assertEqual(result.candidates[0].rank, 1)
        self.assertEqual(result.candidates[0].product_size, len(TEMPLATE))
        self.assertEqual(result.candidates[0].forward_coords.start, 0)
        self.assertEqual(result.candidates[0].reverse_coords.end, len(TEMPLATE) - 1)
        self.assertEqual(result.candidates[0].warnings, [])

    def test_amplicon_mode_can_move_inside_template(self):
        result = design_pcr_primers(
            TEMPLATE,
            product_min=120,
            product_max=320,
            primer_count=3,
            design_mode="amplicon",
        )
        self.assertGreaterEqual(result.candidates[0].product_size, 120)
        self.assertLessEqual(result.candidates[0].product_size, 320)
        self.assertGreater(result.candidates[0].forward_coords.start, 0)

    def test_amplicon_candidates_are_coordinate_diverse(self):
        result = design_pcr_primers(
            TEMPLATE,
            product_min=120,
            product_max=320,
            primer_count=8,
            design_mode="amplicon",
        )
        self.assertGreaterEqual(len(result.candidates), 3)
        for i, first in enumerate(result.candidates):
            for second in result.candidates[i + 1:]:
                same_local_cluster = (
                    abs(first.forward_coords.start - second.forward_coords.start) <= 8
                    and abs(first.forward_coords.end - second.forward_coords.end) <= 8
                    and abs(first.reverse_coords.start - second.reverse_coords.start) <= 8
                    and abs(first.reverse_coords.end - second.reverse_coords.end) <= 8
                    and abs(first.product_size - second.product_size) <= 16
                )
                self.assertFalse(same_local_cluster)

    def test_design_warns_when_no_clean_pair_is_found(self):
        result = design_pcr_primers(
            "ATGC" * 160,
            product_min=120,
            product_max=320,
            primer_count=3,
            design_mode="amplicon",
        )
        self.assertTrue(result.warnings)
        self.assertIn("Relax product size", result.warnings[0])


class ApiWorkflowTests(unittest.TestCase):
    def setUp(self):
        self.client = TestClient(app)

    def test_health(self):
        res = self.client.get("/health")
        self.assertEqual(res.status_code, 200)
        self.assertEqual(res.json()["status"], "ok")

    def test_design_endpoint(self):
        res = self.client.post(
            "/design",
            json={
                "template": TEMPLATE,
                "product_min": 120,
                "product_max": 320,
                "primer_count": 2,
            },
        )
        self.assertEqual(res.status_code, 200)
        body = res.json()
        self.assertGreaterEqual(len(body["candidates"]), 1)
        self.assertIn("explanations", body["candidates"][0])

    def test_rejects_invalid_template(self):
        res = self.client.post(
            "/design",
            json={"template": "ATGCNNNNATGCATGCATGC", "product_min": 40, "product_max": 80},
        )
        self.assertEqual(res.status_code, 400)

    def test_golden_gate_endpoint(self):
        res = self.client.post(
            "/golden-gate",
            json={"overhangs": ["AACG", "AATG", "ATAG"], "enzyme": "BsaI"},
        )
        self.assertEqual(res.status_code, 200)
        self.assertTrue(res.json()["all_unique"])


if __name__ == "__main__":
    unittest.main()
