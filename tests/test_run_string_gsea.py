import unittest
from a373_string_gsea.run_string_gsea import get_ox_fields

class TestRunStringGSEA(unittest.TestCase):

    def test_get_ox_fields(self):
        fasta_content = (">sp|P0DTC2|SPIKE_SARS2 Spike glycoprotein OS=Severe acute respiratory syndrome coronavirus 2 OX=2697049 GN=S"
        "\n"
        ">sp|P0DTC1|NCAP_SARS2 Nucleoprotein OS=Severe acute respiratory syndrome coronavirus 2 OX=2697049 GN=N")
        expected_ox = ['2697049', '2697049']
        result = get_ox_fields(fasta_content.splitlines())
        self.assertEqual(result, expected_ox)


if __name__ == '__main__':
    unittest.main()