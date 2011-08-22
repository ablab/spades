import gzip
import os
import sys
from unittest import TestCase
from sequenceDistanceMeter import SequenceDistanceMeter

PATH_TO_BOWTIE = "/home/arts/Documents/bio/bowtie-0.12.7"
PATH_TO_INPUT = "/home/arts/Documents/input"
PATH_TO_OUTPUT = "/home/arts/Documents/output"
PATH_TO_TEST_INPUT = "/home/arts/Documents/test_input"
PATH_TO_TEST_INPUT_COMPRESSED = PATH_TO_TEST_INPUT + ".gz"
INDEX = "MG1655-K12"

class SequenceDistanceMeterTest(TestCase):
    def setUp(self):
        self.sequence_distance_meter = SequenceDistanceMeter(PATH_TO_BOWTIE, PATH_TO_INPUT, PATH_TO_OUTPUT, INDEX)

    def test_join_with_delimiter(self):
        self.assertEqual(self.sequence_distance_meter.join_with_delimiter(PATH_TO_BOWTIE.split("/")), PATH_TO_BOWTIE + "/")

    def test_unify_path(self):
        self.assertEqual(self.sequence_distance_meter.unify_path(PATH_TO_BOWTIE), PATH_TO_BOWTIE + "/")
        self.assertEqual(self.sequence_distance_meter.unify_path(PATH_TO_BOWTIE + "/"), PATH_TO_BOWTIE + "/")
        self.assertEqual(self.sequence_distance_meter.unify_path(PATH_TO_BOWTIE + "/bowtie"), PATH_TO_BOWTIE + "/")
        self.assertEqual(self.sequence_distance_meter.unify_path(PATH_TO_BOWTIE + "some junk"), None)

    def test_initialize_path_to_bowtie(self):
        self.sequence_distance_meter.initialize_path_to_bowtie(PATH_TO_BOWTIE)
        self.assertEqual(self.sequence_distance_meter.path_to_bowtie, PATH_TO_BOWTIE + "/")
        self.assertRaises(IOError, self.sequence_distance_meter.initialize_path_to_bowtie, "")
        self.assertRaises(IOError, self.sequence_distance_meter.initialize_path_to_bowtie, "junk")

    def test_initialize_input(self):
        # test common input
        try:
            input_to_write = open(PATH_TO_TEST_INPUT, "w")
            lines_to_write = [
                    "ATTTCACTTTTACCCACAGTTCAAGGTGAACAGGCGCT CGCGTAACAAAAGTGTCTATAATCACGGCAGAAAAGTCCACATTGATTATTTGCACGGCG\n"
                    ,
                    "GAAGGCTTCTTTGATATAGACGTTACGCACGGTACAGCCGTTGATCTGGCTATGGCCCTGCGTG TGATGATCGGCATTAATATGGCCGAAAAAGACGA\n"
                    ,
                    "TAATGCGGGCGATGATTTCCGCACTTGTATGGAATCTTTAATTA ACGCCTTCGTAGCCTTCTGACGGCATTAACACTTTCGAGCGGCCAGGTAGCGTG"]
            input_to_write.writelines(lines_to_write)
            input_to_write.close()
        except IOError:
            sys.exit("Cannot create test input file. Terminated")
        self.assertEqual(lines_to_write, self.sequence_distance_meter.initialize_input(PATH_TO_TEST_INPUT).readlines())
        os.remove(PATH_TO_TEST_INPUT)

        # test gzip input
        try:
            input_to_write = gzip.open(PATH_TO_TEST_INPUT_COMPRESSED, "w")
            lines_to_write = [
                    "ATTTCACTTTTACCCACAGTTCAAGGTGAACAGGCGCT CGCGTAACAAAAGTGTCTATAATCACGGCAGAAAAGTCCACATTGATTATTTGCACGGCG\n"
                    ,
                    "GAAGGCTTCTTTGATATAGACGTTACGCACGGTACAGCCGTTGATCTGGCTATGGCCCTGCGTG TGATGATCGGCATTAATATGGCCGAAAAAGACGA\n"
                    ,
                    "TAATGCGGGCGATGATTTCCGCACTTGTATGGAATCTTTAATTA ACGCCTTCGTAGCCTTCTGACGGCATTAACACTTTCGAGCGGCCAGGTAGCGTG"]
            input_to_write.writelines(lines_to_write)
            input_to_write.close()
        except IOError:
            sys.exit("Cannot create test input file. Terminated")
        self.assertEqual(lines_to_write, self.sequence_distance_meter.initialize_input(PATH_TO_TEST_INPUT_COMPRESSED).readlines())
        os.remove(PATH_TO_TEST_INPUT_COMPRESSED)

    def test_align_sequence(self):
        self.sequence_distance_meter.safe_initialize_path_to_bowtie(PATH_TO_BOWTIE)
        sequence = "GAAGGCTTCTTTGATAT"
        self.assertEqual(self.sequence_distance_meter.align_sequence(sequence), [3887212, 1016059, 4059797])

    def test_compute_distances(self):
        self.assertEqual(self.sequence_distance_meter.compute_distances([2, 5], [3, 4]), [1, 1, 2, 2])

    def test_process_all_pairs(self):
        self.sequence_distance_meter.process_all_pairs()


        




        
