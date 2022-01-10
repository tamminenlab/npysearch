from npysearch import read_fasta, write_fasta, blast, cigar_string
from tempfile import NamedTemporaryFile as nt
import unittest


dna_example = {'DNASeq1': 'ATGCATCGGGCGAATT',
               'DNASeq2': 'TAGCTGGTGGACCACC'}

db_example = {'DBSeq1': 'TTGGATGCATCGGGCGAATTAACC',
              'DBSeq2': 'AACCTAGCTGGTGCACCACCGGTT'}

prot_example = {'ProtSeq1': 'LERAQC',
                'ProtSeq2': 'DYMFKW'}

class Tests(unittest.TestCase):

    def test_fasta_io(self):
        with nt(suffix='.fasta') as dna_fasta:
            write_fasta(dna_fasta.name, dna_example)
            dna_dict = read_fasta(dna_fasta.name)
            self.assertEqual(list(dna_dict.keys())[0], "DNASeq1")
            self.assertEqual(list(dna_dict.values())[0], "ATGCATCGGGCGAATT")
            self.assertEqual(list(dna_dict.keys())[1], "DNASeq2")
            self.assertEqual(list(dna_dict.values())[1], "TAGCTGGTGGACCACC")
        with nt(suffix='.fasta') as prot_fasta:
            write_fasta(prot_fasta.name, prot_example)
            prot_dict = read_fasta(prot_fasta.name)
            self.assertEqual(list(prot_dict.keys())[0], "ProtSeq1")
            self.assertEqual(list(prot_dict.values())[0], "LERAQC")
            self.assertEqual(list(prot_dict.keys())[1], "ProtSeq2")
            self.assertEqual(list(prot_dict.values())[1], "DYMFKW")

    def test_blast(self):
        with nt(suffix='.fasta') as dna_fasta, nt(suffix='.fasta') as db_fasta:
            write_fasta(dna_fasta.name, dna_example)
            dna_dict = read_fasta(dna_fasta.name)
            write_fasta(db_fasta.name, db_example)
            db_dict = read_fasta(db_fasta.name)
            blast_res = blast(dna_dict, db_dict)
            self.assertEqual(blast_res['QueryId'][0], 'DNASeq1')
            self.assertEqual(blast_res['TargetId'][0], 'DBSeq1')
            self.assertEqual(blast_res['QueryMatchSeq'][0],  'ATGCATCGGGCGAATT')
            self.assertEqual(blast_res['TargetMatchSeq'][0], 'ATGCATCGGGCGAATT')
            self.assertEqual(blast_res['QueryMatchSeq'][1],  'TAGCTGGTGGACCACC')
            self.assertEqual(blast_res['TargetMatchSeq'][1], 'TAGCTGGTGCACCACC')

    def test_cigar_string(self):
        self.assertEqual(cigar_string("ATGC", "TTGC"), "1X3=")
        with self.assertRaises(AssertionError):
            cigar_string("AATGC", "TTGC")


