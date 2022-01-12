from npysearch import read_fasta, write_fasta, blast, cigar_string, read_csv
from tempfile import NamedTemporaryFile as nt
import unittest, uuid, os


dna_example = {'DNASeq1': 'ATGCATCGGGCGAATT',
               'DNASeq2': 'TAGCTGGTGGACCACC'}

db_example = {'DBSeq1': 'TTGGATGCATCGGGCGAATTAACC',
              'DBSeq2': 'AACCTAGCTGGTGCACCACCGGTT'}

prot_example = {'ProtSeq1': 'LERAQC',
                'ProtSeq2': 'DYMFKW'}

line1 = ("Query1,Ref1,1,16,4,22,ATCGTGTACCAGGATG,ATCGTGTCCCACCAGGATG,"
         "19,16,0,3,0.842,7=3D9=\n")

line2 = ("Query1,Ref2,16,2,3,17,CATCCTGGTACACGA,CATCCTCGTACACGA,15,14"
         ",1,0,0.933,6=1X8=\n")

columns = ['QueryId', 'TargetId', 'QueryMatchStart', 'QueryMatchEnd',
           'TargetMatchStart', 'TargetMatchEnd', 'QueryMatchSeq',
           'TargetMatchSeq', 'NumColumns', 'NumMatches',
           'NumMismatches', 'NumGaps', 'Identity', 'Alignment']

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

    def test_read_csv(self):
        tmp_name = str(uuid.uuid4()) + ".csv"
        try:
            tmp_file = open(tmp_name, "w")
            tmp_file.write(line1)
            tmp_file.write(line2)
        finally:
            tmp_file.close()
            csv_dict = read_csv(tmp_name)
            os.remove(tmp_name)
        self.assertEqual(list(csv_dict.keys()), columns)
        self.assertEqual(csv_dict['QueryId'], ['Query1', 'Query1'])
        self.assertEqual(csv_dict['TargetId'], ['Ref1', 'Ref2'])

