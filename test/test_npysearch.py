from npysearch import read_fasta, write_fasta, blast, cigar_string
import tempfile

def test_read_fasta() -> None:
    dna_iter = read_fasta("data/query.fasta")
    assert list(dna_iter.keys())[0] == "MISEQ:1:1101:13226:2432#TTGGTCTG/1_1"
    assert list(dna_iter.values())[0] == "GTGAGTGATGGTTGAGGTA"
    prot_iter = read_fasta("data/prot.fasta")
    assert list(prot_iter.keys())[0] == "Seq_1"
    assert list(prot_iter.values())[0] == "LERAQCGHDYMFKWTSIVPN"

def test_write_fasta() -> None:
    fasta_dict = {'Seq1': 'ATGC'}
    with tempfile.NamedTemporaryFile(suffix='.fasta') as t:
        write_fasta(t.name, fasta_dict)
        dna_iter = read_fasta(t.name)
        assert list(dna_iter.keys())[0] == "Seq1"
        assert list(dna_iter.values())[0] == "ATGC"

def test_blast() -> None:
    blast_res = blast("data/query.fasta", "data/db.fasta")
    assert list(blast_res.keys()) == ['QueryId', 'TargetId', 'QueryMatchStart',
            'QueryMatchEnd', 'TargetMatchStart', 'TargetMatchEnd',
            'QueryMatchSeq', 'TargetMatchSeq', 'NumColumns',
            'NumMatches', 'NumMismatches', 'NumGaps', 'Identity', 'Alignment']
    assert blast_res['QueryId'][0] == "MISEQ:1:1101:13226:2432#TTGGTCTG/1_1"
    assert blast_res['QueryMatchSeq'][0] == "GTGAGTGATGGTTGAGGTA"
    assert blast_res['TargetMatchSeq'][0] == "GTGAGTGATGGTTGAGGTA"

def test_cigar_string() -> None:
    assert cigar_string("ATGC", "TTGC") == "1X3="

