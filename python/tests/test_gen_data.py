import unittest

from Bio import SeqIO

from gen_data import *

class TestGenData(unittest.TestCase):
    def setUp(self):
        pass

    def teatDown(self):
        pass

    def test_gen_contig(self):
        contig = gen_contig(30)
        self.assertEqual( len(contig), 30 )
        self.assertRegexpMatches( contig, '[ACGT]+' )

    def test_complement(self):
        self.assertEqual( complement('ACCGT'), 'TGGCA' )

    def test_reverse_complement(self):
        self.assertEqual( reverse_complement('GACTAA'), 'TTAGTC' )


    def test_gen_read(self):
        chr20 = SeqIO.read('../srsf6/chr20_chunk.fa', 'fasta')
        chr20_seq = str(chr20.seq)

        # basic SE
        se_basic = gen_read(chr20_seq, 1, 30)
        self.assertEqual( se_basic,  'TTGGCCCTTCACATAGGATGCTGTGTTTAA' )

        se_rev = gen_read(chr20_seq, 1, 30, False)
        # print reverse_complement('TTGGCCCTTCACATAGGATGCTGTGTTTAA')
        self.assertEqual( se_rev, 'TTAAACACAGCATCCTATGTGAAGGGCCAA')

        all_trans = gtf_parse( '../srsf6/srsf6_coords_fixed.gtf' )

        # basic spliced SE
        se_spliced = gen_read(chr20_seq, 104964, 28, trans = all_trans['ENST00000244020'])
        self.assertEqual(se_spliced, 'AGTAGACCTCAAAAATGG' + 'GTACGGCTTC')

        # TODO: generate a paired-end read
        # TODO: generate a spliced paired-end read from trans
        pass
