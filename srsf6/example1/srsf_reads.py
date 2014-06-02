import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sys.path.append('../../python')

from gen_data import *
from gtf_parser import *

def toSeqRecord(seq, name):
    return SeqRecord(Seq(seq), id = name, description = '')

def main():

    chr20 = SeqIO.read('../chr20_chunk.fa', 'fasta')
    chr20_seq = str(chr20.seq)

    trans = gtf_parse( '../srsf6_coords_fixed.gtf' )
    trans1 = trans['ENST00000244020']

    reads = []

    # generate some exonic reads
    reads.append(toSeqRecord(gen_read(chr20_seq, 104830, 50, trans = trans1), 'ex1.1'))
    reads.append(toSeqRecord(gen_read(chr20_seq, 104824, 50, False, trans = trans1), 'ex1.2'))
    reads.append(toSeqRecord(gen_read(chr20_seq, 104840, 50, trans = trans1), 'ex1.3'))
    reads.append(toSeqRecord(gen_read(chr20_seq, 104840, 50, trans = trans1), 'ex1.4'))
    reads.append(toSeqRecord(gen_read(chr20_seq, 104811, 50, False, trans = trans1), 'ex1.5'))

    # generate some reads at the exon-exon junction
    reads.append(toSeqRecord(gen_read(chr20_seq, 104970, 50, trans = trans1), 'ee1'))
    reads.append(toSeqRecord(gen_read(chr20_seq, 104967, 50, False, trans = trans1), 'ee2'))
    reads.append(toSeqRecord(gen_read(chr20_seq, 104975, 50, trans = trans1), 'ee3'))

    # gen some reads at the intron-exon junction
    reads.append(toSeqRecord(gen_read(chr20_seq, 104975, 50), 'ei1'))
    reads.append(toSeqRecord(gen_read(chr20_seq, 104977, 50), 'ei2'))
    reads.append(toSeqRecord(gen_read(chr20_seq, 104977, 50, False), 'ei3'))
    reads.append(toSeqRecord(gen_read(chr20_seq, 104976, 50, False), 'ei4'))

    # gen some intron reads
    reads.append(toSeqRecord(gen_read(chr20_seq, 105030, 50, ), 'in1'))
    reads.append(toSeqRecord(gen_read(chr20_seq, 105022, 50, False), 'in2'))
    reads.append(toSeqRecord(gen_read(chr20_seq, 105055, 50, False), 'in3'))
    reads.append(toSeqRecord(gen_read(chr20_seq, 105062, 50, False), 'in4'))
    reads.append(toSeqRecord(gen_read(chr20_seq, 105078, 50, False), 'in5'))

    # one read at other intron-exon junction
    reads.append(toSeqRecord(gen_read(chr20_seq, 105189, 50), 'ie1'))

    # generate some more exonic reads
    reads.append(toSeqRecord(gen_read(chr20_seq, 105213, 50, trans = trans1), 'ex2.1'))
    reads.append(toSeqRecord(gen_read(chr20_seq, 105210, 50, False, trans = trans1), 'ex2.2'))
    reads.append(toSeqRecord(gen_read(chr20_seq, 105220, 50, trans = trans1), 'ex2.3'))
    reads.append(toSeqRecord(gen_read(chr20_seq, 105207, 50, trans = trans1), 'ex2.4'))
    reads.append(toSeqRecord(gen_read(chr20_seq, 105242, 50, False, trans = trans1), 'ex2.5'))

    with open('reads.fa', 'w') as reads_handle:
        SeqIO.write(reads, reads_handle, 'fasta')

if __name__ == '__main__':
    main()
