import io
import random
import sys

from gtf_parser import *

def gen_contig(contig_len):
    alphabet = ['A', 'C', 'G', 'T']

    return ''.join([alphabet[random.randint(0, 3)] for i in xrange(contig_len)])

def simple_annotation():
    # trans1 = transcript('trans1', )
    pass

def complement(base):
    base = base.upper()
    def complement_base(x):
        if x == 'A':
            return 'T'
        elif x == 'T':
            return 'A'
        elif x == 'G':
            return 'C'
        return 'G'
    return ''.join([complement_base(b) for b in base])

def reverse_complement(seq):
    return complement(seq[::-1])

# given a start site (genome 1 based coords), a contig, 
# and possibly a transcript will generate a read
def gen_read(contig, start_pos, read_len, forward_strand = True, 
        paired = False, insert_len = False, trans = None):

    def continuous_read(start, stop):
        return contig[start:stop]
            

    left_end_pos = start_pos + read_len - 1
    left = None
    if  left_end_pos > len(contig):
        print >> sys.stderr, 'Read is too long to start at position ', start_pos, 
        ' (contig length: ', len(contig), ')'
        raise Exception('Invalid start position (start_pos + read_len - 1 > len(contig)')

    if trans is None:
        left = continuous_read(start_pos - 1, left_end_pos)
    else:
        left_start_compat = trans.compatible(start_pos)
        left_end_compat = trans.compatible(left_end_pos)

        if left_start_compat == -1:
            raise Exception('Start site does not belong to transcript')

        if left_start_compat != left_end_compat:
            l1 = continuous_read(start_pos - 1, trans.exons[left_start_compat][1])
            remaining_len = read_len - (trans.exons[left_start_compat][1] - start_pos + 1)
            l2 = continuous_read(trans.exons[left_start_compat + 1][0] - 1, 
                    trans.exons[left_start_compat + 1][0] + remaining_len - 1)
            left = l1 + l2
            # TODO: generate spliced read
        else:
            left = continuous_read(start_pos - 1, left_end_pos)

    if not paired:
        if not forward_strand:
            left = reverse_complement(left)
        return left

    # TODO: deal w/ paired end
    
    return None 
