import io
import sys

from Bio import SeqIO 

def main():
    genome_browser = {}
    for rec in SeqIO.parse('srsf6_sequence.fa', 'fasta'):
        genome_browser[rec.id] = rec.seq
    
    sanity_check = {}
    for rec in SeqIO.parse('./srsf6_sanity_sequence.fa', 'fasta'):
        sanity_check[rec.id] = rec.seq
    
    n_fail = 0
    n_total = 0
    for key in sanity_check:
        if str(genome_browser[key]) == str(sanity_check[key]):
            print key, ' match'
        else:
            print key, ' FAIL!'
            n_fail += 1
        n_total += 1

    print '********************************************************************************'
    print "Number of failures: ", n_fail
    print "Total: ", n_total

if __name__ == "__main__":
    main()
