import io
import sys

from gtf_parser import *


def main():
    contig_start = 41981800
    srsf6_gtf = gtf_parse('./srsf6_ensembl_exon.gtf')
    for key in srsf6_gtf:
        cur_trans = srsf6_gtf[key]
        for i in xrange(len(cur_trans.exons)):
            cur_trans.exons[i] = (cur_trans.exons[i][0] - contig_start + 1,
                    cur_trans.exons[i][1] - contig_start + 1)

    gtf_write(srsf6_gtf, sys.stdout)

if __name__ == '__main__':
    main()
