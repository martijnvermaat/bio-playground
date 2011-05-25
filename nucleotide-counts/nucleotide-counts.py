#!/usr/bin/env python

# Count the number of occurrences per nucleotide in a fasta file.
#
# Usage:
#   ./nucleotide-counts.py sequence.fasta
#
# 2011-05-25, Martijn Vermaat <m.vermaat.hg@lumc.nl>


from __future__ import division
import sys
from Bio import SeqIO


def print_nucleotide_counts(sequence_file):
    """
    Count.
    """
    records = SeqIO.parse(open(sys.argv[1], 'r'), 'fasta')

    for record in records:
        # It doesn't feel right that we need this case hack...
        counts = dict([(n, record.seq.count(n) + record.seq.count(n.lower()))
                       for n in 'ATCGN'])
        total = len(record.seq)
        format = '%%s: %%%dd (%%6.3f%%%%)' % len(str(total))

        print record.name
        for n, count in counts.items():
            print format % (n, count, count / total * 100)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print """Count the number of occurrences per nucleotide.

Usage:
  {command} sequence.fasta""".format(command=sys.argv[0])
        sys.exit(1)
    print_nucleotide_counts(sys.argv[1])
