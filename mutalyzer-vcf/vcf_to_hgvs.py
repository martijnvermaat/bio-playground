#!/usr/bin/env python

# Create HGVS descriptions from a VCF file.
#
# Usage:
#   ./vcf_to_hgvs.py sample.vcf
#
# The VCF file must be in VCFv4.1 format (as created by Samtools 0.1.16
# for example).
#
# Copyright (c) 2011 Leiden University Medical Center <humgen@lumc.nl>
# Copyright (c) 2011 Martijn Vermaat <m.vermaat.hg@lumc.nl>


import sys


def main(vcf_file):
    """
    Read lines from VCF file and print HGVS descriptions.
    """
    try:
        vcf = open(vcf_file, 'r')
    except IOError as (_, message):
        print 'Could not read VCF file: %s' % vcf_file
        sys.exit(1)

    while True:
        line = vcf.readline()
        if not line:
            break

        if line.startswith('#'):
            continue

        parts = line.split()

        if len(parts) < 5:
            print 'Could not interpret line: %s' % line
            sys.exit(1)

        try:
            chromosome = 'chr%d' % int(parts[0])
        except ValueError:
            chromosome = parts[0]

        try:
            position = int(parts[1])
        except:
            print 'Could not read position: %s' % line

        reference = parts[3].upper()
        alternates = parts[4].upper()

        for alternate in alternates.split(','):

            if len(reference) == len(alternate) == 1:
                # SNP
                print '%s:g.%d%s>%s' % \
                      (chromosome, position, reference, alternate)
            else:
                # Consider this an indel
                print '%s:g.%d_%ddel%sins%s' % \
                      (chromosome, position, position + len(reference) - 1,
                       reference, alternate)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print """Create HGVS descriptions from a VCF file.

Usage:
  {command} sample.vcf

The VCF file must be in VCFv4.1 format (as created by Samtools 0.1.16
for example).""".format(command=sys.argv[0])
        sys.exit(1)
    main(sys.argv[1])
