#!/usr/bin/env python

# Annotate variants using SeattleSeq Annotation.
#
# Given a VCF file with variants, submit the file to the SeattleSeq Annotation
# web interface [1]. The annotation result is retrieved as plain-text
# tab-separated file.
#
# Usage:
#   ./seattle_seq_annotation.py snp|indel <input.vcf> <mail@domain.com>
#
# The result is written to disk as <input.vcf.annotation>. Temporary data is
# written to <input.vcf.part.*> files and removed afterwards.
#
# Requires the poster Python library [2].
#
# [1] http://snp.gs.washington.edu/SeattleSeqAnnotation131/
# [2] http://atlee.ca/software/poster/
#
# 2011-06-14, Martijn Vermaat <m.vermaat.hg@lumc.nl>


# Seconds to wait between polling for job result
WAIT_STEP = 120

# Maximum seconds to wait for job result
WAIT_MAX = 43200

# Print debugging information
DEBUG = True

# Columns to include
COLUMNS = ['sampleAlleles',
#           'dbSNPGenotype',
           'allelesDBSNP',
           'scorePhastCons',
           'consScoreGERP',
           'chimpAllele',
           'CNV',
           'geneList',
           'HapMapFreq',
           'hasGenotypes',
           'dbSNPValidation',
           'repeats',
           'proteinSequence',
#           'cDNAPosition',
           'polyPhen',
           'clinicalAssociation'
#           'distanceToSplice',
#           'microRNAs',
#           'grantham'
           ]

# Maximum number of variants SeattleSeq accepts
MAX_VARIANTS = 900000   # To be sure, actually 1000000

# SeattleSeq Annotation location
BASE_URL = 'http://snp.gs.washington.edu/SeattleSeqAnnotation131/'
POST_URL = BASE_URL + 'BatchQueryServlet'


import sys
import os
import time
from collections import defaultdict
import uuid
from poster.encode import multipart_encode
from poster.streaminghttp import register_openers
import urllib
import urllib2


def seattle_seq_annotation(mode, vcf_file, address):
    """
    Submit a VCF file to the SeattleSeq Annotation web interface. The
    annotation result is retrieved as plain-text tab-separated and written
    to disk.
    """
    part_files = create_parts(vcf_file)

    submissions = [submit_part(mode, p, address) for p in part_files]

    waiting = 0
    first = True
    versions = set()
    summary = defaultdict(int)

    # Todo: perhaps first truncate the annotation file

    for monitor_url, result_url in submissions:
        waiting += wait_for_result(monitor_url, waiting)
        append_result(result_url, vcf_file + '.annotation', versions,
                      summary, discard_header=not first)
        first = False

    append_summary(vcf_file + '.annotation', versions, summary)


def submit_part(mode, part_file, address):
    """
    Submit a VCF file to the SeattleSeq Annotation server and return as a
    tuple the url to monitor the job and the url to get the result.

    The {mode} argument can be 'snp' (default) or 'indel'.
    """
    format = 'VCFIndel' if mode == 'indel' else 'VCF'

    try:
        part = open(part_file, 'r')
    except IOError as (_, message):
        fatal_error(message)

    parameters = [('genotypeSource',   'FileInput'),
                  ('EMail',            address),
                  ('GenotypeFile',     part),
                  ('fileFormat',       format),
                  ('outputFileFormat', 'original'),
                  ('geneData',         'NCBI'),
                  ('HapMapFreqType',   'HapMapFreqMinor'),
                  ('gFetch',           'Display Genotypes')]
    for column in COLUMNS:
        parameters.append( ('columns', column) )

    debug('Submitting for %s annotation: %s' % (mode, part_file))

    # Response contains result url and monitor url separated by a comma
    response = post_multipart(POST_URL, parameters)
    urls = response.readline().strip().split(',')
    response.close()

    part.close()
    os.unlink(part_file)

    if not len(urls) == 2:
        fatal_error('Could not read urls from submit response.')

    result_url, monitor_url = urls

    debug('Monitor url: %s' % monitor_url)
    debug('Result url: %s' % result_url)

    return monitor_url, result_url


def create_parts(vcf_file):
    """
    Split the VCF file in parts that contain less variants than the maximum
    number allowed by SeattleSeq.

    Make sure to include the header lines in each part and prepend a line
    needed by the SeattleSeq Annotation server.

    Return the filenames of the generated parts.
    """
    try:
        vcf = open(vcf_file, 'r')
    except IOError as (_, message):
        fatal_error(message)

    # Header line for automated processing
    header = '# autoFile vcfAuto.txt\n'

    # Filenames for parts
    part_files = []

    # Open a new part file
    def open_part():
        # Nicer would be to use something like mktemp()
        part_file = vcf_file + '.part.' + uuid.uuid1().hex
        try:
            part = open(part_file, 'w')
        except IOError as (_, message):
            fatal_error(message)
        part_files.append(part_file)
        part.write(header)
        return part

    part = None

    # Start with MAX_VARIANTS so we first open a new part file
    line_count = MAX_VARIANTS

    # Buffered read/write
    while True:
        line = vcf.readline()
        if not line:
            break

        # Read the original header lines
        if line.startswith('#'):
            header += line
            continue

        line_count += 1
        if line_count > MAX_VARIANTS:
            if part:
                part.close()
            part = open_part()
            line_count = 0

        part.write(line)

    # Call for help in case of an empty VCF file (SeattleSeq does not handle
    # this gracefully)
    if line_count == 0:
        fatal_error('VCF file contains no variants.')

    debug('Splitted VCF file to temporary files.')

    vcf.close()
    return part_files


def wait_for_result(monitor_url, waiting):
    """
    Monitor the job for progress. Return the number of seconds it took when
    the job is completed.
    """
    while True:
        if waiting > WAIT_MAX:
            fatal_error('Job took too long.')

        time.sleep(WAIT_STEP)
        waiting += WAIT_STEP

        # Response contains number of variants processed and total number of
        # variants separated by a comma
        response = get(monitor_url)
        monitor = response.readline().strip().split(',')
        response.close()

        if not len(monitor) == 2:
            fatal_error('Could not read progress from monitor.')

        processed, total = monitor

        # Todo: If we submit a file with 0 variants, the monitor page keeps
        # saying '0 variations in your file have been processed', even after
        # the job has finished (and an email has been sent).
        # So from the monitor we cannot detect completion in this case. The
        # only solution I can think of is checking the input file for this
        # special case.

        # If no variants have yet been processed, we always get 0,0
        if total == '0':
            debug('Waiting for part: 0%')
            continue

        debug('Waiting for part: %s / %s' % (processed, total))

        # See if we are done
        if processed == total:
            break

    return waiting


def append_result(result_url, output_file, versions, summary,
                  discard_header=False):
    """
    Get plain-text result from server and write it to a file.
    """
    try:
        output = open(output_file, 'ab')
    except IOError as (_, message):
        fatal_error(message)

    response = get(result_url)

    in_header = True

    # Buffered read/write
    while True:
        line = response.readline()
        if not line:
            break

        if line.startswith('#'):
            # Comment lines
            if in_header and not discard_header:
                output.write(line)
            if not in_header:
                add_to_summary(versions, summary, line)
        else:
            in_header = False
            output.write(line)

    debug('Result written to: %s' % output_file)

    response.close()
    output.close()


def append_summary(annotation_file, versions, summary):
    """
    Add summary lines to annotation file.
    """
    try:
        annotation = open(annotation_file, 'ab')
    except IOError as (_, message):
        fatal_error(message)

    annotation.write("""
# The following summary is possibly calculated from different result parts
# by the seattle_seq_annotation.py script.
#
""".lstrip())

    for description in versions:
        annotation.write('# %s\n' % description)

    annotation.write('#\n')

    for description, count in summary.items():
        annotation.write('# %s = %d\n' % (description, count))

    annotation.close()


def add_to_summary(versions, summary, line):
    """
    Parse the line and add the results to the summary dictionary or versions
    set.

    Example summary lines:

      # geneDataSource NCBI_hg19 SeattleSeqAnnotation131Version_6.14
      # HapMapFreqType HapMapFreqMinor
      #
      # Count Missense SNPs = 1
      # Count Nonsense SNPs = 0
      # Count SNPs in Splice Sites = 0
      # Count SNPs in Coding Synonymous = 0
      # Count SNPs in Coding (not mod 3) = 0
      # Count SNPs in a UTR = 0
      # Count SNPs near a gene = 8
      # Count SNPs in Introns = 0
      # Count Intergenic SNPs = 19
      #
      # number SNPs in MicroRNAs = 0
      #
      # number accessions coding-synonymous NCBI = 0
      # number accessions missense NCBI = 1
      # number accessions nonsense NCBI = 0
      # number accessions splice-site NCBI = 0
      # number SNPs in dbSNP = 4
      # number SNPs not in dbSNP = 24
      # number SNPs total = 28

    The lines are description=count pairs where description is the key whose
    value in the summary dictionary is incremented by count. Lines without a
    =count suffix are added to the version set.
    """
    parts = line[1:].split('=')
    if len(parts) == 1:
        description = parts[0].strip()
        if description:
            versions.add(description)
    else:
        description, count = '='.join(parts[0:-1]).strip(), int(parts[-1])
        summary[description] += count


def get(url):
    """
    Do a HTTP GET request and return file-like response object.
    """
    request = urllib2.Request(url)
    response = urllib2.urlopen(request)
    return response


def post_multipart(url, parameters={}):
    """
    Do a HTTP POST request encoded as multipart/form-data and return file-like
    response object.
    """
    register_openers()
    datagen, headers = multipart_encode(parameters)
    request = urllib2.Request(POST_URL, datagen, headers)
    response = urllib2.urlopen(request)
    return response


def debug(message):
    if DEBUG: print message


def fatal_error(message):
    print 'Error: %s' % message
    sys.exit(1)


if __name__ == '__main__':
    if len(sys.argv) < 4 or not sys.argv[1] in ['snp', 'indel']:
        print """Annotate variants using SeattleSeq Annotation.

Given a VCF file with variants, submit the file to the SeattleSeq Annotation
web interface [1]. The annotation result is retrieved as tab-separated file.

Usage:
  {command} snp|indel <input.vcf> <mail@domain.com>

The result is written to disk as <input.vcf.annotation>. Temporary data is
written to <input.vcf.part.*> files and removed afterwards.

[1] {url}""".format(command=sys.argv[0], url=BASE_URL)
        sys.exit(1)
    seattle_seq_annotation(sys.argv[1], sys.argv[2], sys.argv[3])
