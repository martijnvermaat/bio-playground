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
# written to <input.vcf.annotation.tmp> and removed afterwards.
#
# Requires the poster Python library [2].
#
# Todo: Split VCF file if it contains more than the maximum number of variants
# accepted by SeattleSeq and submit the parts separately.
#
# [1] http://snp.gs.washington.edu/SeattleSeqAnnotation131/
# [2] http://atlee.ca/software/poster/
#
# 2011-02-10, Martijn Vermaat <m.vermaat.hg@lumc.nl>


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
           'polyPhen',
#           'distanceToSplice',
#           'microRNAs',
           'clinicalAssociation']

# Maximum number of variants SeattleSeq accepts
MAX_VARIANTS = 2000000

# SeattleSeq Annotation location
BASE_URL = 'http://snp.gs.washington.edu/SeattleSeqAnnotation131/'
POST_URL = BASE_URL + 'BatchQueryServlet'


import sys
import os
import time
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
    monitor_url, result_url = submit_vcf_file(mode, vcf_file, address)
    wait_for_result(monitor_url)
    write_result_file(result_url, vcf_file + '.annotation')


def submit_vcf_file(mode, vcf_file, address):
    """
    Submit a VCF file to the SeattleSeq Annotation server and return as a
    tuple the url to monitor the job and the url to get the result.

    The {mode} argument can be 'snp' (default) or 'indel'.
    """
    format = 'VCFIndel' if mode == 'indel' else 'VCF'

    auto_file = vcf_file + '.annotation.tmp'
    create_auto_file(vcf_file, auto_file)

    try:
        auto = open(auto_file, 'r')
    except IOError as (_, message):
        fatal_error(message)

    parameters = [('genotypeSource',   'FileInput'),
                  ('EMail',            address),
                  ('GenotypeFile',     auto),
                  ('fileFormat',       format),
                  ('outputFileFormat', 'original'),
                  ('geneData',         'NCBI'),
                  ('HapMapFreqType',   'HapMapFreqMinor'),
                  ('gFetch',           'Display Genotypes')]
    for column in COLUMNS:
        parameters.append( ('columns', column) )

    debug('Submitting for %s annotation: %s' % (mode, vcf_file))

    # Response contains result url and monitor url separated by a comma
    response = post_multipart(POST_URL, parameters)
    urls = response.readline().strip().split(',')
    response.close()

    auto.close()
    os.unlink(auto_file)

    if not len(urls) == 2:
        fatal_error('Could not read urls from submit response.')

    result_url, monitor_url = urls

    debug('Monitor url: %s' % monitor_url)
    debug('Result url: %s' % result_url)

    return monitor_url, result_url


def create_auto_file(vcf_file, auto_file):
    """
    Copy the VCF file to another file, prepending it with a line needed by
    the SeattleSeq Annotation server.
    """
    try:
        vcf = open(vcf_file, 'r')
    except IOError as (_, message):
        fatal_error(message)

    try:
        auto = open(auto_file, 'w')
    except IOError as (_, message):
        fatal_error(message)

    # Header line for automated processing
    auto.write('# autoFile testAuto.txt\n')

    line_count = 0

    # Buffered read/write
    while True:
        line = vcf.readline()
        if not line:
            break
        auto.write(line)
        line_count += 1
        if line_count > MAX_VARIANTS:
            fatal_error('VCF file contains more variants than accepted by '
                        'SeattleSeq (%d).' % MAX_VARIANTS)

    # Todo: Ignore header lines (but this will be fixed when we implement
    # splitting of large files anyway).
    if line_count == 0:
        fatal_error('VCF file contains no variants.')

    debug('Copied VCF file to temporary file.')

    vcf.close()
    auto.close()


def wait_for_result(monitor_url):
    """
    Monitor the job for progress. Return when the job is completed.
    """
    waiting = 0
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
            debug('Waiting: 0%')
            continue

        debug('Waiting: %s / %s' % (processed, total))

        # See if we are done
        if processed == total:
            break


def write_result_file(result_url, output_file):
    """
    Get plain-text result from server and write it to a file.
    """
    try:
        output = open(output_file, 'wb')
    except IOError as (_, message):
        fatal_error(message)

    response = get(result_url)

    # Buffered read/write
    while True:
        line = response.readline()
        if not line:
            break
        output.write(line)

    debug('Result written to: %s' % output_file)

    response.close()
    output.close()


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
written to <input.vcf.annotation.tmp> and removed afterwards.

[1] {url}""".format(command=sys.argv[0], url=BASE_URL)
        sys.exit(1)
    seattle_seq_annotation(sys.argv[1], sys.argv[2], sys.argv[3])
