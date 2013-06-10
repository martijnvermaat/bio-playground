#!/bin/bash

# One line of output for every failed output file, no output if everything is
# fine.

for TEST in $(ls | grep -v run.sh); do
    pushd $TEST > /dev/null
    if [ -e gzipped ]; then
        SYNCED_1=$(mktemp --suffix .fq.gz)
        SYNCED_2=$(mktemp --suffix .fq.gz)
        ../../sync_paired_end_reads.py orig.fq.gz filtered_1.fq.gz filtered_2.fq.gz $SYNCED_1 $SYNCED_2 > /dev/null
        diff -q <(zcat synced_1.fq.gz) <(zcat $SYNCED_1) > /dev/null || echo "Failed: $TEST/synced_1.fq.gz"
        diff -q <(zcat synced_2.fq.gz) <(zcat $SYNCED_2) > /dev/null || echo "Failed: $TEST/synced_2.fq.gz"
    else
        SYNCED_1=$(mktemp --suffix .fq)
        SYNCED_2=$(mktemp --suffix .fq)
        ../../sync_paired_end_reads.py orig.fq filtered_1.fq filtered_2.fq $SYNCED_1 $SYNCED_2 > /dev/null
        diff -q synced_1.fq $SYNCED_1 > /dev/null || echo "Failed: $TEST/synced_1.fq"
        diff -q synced_2.fq $SYNCED_2 > /dev/null || echo "Failed: $TEST/synced_2.fq"
    fi
    rm $SYNCED_1 $SYNCED_2
    popd > /dev/null
done
