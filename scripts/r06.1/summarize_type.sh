#!/bin/bash

MIN_FREQ=0

WIDTH=$1
BETA=$2
N=$3
MIN_FREQ=$4

OUTDIR=w${WIDTH}/b${BETA}
mkdir -p $OUTDIR

#####
# 1. summary stats by variant type:
#####
OUTFILE=$OUTDIR/summary_stats_type.min_frq_${MIN_FREQ}.txt

echo -e "type\twidth\tbeta\tnum_pairs\tD_mean\tD_stdev" > $OUTFILE
# > $OUTDIR/summary_stats_type.txt
for TYPE in all up_up up_down down_down up_neutral down_neutral neutral_neutral
do
    for x in `ls ${OUTDIR}/*.ld.mt.ac.txt.gz`
    do
        if [[ $TYPE == all ]]
        then
            zcat $x | grep -v "^#" \
		| awk -v MIN_FREQ=${MIN_FREQ} -v N=$N '$NF>=MIN_FREQ * N' \
		| cut -f 10 | zstats > $OUTDIR/stats.txt
        elif [[ $TYPE == up_up ]]
        then
            zcat $x | grep -v "^#" | awk '$12>0 && $13>0' \
		| awk -v MIN_FREQ=${MIN_FREQ} -v N=$N '$NF>=MIN_FREQ * N' \
		| cut -f 10 | zstats > $OUTDIR/stats.txt
        elif [[ $TYPE == up_down ]]
        then
            zcat $x | grep -v "^#" | awk '($12*$13 < 0)' \
		| awk -v MIN_FREQ=${MIN_FREQ} -v N=$N '$NF>=MIN_FREQ * N' \
		| cut -f 10 | zstats > $OUTDIR/stats.txt
        elif [[ $TYPE == down_down ]]
        then
            zcat $x | grep -v "^#" | awk '$12<0 && $13<0' \
		| awk -v MIN_FREQ=${MIN_FREQ} -v N=$N '$NF>=MIN_FREQ * N' \
		| cut -f 10 | zstats > $OUTDIR/stats.txt

        elif [[ $TYPE == up_neutral ]]
        then
            zcat $x | grep -v "^#" | awk '($12>0 && $13==0) || ($12==0 && $13>0)' \
		| awk -v MIN_FREQ=${MIN_FREQ} -v N=$N '$NF>=MIN_FREQ * N' \
		| cut -f 10 | zstats > $OUTDIR/stats.txt
        elif [[ $TYPE == down_neutral ]]
        then
            zcat $x | grep -v "^#" | awk '($12<0 && $13==0) || ($12==0 && $13<0)' \
		| awk -v MIN_FREQ=${MIN_FREQ} -v N=$N '$NF>=MIN_FREQ * N' \
		| cut -f 10 | zstats > $OUTDIR/stats.txt
        elif [[ $TYPE == neutral_neutral ]]
        then
            zcat $x | grep -v "^#" | awk '$12==0 && $13==0' \
		| awk -v MIN_FREQ=${MIN_FREQ} -v N=$N '$NF>=MIN_FREQ * N' \
		| cut -f 10 | zstats > $OUTDIR/stats.txt
        fi

        # check if gene is pseudogenized
        x_base=`basename $x .ld.mt.ac.txt.gz`
        x_vcf=${OUTDIR}/${x_base}.vcf.gz

        NLINES=`cat $OUTDIR/stats.txt | grep "num lines" | awk '{ print $NF } END { if (NR==0) { print 0 } }'`
        MEAN=`cat $OUTDIR/stats.txt | grep arith | awk '{ print $NF } END { if (NR==0) { print "NA" } }'`
        STDEV=`cat $OUTDIR/stats.txt | grep stdev | awk '{ print $NF } END { if (NR==0) { print "NA" } }'`
        echo -e "${TYPE}\t${WIDTH}\t${BETA}\t${NLINES}\t${MEAN}\t${STDEV}"
    done
done >> $OUTFILE

######
# 2. Summarize covariance stats
######
COVFILE=$OUTDIR/summary_cov.min_frq_${MIN_FREQ}.txt

# get header
# w2/b1e-1/3043.cl.bin_covar.csv
head -n 1 `ls ${OUTDIR}/*.cl.bin_covar.csv | head -n 1` > $COVFILE
for x in `ls ${OUTDIR}/*.cl.bin_covar.csv`
do
    cat ${x} | sed 1d
done >> $COVFILE

