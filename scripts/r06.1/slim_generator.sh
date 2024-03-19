#!/bin/bash

############################################################
#  Program:
#  Author :
############################################################


## BEGIN SCRIPT
usage()
{
    cat <<EOF

usage: $0 OPTIONS

OPTIONS:
    -h      Show this message
    -N      effective population size
    -g      number of generations to run simulation for
    -R      region size (bp)
    -a      effect size (0 is neutral)
    -w      stabilizing selection dist width
    -d      dominance coefficient [0, 0.5]
    -u      mutation rate
    -r      recombination rate
    -o      output dir
    -i      taskid (SLURM)
    -v      Verbose (boolean)

EOF
}

# Show usage when there are no arguments.
if test -z "$1"
then
    usage
    exit
fi

VERBOSE=
FILENAME=

# Check options passed in.
while getopts "h b:w:u:r:o:N:R:i:d:v" OPTION
do
    case $OPTION in
	h)
            usage
            exit 1
            ;;
	R)
	    REG_SIZE=$OPTARG
	    ;;
	b)
            BETA=$OPTARG
            ;;
	w)
	    WIDTH=$OPTARG
	    ;;
	d)
	    DOM_COEFF=$OPTARG
	    ;;
	u)
            MU=$OPTARG
            ;;
	r)
            RECOMBINATION=$OPTARG
            ;;
	N)
            POP_SIZE=$OPTARG
            ;;
	o)
            OUTDIR=$OPTARG
            ;;
	i)
	    TASK_ID=$OPTARG
	    ;;
        v)
            VERBOSE=1
            ;;
        ?)
            usage
            exit
            ;;
    esac
done

OUTPREFIX=${OUTDIR}/$TASK_ID

# run slim
time slim \
    -s $TASK_ID \
    -d reg_length=${REG_SIZE} \
    -d dom_coeff=${DOM_COEFF} \
    -d recomb_rate=${RECOMBINATION} \
    -d mu=${MU} \
    -d pop_size=${POP_SIZE} \
    -d width=${WIDTH} \
    -d beta=${BETA} \
    -d "out_prefix='${OUTPREFIX}'" \
    slim_cmd.slim

cat ${OUTPREFIX}.vcf | vawk --header '{ $3=I$MID ; print }' | bgzip -c > ${OUTPREFIX}.vcf.gz && rm ${OUTPREFIX}.vcf
tabix -f -p vcf ${OUTPREFIX}.vcf.gz;

emeraLD --in ${OUTPREFIX}.vcf.gz --stdout --phased --extra --dstats | gzip -c > ${OUTPREFIX}.ld.txt.gz

zcat ${OUTPREFIX}.ld.txt.gz \
    | zjoin -r -a stdin -b <(zcat ${OUTPREFIX}.vcf.gz | vawk '{ print $3, I$S }' OFS="\t") -1 3 -2 1 \
    | cut -f -11,13 \
    | zjoin -a stdin -b <(zcat ${OUTPREFIX}.vcf.gz | vawk '{ print $3, I$S }' OFS="\t") -1 6 -2 1 \
    | cut -f -12,14 \
    | zjoin -r -a stdin -b <(zcat ${OUTPREFIX}.vcf.gz | vawk '{ print $3, I$AC }' OFS="\t") -1 3 -2 1 \
    | cut -f -13,15 \
    | zjoin -r -a stdin -b <(zcat ${OUTPREFIX}.vcf.gz | vawk '{ print $3, I$AC }' OFS="\t") -1 6 -2 1 \
    | cut -f -14,16 \
    | gzip -c \
    > ${OUTPREFIX}.ld.mt.ac.txt.gz

# calculate snp-pair effect
zcat ${OUTPREFIX}.vcf.gz \
    | ./parse_vcf_gravel.py \
    /dev/stdin \
    ${OUTPREFIX}
