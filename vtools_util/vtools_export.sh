#!/usr/bin/env bash

SAMPLE=${1-"%"}
CODEDIR=$HOME/mutation-analysis/vtools_util
HEADER=$CODEDIR/formats/mutect_report.header
FORMAT=$CODEDIR/formats//mutect2_report.fmt
OUTPUT="mutect-polyp.report"
SAMOUTPUT="mutect-polyp-sample.report"
PAIRS="../pairs.txt"

if [ ! -f "$OUTPUT" ]; then
	echo "vtools export..."
	less $HEADER | vtools export variant --format $FORMAT --header - --output $OUTPUT --samples 'sample_name like "%"'
fi

if [ "$?" != 0 ]; then
	echo "ERROR: vtools export fail."; exit 1
fi

echo "Transform to sample mutation..."
python2.7 $CODEDIR/mutect2ReportBySampleMut.py -i $OUTPUT > $SAMOUTPUT
head -1 $SAMOUTPUT > tmp
grep -w -f <(cut -f2 $PAIRS) $SAMOUTPUT >> tmp
mv tmp $SAMOUTPUT
