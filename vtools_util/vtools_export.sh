


CALLER=$1
VARTABLE=${2-variant}
SAMPLE=$3

less ~/mutation-analysis/vtools_util/formats/mutect_report.header | vtools export variant --format ~/mutation-analysis/vtools_util/formats/mutect2_report.fmt --header - --output mutect-polyp.report --samples 'sample_name like "$SAMPLE"'

python ~/mutation-analysis/vtools_util/mutect2ReportBySampleMut.py -i mutect-polyp.report > mutect-polyp-sample.report
