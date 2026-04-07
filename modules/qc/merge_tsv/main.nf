// merge qc metrics in a single table

process merge_qc_metrics {

    cpus 1
	time { 1.minute * task.attempt }
	memory { 1.GB * task.attempt }
    publishDir "${params.output}/qc", mode: params.publish,  pattern: "*.tsv"

    input:
    path qc_files

    output:
    path "all_qc_metrics.tsv", emit: TSV

    script:
    """
    # take header from first file
    head -n 1 ${qc_files[0]} > all_qc_metrics.tsv

    # append all rows except headers
    for f in ${qc_files}; do
        tail -n +2 \$f >> all_qc_metrics.tsv
    done
    """
}


process merge_xy_intensities {

    cpus 1
	time { 1.minute * task.attempt }
	memory { 1.GB * task.attempt }
    //publishDir "${params.output}/qc", mode: params.publish,  pattern: "*.tsv"

    input:
    path qc_files

    output:
    path "all_xy_intensities.tsv"

    script:
    """
    # take header from first file
    head -n 1 ${qc_files[0]} > all_xy_intensities.tsv

    # append all rows except headers
    for f in ${qc_files}; do
        tail -n +2 \$f >> all_xy_intensities.tsv
    done
    """
}