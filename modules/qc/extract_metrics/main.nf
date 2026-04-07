// Extract QC metrics main

process extract_qc_metrics {

    cpus 1
	time { 5.minute * task.attempt }
	memory { 1.GB * task.attempt }
	publishDir "${params.output}/qc/sample_metrics", mode: params.publish, pattern: "*.tsv"

    input:
    tuple val(sample_id), path(meth_rds)

    output:
    tuple val(sample_id), path("${sample_id}_qc_metrics_output.tsv") //Store index to link files

    script:
    template "extract_qc_metrics.R"
}

