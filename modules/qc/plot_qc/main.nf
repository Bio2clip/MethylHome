// Plot QC metrics main

process plot_qc {

    cpus 1
	time { 5.minute * task.attempt }
	memory { 1.GB * task.attempt }
	publishDir "${params.output}/qc/sample_plots", mode: params.publish, pattern: "*.pdf"
	publishDir "${params.output}/qc/sample_metrics", mode: params.publish, pattern: "*.tsv"

    input:
    tuple val(sample_id), val(sex), path(meth_rds)   
    path qc_ref_set
    val quality_threshold

    output:
    path "${sample_id}_qc_metrics_output.tsv", emit: TSV
    path "${sample_id}_qc_plot.pdf", emit: PDF

    script:
    template "plot_qc.R" 
}

