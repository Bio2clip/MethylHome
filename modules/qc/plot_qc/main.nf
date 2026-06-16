// Plot QC metrics main

process plot_qc {

    cpus 1
	time { 5.minute * task.attempt }
	memory { 1.GB * task.attempt }
	publishDir "${params.output}/qc/sample_plots", mode: params.publish, pattern: "*.pdf"

    input:
    tuple val(sample_id), path(qc_tsv), path(meth_rds)   
    path qc_ref_set

    output:
    path "${sample_id}_qc_plot.pdf", emit: PDF

    script:
    template "plot_qc.R" 
}

