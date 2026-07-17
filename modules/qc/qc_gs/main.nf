//Compute qc GS 

process compute_qc_gs {

    cpus 1
	time { 5.minute * task.attempt }
	memory { 1.GB * task.attempt }
	//publishDir "${params.output}/qc/sample_gs", mode: params.publish, pattern: "*.tsv"

    input: 
        tuple val(sample_id), val(sex), path(meth_rds)

    output:
        path ("${sample_id}_qc_metrics_gs.tsv")

    script:
    template "compute_qc_gs.R"
}


//Plot qc GS

process plot_qc_gs {

    cpus 1
	time { 5.minute * task.attempt }
	memory { 1.GB * task.attempt }
	publishDir "${params.output}/qc", mode: params.publish, pattern: "*.pdf"

    input:
    path(control_metrics_gs) 

    output:
    path "genome_studio_like_plot.pdf", emit: PDF

    script:
    template "plot_qc_gs.R" 
}