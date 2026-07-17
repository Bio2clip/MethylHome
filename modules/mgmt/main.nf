
process mgmt {

    cpus 1
	time { 4.minute * task.attempt }
	memory { 3.GB * task.attempt }
	publishDir "${params.output}/mgmt", mode: params.publish, pattern: "*.pdf"
	publishDir "${params.output}/mgmt", mode: params.publish, pattern: "*.tsv"

    input:
    tuple val(sample_id), val(sex), path(rgset_rds)

    output:
    path "MGMT_plot_minfi_${sample_id}.pdf", emit: PDF
    path "${sample_id}_mgmt.tsv", emit: TSV

    script:
    template "mgmt.R" 
}
