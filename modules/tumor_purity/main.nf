
process tumor_purity {

    cpus 1
	time { 5.minute * task.attempt }
	memory { 2.GB * task.attempt }
	publishDir "${params.output}/tumor_purity", mode: params.publish, pattern: "*.pdf"
    publishDir "${params.output}/tumor_purity", mode: params.publish, pattern: "*.tsv"


    input:
    tuple val(sample_id), path(rgset_rds)
    
    output:
    path "${sample_id}_tumor_purity.pdf", emit: PDF
    path "${sample_id}_tumor_purity.tsv", emit: TSV

    script:
    template "compute_tumor_purity.R" 
}
