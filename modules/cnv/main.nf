// CNV main

process compute_cnv {

    cpus 1
	time { 8.minute * task.attempt }
	memory { 5.GB * task.attempt }
	publishDir "${params.output}/cnv/plots", mode: params.publish, pattern: "*_AllChr.png" 
    publishDir "${params.output}/cnv/plots/chr", mode: params.publish, pattern: "*_chr*.png" 
    publishDir "${params.output}/cnv/plots/gene", mode: params.publish, pattern: "*_gene.png" 

    publishDir "${params.output}/cnv/tables", mode: params.publish, pattern: "*.{tsv,igv,seg}" 

    input: 
        path sample_sheet
        tuple val(sample_id), path(rgset_rds)
        path ref_f
        path ref_m
        path ref_mf
        path anno
        val(CNV_focal)

    output:
        tuple path("${sample_id}_*.png"), path("${sample_id}_*_CNVdetail.tsv"), path("${sample_id}_*_metrics.tsv"), path("${sample_id}_*.igv"), path("${sample_id}_*.seg")

    script:
    template "compute_cnv.R"
}

