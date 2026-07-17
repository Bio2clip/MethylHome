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
        tuple val(sample_id), val(sex), path(rgset_rds)
        path ref_f
        path ref_m
        path ref_mf
        path anno
        val(CNV_focal)

    output:
        path("${sample_id}_*.png"), emit: PNG 
        path("${sample_id}_*_CNVdetail.tsv"), emit: DETAIL
        path("${sample_id}_*_metrics.tsv"), emit: METRICS
        path("${sample_id}_*.igv"), emit: IGV
        path("${sample_id}_*.seg"), emit : SEGMENTS

    script:
    template "compute_cnv.R"
}

