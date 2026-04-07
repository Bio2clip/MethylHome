// Control sex main

process extract_xy_intensities {

    cpus 1
	time { 5.minute * task.attempt }
	memory { 1.GB * task.attempt }
	publishDir "${params.output}/qc/sample_sex", mode: params.publish, pattern: "*.tsv"

    input: 
        path sample_sheet
        tuple val(sample_id), path(meth_rds)

    output:
        path "${sample_id}_xy_intensities.tsv"

    script:
    template "compute_xy_intensities.R"
}


process control_sex {
    cpus 1
	time { 2.minute * task.attempt }
	memory { 1.GB * task.attempt }
	publishDir "${params.output}/qc", mode: params.publish

    input: 
        path sex_info

    output:
        path "control_sex_report.pdf", emit: PDF
        path "all_predicted_sex.tsv", emit: TSV

    script:
    template "control_sex.R"

}
