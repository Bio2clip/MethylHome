// merge qc reports in a single file

process merge_qc_reports {

    cpus 1
	time { 5.minute * task.attempt }
	memory { 1.GB * task.attempt }
    publishDir "${params.output}/qc", mode: params.publish,  pattern: "*.pdf"

    input:
    path qc_reports

    output:
    path "all_qc_reports.pdf", emit: PDF

    script:
    """
    pdftk ${qc_reports} cat output all_qc_reports.pdf
    """
}