#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// The container used to run the script
params.container = "MethylHome_QC-v0.1_20260322-144237.sif"

// The default values for all parameters
params.output = ""
params.sample_sheet = ""
params.database = 'data/qc_metrics_output_db.tsv'
params.publish = 'copy'

include { load_idats }                    from "./modules/load_idats"
include { extract_qc_metrics }            from "./modules/qc/extract_metrics"
include { merge_qc_metrics }              from "./modules/qc/merge_tsv"
include { merge_xy_intensities }          from "./modules/qc/merge_tsv"
include { plot_qc }                       from "./modules/qc/plot_qc"
include { extract_xy_intensities }        from "./modules/qc/control_sex"
include { control_sex }                   from "./modules/qc/control_sex"
include { merge_qc_reports }              from "./modules/qc/merge_qc_reports"

workflow {

    container = params.container

    // Raise an error if the expected parameters were not found
    if (!params.output){error "Must provide parameter 'output'"}

    // Raise an error if the sample sheet file can't be found
    sample_sheet = file("${params.sample_sheet}", checkIfExists: true)

    // Read first line to check if there is a header
    first_line = Channel.fromPath(sample_sheet)
                    .splitCsv()
                    .map { row -> row[0] }
                    .first()
                    .view()

    //lines = file(sample_sheet).text.readLines()

    //data_index = lines.findIndexOf { it.contains('Sample_Name') }

    if (first_line.contains('[Header]')){
        // Create channel with all samples
        samples_ch = Channel.fromPath(sample_sheet)
                            .splitCsv(header : true, skip: 7)
                            .map { row -> tuple(row.Sample_Name, file("${row.file_path}_Grn.idat"),file("${row.file_path}_Red.idat"))}
                            .view()
    } else {
        // Create channel with all samples
        samples_ch = Channel.fromPath(sample_sheet)
                            .splitCsv(header : true)
                            .map { row -> tuple(row.Sample_Name, file("${row.file_path}_Grn.idat"),file("${row.file_path}_Red.idat"))}
    }

    // Read idats 
    load_idats(samples_ch)

    // Extract QC metrics
    qc_results = extract_qc_metrics(load_idats.out)

    // Join outputs to ensure they are passed on together to the plot_QC process
    plot_input = qc_results.join(load_idats.out)

    // Generate pdf report with figures
    qc_reports = plot_qc(plot_input, file(params.database)) 

    // Control sex
    sex_info = extract_xy_intensities(file(params.sample_sheet), load_idats.out)
    all_sex_info = sex_info.collect()
    all_sex_info_tsv = merge_xy_intensities(all_sex_info)
    control_sex(all_sex_info_tsv)

    all_qc_tsvs = qc_results.map { sample_id, tsv_file -> tsv_file }  // extract only the CSV file
                            .collect()

    merge_qc_metrics(all_qc_tsvs)
}