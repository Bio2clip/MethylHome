#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Container used to run the script
//params.container = "MethylHome_QC-v0.1_20260322-144237.sif"

// Default values for all parameters
params.output = ""
params.sample_sheet = ""
params.qc_ref_set = 'data/qc_metrics_output_db.tsv'
params.ref_f = 'data/epic_geo_ref_f.Rdata'
params.ref_m = 'data/epic_geo_ref_m.Rdata'
params.ref_mf = 'data/epic_geo_ref_mf.Rdata'
params.anno = 'data/annoXY_epic.Rdata'
params.CNV_focal = 'FALSE'
params.publish = 'copy'

include { load_idats }                    from "./modules/load_idats"
include { load_idats_minfi }              from "./modules/load_idats"
include { extract_qc_metrics }            from "./modules/qc/extract_metrics"
include { merge_qc_metrics }              from "./modules/qc/merge_tsv"
include { merge_xy_intensities }          from "./modules/qc/merge_tsv"
include { merge_qc_metrics_gs }           from "./modules/qc/merge_tsv"
include { merge_cnv_detail }              from "./modules/qc/merge_tsv"
include { merge_cnv_metrics }             from "./modules/qc/merge_tsv"
include { merge_cnv_segment }             from "./modules/qc/merge_tsv"
include { merge_tumor_purity }            from "./modules/qc/merge_tsv"
include { merge_mgmt }                    from "./modules/qc/merge_tsv"
include { compute_qc_gs }                 from "./modules/qc/qc_gs"
include { plot_qc_gs }                    from "./modules/qc/qc_gs"
include { plot_qc }                       from "./modules/qc/plot_qc"
include { extract_xy_intensities }        from "./modules/qc/control_sex"
include { control_sex }                   from "./modules/qc/control_sex"
include { merge_qc_reports }              from "./modules/qc/merge_qc_reports"
include { compute_cnv }                   from "./modules/cnv"
include { tumor_purity }                  from "./modules/tumor_purity"
include { mgmt }                          from "./modules/mgmt"


workflow {

    //container = params.container

    // Raise an error if the expected parameters were not found
    if (!params.output){error "Must provide parameter 'output'"}

    // Raise an error if the sample sheet file can't be found
    sample_sheet = file("${params.sample_sheet}", checkIfExists: true)

    // Read all lines
    lines = sample_sheet.readLines()

    // Find the index of the row containing "Sample_Name"
    //sep: '\t' for tsv file
    header_index = lines.findIndexOf { it.split(',').contains('Sample_Name') }
    sep = ','

    // Check if file contains Sample_Name
    if (header_index == -1) {
        
        header_index = lines.findIndexOf { it.split('\t').contains('Sample_Name') }
        sep = '\t'
        if (header_index == -1) {
            error "Column 'Sample_Name' not found in ${params.sample_sheet}"
        }
    }

    // Create channel with all samples using dynamic skip
    samples_ch = Channel.fromPath(sample_sheet)
        .splitCsv(header: true, skip: header_index, sep: sep)
        .map { row -> tuple(row.Sample_Name,file("${row.file_path}_Grn.idat"),file("${row.file_path}_Red.idat")) }

    // Collect sample names
    names_ch = samples_ch.map { it[0] }.collect()

    // Verify uniqueness in sample names
    names_ch.subscribe { names ->
        if (names.size() != names.unique().size()) {
            error "Please enter unique names in ${params.sample_sheet}"
        }
    }

    // ***********************************************************************************************************************************
    // Load data

    // Read idats 
    load_idats(samples_ch)

    // Read idats with minfi 
    load_idats_minfi(samples_ch)

    // ***********************************************************************************************************************************

    // QC module 

    // Extract QC metrics (BeadArray computation)
    qc_results = extract_qc_metrics(load_idats.out)

    // Merge BeadArry-like QC metrics
    all_qc_tsvs = qc_results.map { sample_id, tsv_file -> tsv_file }  // extract only the TSV file
                            .collect()

    merge_qc_metrics(all_qc_tsvs)

    // Join outputs to ensure they are passed on together to the plot_QC process
    plot_input = qc_results.join(load_idats.out)

    // Generate pdf report with figures
    qc_reports = plot_qc(plot_input, file(params.qc_ref_set)) 

    // Control sex
    sex_info = extract_xy_intensities(file(params.sample_sheet), load_idats.out)
    all_sex_info = sex_info.collect()
    all_sex_info_tsv = merge_xy_intensities(all_sex_info)
    control_sex(all_sex_info_tsv)

    // Extract raw QC (Genome Studio)
    qc_info = compute_qc_gs(load_idats.out)

    // Merge tsv files by channel
    all_qc_info = qc_info.collect()
    all_qc_gs_tsv = merge_qc_metrics_gs(all_qc_info)

    // Plot raw quality control metrics (Genome Studio)
    plot_qc_gs(all_qc_gs_tsv)


    // *********************************************************************************************************************************

    // CNV computation
    cnv = compute_cnv(file(params.sample_sheet), load_idats_minfi.out, file(params.ref_f), file(params.ref_m), file(params.ref_mf), file(params.anno), params.CNV_focal)
    
    detail = cnv.map { png, detail, metrics, igv, segments -> detail } 
                .collect()

    merge_cnv_detail(detail)

    segments = cnv.map { png, detail, metrics, igv, segments -> segments } 
                  .collect()

    merge_cnv_segment(segments)

    metrics = cnv.map { png, detail, metrics, igv, segments -> metrics } 
                 .collect()
    
    merge_cnv_metrics(metrics)

    // *********************************************************************************************************************************

    // tumor purity 
    out_tumor_purity = tumor_purity(load_idats_minfi.out)
    tumor_purity_tsv = out_tumor_purity.TSV 
                                       .collect()
    merge_tumor_purity(tumor_purity_tsv)

    // *********************************************************************************************************************************

    // MGMT computation
    out_mgmt = mgmt(load_idats_minfi.out)
    mgmt_tsv = out_mgmt.TSV 
                       .collect()
    merge_mgmt(mgmt_tsv)

}

