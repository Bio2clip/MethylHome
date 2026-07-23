#!/usr/bin/env nextflow

// Sample sheet 
params.sample_sheet = ""
if(params.sample_sheet == "") error "ERROR: --sample_sheet must be provided"

// Samples quality parameters 
params.quality_threshold = 1000
params.purity_quality_threshold = 5

// Reference files parameters
params.qc_ref_set = "$projectDir/data/qc_metrics_output_db.tsv"
params.ref_f = "$projectDir/data/epic_geo_ref_f.Rdata"
params.ref_m = "$projectDir/data/epic_geo_ref_m.Rdata"
params.ref_mf = "$projectDir/data/epic_geo_ref_mf.Rdata"
params.anno = "$projectDir/data/annoXY_epic.Rdata"
params.anno_cosmic = "$projectDir/data/anno_epic_cosmic.Rdata"

// Wether to compute focal CNV and so use internet
params.CNV_focal = false

// Output folders and concatenation files names
params.output = "output"
params.all_qc_metrics_file = "all_qc_metrics"
params.all_qc_metrics_gs_file = "all_qc_metrics_gs"
params.all_CNV_detail_file = "all_CNV_detail"
params.all_CNV_detail_cosmic_file = "all_CNV_detail_cosmic"
params.all_CNV_metrics_file = "all_CNV_metrics"
params.all_CNV_segment_file = "all_CNV_segment"
params.all_tumor_purity_file = "all_tumor_purity"
params.all_mgmt_file = "all_mgmt"

params.publish = 'copy'

include { load_idats }                                      from "./modules/load_idats"
include { load_idats_minfi }                                from "./modules/load_idats"
include { sample_sheet_ch }                                 from "./modules/sample_sheet"
include { merge_tsv_files as merge_qc_metrics }             from "./modules/merge_tsv"
include { merge_xy_intensities }                            from "./modules/merge_tsv"
include { merge_tsv_files as merge_qc_metrics_gs }          from "./modules/merge_tsv"
include { merge_tsv_files as merge_cnv_detail }             from "./modules/merge_tsv"
include { merge_tsv_files as merge_cnv_detail_cosmic}       from "./modules/merge_tsv"
include { merge_tsv_files as merge_cnv_metrics }            from "./modules/merge_tsv"
include { merge_tsv_files as merge_cnv_segment }            from "./modules/merge_tsv"
include { merge_tsv_files as merge_tumor_purity }           from "./modules/merge_tsv"
include { merge_tsv_files as merge_mgmt }                   from "./modules/merge_tsv"
include { compute_qc_gs }                                   from "./modules/qc/qc_gs"
include { plot_qc_gs }                                      from "./modules/qc/qc_gs"
include { plot_qc }                                         from "./modules/qc/plot_qc"
include { extract_xy_intensities }                          from "./modules/qc/control_sex"
include { control_sex }                                     from "./modules/qc/control_sex"
include { merge_qc_reports }                                from "./modules/qc/merge_qc_reports"
include { compute_cnv }                                     from "./modules/cnv"
include { tumor_purity }                                    from "./modules/tumor_purity"
include { mgmt }                                            from "./modules/mgmt"


workflow {

    // Create channel with all samples given by the user
    samples_ch = sample_sheet_ch("${params.sample_sheet}")

    // ***********************************************************************************************************************************
    // Load data

    // Read idats 
    load_idats(samples_ch)

    // Read idats with minfi 
    load_idats_minfi(samples_ch)

    // ***********************************************************************************************************************************

    // QC module 

    // Generate pdf report with figures
    qc_reports = plot_qc(load_idats.out, 
                         file(params.qc_ref_set),
                         params.quality_threshold) 

    merge_qc_metrics(qc_reports.TSV.collect(), 
                     params.all_qc_metrics_file,
                     "qc",
                     1)
    // Control sex 
    sex_info = extract_xy_intensities(load_idats.out,
                                      params.quality_threshold)
    all_sex_info = sex_info.collect()
    all_sex_info_tsv = merge_xy_intensities(all_sex_info)
    control_sex(all_sex_info_tsv)

    // Extract raw QC (Genome Studio)
    qc_info = compute_qc_gs(load_idats.out)

    // Merge tsv files by channel
    all_qc_info = qc_info.collect()
    all_qc_gs_tsv = merge_qc_metrics_gs(all_qc_info, 
                                        params.all_qc_metrics_gs_file,
                                        "qc",
                                        1)

    // Plot raw quality control metrics (Genome Studio)
    plot_qc_gs(all_qc_gs_tsv)


    // *********************************************************************************************************************************

    // CNV computation
    cnv = compute_cnv(file(params.sample_sheet), 
                      load_idats_minfi.out, 
                      file(params.ref_f), 
                      file(params.ref_m), 
                      file(params.ref_mf), 
                      file(params.anno),
                      file(params.anno_cosmic), 
                      params.CNV_focal)

    merge_cnv_detail(cnv.DETAIL.collect(),
                     params.all_CNV_detail_file,
                     "cnv",
                     8)
    
    merge_cnv_detail_cosmic(cnv.DETAIL_COSMIC.collect(),
                            params.all_CNV_detail_cosmic_file,
                            "cnv",
                            8)

    merge_cnv_segment(cnv.SEGMENTS.collect(),
                      params.all_CNV_segment_file,
                      "cnv",
                      10)
    
    merge_cnv_metrics(cnv.METRICS.collect(),
                      params.all_CNV_metrics_file,
                      "cnv",
                      1)

    // *********************************************************************************************************************************

    // tumor purity 
    out_tumor_purity = tumor_purity(load_idats_minfi.out,
                                    params.purity_quality_threshold)
    merge_tumor_purity(out_tumor_purity.TSV.collect(),
                       params.all_tumor_purity_file,
                       "tumor_purity",
                       1)

    // *********************************************************************************************************************************

    // MGMT computation
    out_mgmt = mgmt(load_idats_minfi.out)
    merge_mgmt(out_mgmt.TSV.collect(),
               params.all_mgmt_file,
               "mgmt",
               9)

}

