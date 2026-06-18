// merge qc metrics in a single table

process merge_qc_metrics {

    cpus 1
	time { 1.minute * task.attempt }
	memory { 1.GB * task.attempt }
    publishDir "${params.output}/qc", mode: params.publish,  pattern: "*.tsv"

    input:
    path qc_files

    output:
    path "all_qc_metrics.tsv", emit: TSV

    script:
    """
    # take header from first file
    head -n 1 ${qc_files[0]} > all_qc_metrics.tsv

    # append all rows except headers
    for f in ${qc_files}; do
        tail -n +2 \$f >> all_qc_metrics.tsv
    done
    
    # Sort file according to Sample_Name column (first one)
    (head -n 1 all_qc_metrics.tsv && tail -n +2 all_qc_metrics.tsv | sort -k1) > tmp_metrics.tsv

    mv tmp_metrics.tsv all_qc_metrics.tsv
    """
}


process merge_xy_intensities {

    cpus 1
	time { 1.minute * task.attempt }
	memory { 1.GB * task.attempt }
    //publishDir "${params.output}/qc", mode: params.publish,  pattern: "*.tsv"

    input:
    path qc_files

    output:
    path "all_xy_intensities.tsv"

    script:
    """
    # take header from first file
    head -n 1 ${qc_files[0]} > all_xy_intensities.tsv

    # append all rows except headers
    for f in ${qc_files}; do
        tail -n +2 \$f >> all_xy_intensities.tsv
    done

    # Sort file according to Sample_Name column (first one)
    (head -n 1 all_xy_intensities.tsv && tail -n +2 all_xy_intensities.tsv | sort -k1) > tmp_xy.tsv

    mv tmp_xy.tsv all_xy_intensities.tsv
    """
}


process merge_qc_metrics_gs {

    cpus 1
	time { 1.minute * task.attempt }
	memory { 1.GB * task.attempt }
    publishDir "${params.output}/qc", mode: params.publish,  pattern: "*.tsv"

    input:
    path qc_files

    output:
    path "all_qc_metrics_gs.tsv", emit: TSV

    script:
    """
    # take header from first file
    head -n 1 ${qc_files[0]} > all_qc_metrics_gs.tsv

    # append all rows except headers
    for f in ${qc_files}; do
        tail -n +2 \$f >> all_qc_metrics_gs.tsv
    done

    # Sort file according to Sample_Name column 
    (head -n 1 all_qc_metrics_gs.tsv && tail -n +2 all_qc_metrics_gs.tsv | sort -k1) > tmp.tsv

    mv tmp.tsv all_qc_metrics_gs.tsv    
    
    """
}


process merge_cnv_detail {

    cpus 1
	time { 1.minute * task.attempt }
	memory { 1.GB * task.attempt }
    publishDir "${params.output}/cnv", mode: params.publish,  pattern: "*.tsv" 

    input:
    path cnv_detail

    output:
    path "all_CNV_detail.tsv", emit: TSV

    script:
    """
    # take header from first file
    head -n 1 ${cnv_detail[0]} > all_CNV_detail.tsv

    # append all rows except headers
    for f in ${cnv_detail}; do
        tail -n +2 \$f >> all_CNV_detail.tsv
    done

    # Sort file according to Sample_Name column 
    (head -n 1 all_CNV_detail.tsv && tail -n +2 all_CNV_detail.tsv | sort -k8) > tmp.tsv

    mv tmp.tsv all_CNV_detail.tsv
    
    """
}

process merge_cnv_metrics {

    cpus 1
	time { 1.minute * task.attempt }
	memory { 1.GB * task.attempt }
    publishDir "${params.output}/cnv", mode: params.publish,  pattern: "*.tsv" 

    input:
    path cnv_metrics

    output:
    path "all_CNV_metrics.tsv", emit: TSV

    script:
    """
    # take header from first file
    head -n 1 ${cnv_metrics[0]} > all_CNV_metrics.tsv

    # append all rows except headers
    for f in ${cnv_metrics}; do
        tail -n +2 \$f >> all_CNV_metrics.tsv
    done

    # Sort file according to Sample_Name column 
    (head -n 1 all_CNV_metrics.tsv && tail -n +2 all_CNV_metrics.tsv | sort -k1) > tmp.tsv

    mv tmp.tsv all_CNV_metrics.tsv
    
    """
}

process merge_cnv_segment {

    cpus 1
	time { 1.minute * task.attempt }
	memory { 1.GB * task.attempt }
    publishDir "${params.output}/cnv", mode: params.publish,  pattern: "*.seg" 

    input:
    path cnv_segment

    output:
    path "all_CNV_segment.seg", emit: SEG

    script:
    """
    # take header from first file
    head -n 1 ${cnv_segment[0]} > all_CNV_segment.seg

    # append all rows except headers
    for f in ${cnv_segment}; do
        tail -n +2 \$f >> all_CNV_segment.seg
    done

    # Sort file according to Sample_Name column 
    (head -n 1 all_CNV_segment.seg && tail -n +2 all_CNV_segment.seg | sort -k10) > tmp.tsv

    mv tmp.tsv all_CNV_segment.seg
    
    """
}


process merge_mgmt {

    cpus 1
	time { 1.minute * task.attempt }
	memory { 1.GB * task.attempt }
    publishDir "${params.output}/mgmt", mode: params.publish,  pattern: "*.tsv" 

    input:
    path mgmt

    output:
    path "all_mgmt.tsv", emit: TSV

    script:
    """
    # take header from first file
    head -n 1 ${mgmt[0]} > all_mgmt.tsv

    # append all rows except headers
    for f in ${mgmt}; do
        tail -n +2 \$f >> all_mgmt.tsv
    done

    # Sort file according to Sample_Name column 
    (head -n 1 all_mgmt.tsv && tail -n +2 all_mgmt.tsv | sort -k9) > tmp.tsv

    mv tmp.tsv all_mgmt.tsv
    
    """
}


process merge_tumor_purity {

    cpus 1
	time { 1.minute * task.attempt }
	memory { 1.GB * task.attempt }
    publishDir "${params.output}/tumor_purity", mode: params.publish,  pattern: "*.tsv" 

    input:
    path tumor_purity

    output:
    path "all_tumor_purity.tsv", emit: TSV

    script:
    """
    # take header from first file
    head -n 1 ${tumor_purity[0]} > all_tumor_purity.tsv

    # append all rows except headers
    for f in ${tumor_purity}; do
        tail -n +2 \$f >> all_tumor_purity.tsv
    done

    # Sort file according to Sample_Name column 
    (head -n 1 all_tumor_purity.tsv && tail -n +2 all_tumor_purity.tsv | sort -k1) > tmp.tsv

    mv tmp.tsv all_tumor_purity.tsv
    
    """
}
