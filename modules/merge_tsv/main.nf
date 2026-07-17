// merge qc metrics in a single table

process merge_tsv_files {

    cpus 1
	time { 1.minute * task.attempt }
	memory { 1.GB * task.attempt }
    publishDir "${params.output}/${output_folder}", mode: params.publish,  pattern: "*.tsv"

    input:
    path qc_files
    val outfile
    val output_folder
    val colsort

    output:
    path "${outfile}.tsv", emit: TSV

    script:
    """
    # take header from first file
    head -n 1 ${qc_files[0]} > tmp.tsv

    # append all rows except headers
    for f in ${qc_files}; do
        tail -n +2 \$f >> tmp.tsv
    done
    
    # Write header
    head -n 1 tmp.tsv > "${outfile}.tsv"
    
    # Sort data rows acoording to colsort and write them 
    tail -n +2 tmp.tsv | sort -k"${colsort}" >> "${outfile}.tsv"

    """
}


process merge_xy_intensities {

    cpus 1
	time { 1.minute * task.attempt }
	memory { 1.GB * task.attempt }

    input:
    path qc_files

    output:
    path "all_xy_intensities.tsv"

    script:
    """
    # take header from first file
    head -n 1 ${qc_files[0]} > tmp.tsv

    # append all rows except headers
    for f in ${qc_files}; do
        tail -n +2 \$f >> tmp.tsv
    done

    # Sort file according to Sample_Name column (first one)
    (head -n 1 tmp.tsv && tail -n +2 tmp.tsv | sort -k1) > all_xy_intensities.tsv
    # Write header
    head -n 1 tmp.tsv > all_xy_intensities.tsv

    # Sort data rows acoording to colsort and write them 
    tail -n +2 tmp.tsv | sort -k1 >> all_xy_intensities.tsv

    """
}

