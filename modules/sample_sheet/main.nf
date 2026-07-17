def sample_sheet_ch(sample_sheet_path) {
    // Raise an error if the sample sheet file can't be found
    sample_sheet = file("${sample_sheet_path}")

    // Read all lines
    lines = sample_sheet.readLines()

    // Find the index of the row containing "Sample_Name"
    header_index = lines.findIndexOf { it.split(',').contains('Sample_Name') }
    sep = ','

    // if file is not csv and does not contains Sample_Name 
    if (header_index == -1) {
        //sep: '\t' for tsv file
        header_index = lines.findIndexOf { it.split('\t').contains('Sample_Name') }
        sep = '\t'
        // if file is not csv nor tsv and does not contains Sample_Name 
        if (header_index == -1) {
            error "Column 'Sample_Name' not found in ${params.sample_sheet}"
        }
    }

    // Create channel with all samples using dynamic skip
    samples_ch = Channel.fromPath(sample_sheet)
                        .splitCsv(header: true, skip: header_index, sep: sep)
                        .map { row -> tuple(row.Sample_Name,
                                            row.Gender ?: row.Sex, //Elvis operator, if Gender does not exist take Sex
                                            file("${row.file_path}_Grn.idat"),
                                            file("${row.file_path}_Red.idat")) }

    // Collect sample names
    names_ch = samples_ch.map { it[0] }.collect()

    // Verify uniqueness in sample names
    names_ch.subscribe { names ->
        if (names.size() != names.unique().size()) {
            error "Please enter unique names in ${params.sample_sheet}"
        }
    }

    return samples_ch
}