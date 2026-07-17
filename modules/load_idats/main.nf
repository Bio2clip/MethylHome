// Load idats

process load_idats {

    cpus 1
	time { 5.minute * task.attempt }
	memory { 1.GB * task.attempt }

    input:
    tuple val(sample_id), val(sex), path(idat_green), path(idat_red)

    output:
    tuple val(sample_id), val(sex), path("${sample_id}.rds") 

    script:
    template "load_idats.R" 
}


// Load idats using minfi

process load_idats_minfi {

    cpus 1
	time { 5.minute * task.attempt }
	memory { 1.5.GB * task.attempt }

    input:
    tuple val(sample_id), val(sex), path(idat_green), path(idat_red)

    output:
    tuple val(sample_id), val(sex), path("${sample_id}_minfi.rds") 

    script:
    template "load_idats_minfi.R" 
}