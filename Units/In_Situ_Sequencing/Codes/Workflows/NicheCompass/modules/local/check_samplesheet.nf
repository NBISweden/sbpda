process CHECK_SAMPLESHEET {
    conda '/home/nima/miniconda3/envs/nichecompass'

    input:
    tuple val(meta), path(samplesheet)
    path(data_path)

    output:
    tuple val(meta), path("validated_samplesheet.csv"), emit: validatedSheet

    script:
    """
    check_samples.py --samplesheet $samplesheet --data_path $data_path --output validated_samplesheet.csv
    """
}