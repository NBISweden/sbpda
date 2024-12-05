process checkSamplesheet {
    conda "${params.env}"

    publishDir params.tempDir, mode: 'copy'

    input:
    path samplesheet

    output:
    tuple val(meta), path("validated_samplesheet.csv"), emit: validatedSheet

    script:
    """
    check_samples.py --samplesheet $samplesheet --output validated_samplesheet.csv
    """
}