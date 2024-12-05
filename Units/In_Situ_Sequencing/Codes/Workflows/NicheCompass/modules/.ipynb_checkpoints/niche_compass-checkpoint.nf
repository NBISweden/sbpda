process NicheCompass {
    conda "${env}"

    publishDir outputDir, mode: 'copy'

    input:
    tuple val(meta), path(validatedSheet)

    output:
    path "figure_folder_path", emit: figureFolderPath
    path "model_folder_path", emit: modelFolderPath

    script:
    """
    python NicheCompass.py \
        --samplesheet $validatedSheet \
        --base_path ${baseDir} \
        --data_path ${dbDir} \
        --tempDir ${tempDir} \
        --species ${species}
    """
}