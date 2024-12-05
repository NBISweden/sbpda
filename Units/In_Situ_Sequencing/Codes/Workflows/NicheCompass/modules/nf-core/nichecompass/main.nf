nextflow.enable.dsl=2

// Import necessary Groovy classes
import java.text.SimpleDateFormat
import java.util.Date

// Define base directory
baseDir = '/mypath'

// Get the current timestamp
current_timestamp = new SimpleDateFormat('ddMMyyyy_HHmmss').format(new Date())



process checkSamplesheet {
    input:
    path samplesheet from file(params.samplesheet)

    output:
    tuple val(meta), path('validated_samplesheet.csv') emit: validatedSheet

    conda '${params.env}'

    publishDir params.tempDir, mode: 'copy'  

    script:
    '''
    source activate ${params.env} && \
    python check_samples.py --samplesheet $samplesheet --output validated_samplesheet.csv
    '''
}

process NicheCompass {
    input:
    tuple val(meta), path(validatedSheet)

    output: 
    path 'figure_folder_path' emit: figureFolderPath
    path 'model_folder_path' emit: modelFolderPath

    publishDir params.outputDir, mode: 'copy'   

    conda '${params.env}'

    script:
    '''
    conda activate ${params.env} && \
    python NicheCompass.py \
                        --samplesheet $validatedSheet \
                        --base_path ${params.baseDir} \
                        --data_path ${params.dbDir} \
                        --tempDir ${params.tempDir} \
                        --species ${params.species}
    '''
}
