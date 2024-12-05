process NICHE_COMPASS {
    conda '/home/nima/miniconda3/envs/nichecompass'

    input:
    tuple val(meta), path(validatedSheet)

    // Paths
    val baseDir
    val dbDir
    val tempDir
    val species
     
    val spatial_key 
    val n_neighbors 
    
    // AnnData keys
    val cat_covariates_keys 
    
    // Architecture
    val cat_covariates_embeds_nums 
    val conv_layer_encoder 
    val active_gp_thresh_ratio 
    
    // Trainer
    val n_epochs 
    val n_epochs_all_gps 
    val lr 
    val lambda_edge_recon 
    val lambda_gene_expr_recon 
    val lambda_l1_masked 
    val lambda_l1_addon 
    val edge_batch_size 
    val n_sampled_neighbors 
    val use_cuda_if_available 
    
    // Analysis
    val cell_type_key 
    val latent_leiden_resolution 
    val latent_cluster_key 
    val sample_key 
    val spot_size 

    output:
    path "figure_folder_path.txt", emit: figureFolderPath
    path "model_folder_path.txt", emit: modelFolderPath

    script:
    """
    echo "Base Dir: ${baseDir}"
    echo "DB Dir: ${dbDir}"
    echo "Temp Dir: ${tempDir}"
    echo "Species: ${species}"
    echo "Validated Sheet: ${validatedSheet}"

    NicheCompass.py \
        --samplesheet ${validatedSheet} \
        --base_path "${baseDir}" \
        --data_path "${dbDir}" \
        --tempDir "${tempDir}" \
        --species ${species} \
        --spatial_key 	"${spatial_key}" \
        --n_neighbors 	"${n_neighbors}" \
        --cat_covariates_keys 	"${cat_covariates_keys}" \
        --cat_covariates_embeds_nums 	"${cat_covariates_embeds_nums}" \
        --conv_layer_encoder 	"${conv_layer_encoder}" \
        --active_gp_thresh_ratio 	"${active_gp_thresh_ratio}" \
        --n_epochs 	"${n_epochs}" \
        --n_epochs_all_gps 	"${n_epochs_all_gps}" \
        --lr 	"${lr}" \
        --lambda_edge_recon 	"${lambda_edge_recon}" \
        --lambda_gene_expr_recon 	"${lambda_gene_expr_recon}" \
        --lambda_l1_masked 	"${lambda_l1_masked}" \
        --lambda_l1_addon 	"${lambda_l1_addon}" \
        --edge_batch_size 	"${edge_batch_size}" \
        --n_sampled_neighbors 	"${n_sampled_neighbors}" \
        --use_cuda_if_available 	"${use_cuda_if_available}" \
        --cell_type_key 	"${cell_type_key}" \
        --latent_leiden_resolution 	"${latent_leiden_resolution}" \
        --latent_cluster_key 	"${latent_cluster_key}" \
        --sample_key 	"${sample_key}" \
        --spot_size 	"${spot_size}"


    # Verify the contents of the output files
    echo "Contents of figure_folder_path.txt: "
    cat figure_folder_path.txt
    echo "Contents of model_folder_path.txt:"
    cat model_folder_path.txt
    """
}
