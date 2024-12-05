include { CHECK_SAMPLESHEET } from './modules/local/check_samplesheet'
include { NICHE_COMPASS } from './modules/local/niche_compass'
include { CHECK_MEBOCOST_FILES } from './modules/local/check_mebocost_files'

//params.mebocost_dir = params.mebocost_dir ?: "${projectDir}/assets/gene_programs/metabolite_enzyme_sensor_gps"

workflow {
    CHECK_MEBOCOST_FILES(params.mebocost_dir)
    CHECK_SAMPLESHEET(
        channel.fromPath(params.samplesheet).map{samplesheet -> tuple([id: samplesheet.baseName], samplesheet)},
        params.data_path
    )
    
    CHECK_SAMPLESHEET.out.validatedSheet.view { meta, file ->
        "Validated samplesheet: ${meta.id} - ${file}"
    }
    NICHE_COMPASS(
        CHECK_SAMPLESHEET.out.validatedSheet,
        params.myDir,
        params.dbDir,
        params.tempDir,
        params.species,
        params.spatial_key,
        params.n_neighbors,
        params.cat_covariates_keys,
        params.cat_covariates_embeds_nums,
        params.conv_layer_encoder,
        params.active_gp_thresh_ratio,
        params.n_epochs,
        params.n_epochs_all_gps,
        params.lr,
        params.lambda_edge_recon,
        params.lambda_gene_expr_recon,
        params.lambda_l1_masked,
        params.lambda_l1_addon,
        params.edge_batch_size,
        params.n_sampled_neighbors,
        params.use_cuda_if_available,
        params.cell_type_key,
        params.latent_leiden_resolution,
        params.sample_key,
        params.spot_size,
        params.enable_tissuumaps = 'True',
        params.mebocost_dir
        
    )

    // Print the returned paths
    NICHE_COMPASS.out.figureFolderPath.subscribe { analysisOutput ->
        println("Figure folder path: ${analysisOutput}")
    }
    NICHE_COMPASS.out.modelFolderPath.subscribe { analysisOutput ->
        println("Model folder path: ${analysisOutput}")
    }
}
