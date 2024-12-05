## NicheCompass Nextflow Pipeline

### Overview
This workflow runs NicheCompass which is a graph deep learning method designed to analyze spatial omics data by characterizing cell niches through cellular communication principles. It integrates knowledge of inter- and intracellular interaction pathways to learn an interpretable latent space of cells across multiple tissue samples, facilitating the construction and querying of spatial reference atlases.
For more details check [NicheCompass page]  (https://nichecompass.readthedocs.io/en/latest/index.html).  

#### Step-by-step description of your NicheCompass workflow:

##### CHECK_SAMPLESHEET  

**Input:** A samplesheet CSV file containing information about the input data.  
**Purpose:** Validates the samplesheet to ensure all specified files exist and are accessible.  
**Output:** A validated samplesheet CSV file.  

##### NICHE_COMPASS

**Input:**
- Validated samplesheet  
- Various parameters including paths, species, spatial key, and analysis settings  
**Purpose:** Runs the main NicheCompass analysis, which includes: a. Loading and preprocessing the input data b. Extracting and combining gene programs c. Initializing and training the NicheCompass model d. Identifying niches and performing clustering e. Generating visualizations and analysis results  
**Output:**  
- Paths to figure and model folders containing analysis results  

**Sample Validation:** The pipeline starts by validating the input samplesheet. This ensures that all specified data files exist and are accessible before proceeding with the analysis.

**Data Loading and Preprocessing:** The validated sample data is loaded into AnnData objects. Spatial neighborhood graphs are computed for each sample and combined into a single graph structure.

**Gene Program Extraction:** Gene programs are extracted from various sources (OmniPath, NicheNet, MEBOCOST) and combined. These gene programs represent potential cellular communication pathways.

**Model Initialization:** The NicheCompass model is initialized with the preprocessed data and extracted gene programs.

**Model Training:** The model is trained using the specified parameters, learning to represent the spatial relationships and gene expression patterns in the data.

**Niche Identification:** After training, the model is used to identify cellular niches. This involves clustering cells based on their learned representations.

**Analysis and Visualization: Various analyses are performed on the identified niches, including:**
Plotting niches in latent and physical space
Analyzing niche composition
Generating UMAP embeddings of the data

**Results Output:** The analysis results, including figures and the trained model, are saved to specified output directories.
**Optional TissUUmaps Preparation:** If enabled, the data is prepared for visualization with TissUUmaps by adding necessary grouping information.
### Quick Start
```
cd nf-core
nextflow run main.nf --input samplesheet.csv --outdir results
```

### Parameters   

#### Input/Output Options  
--samplesheet                  Path to the CSV file containing sample information [default: "${myDir}/data/samplesheet_2.csv"]  
--data_path                    Path to the directory containing input data files [default: "${myDir}/data/"]  
--myDir                        Base directory for the analysis [default: '/home/nima/test_base/']  
--tempDir                      Directory for intermediate files [default: "${myDir}/intermediate/${current_timestamp}/"]  
--dbDir                        Directory for database files [default: "${myDir}/data/${current_timestamp}/"]  
--outputDir                    Directory for final results [default: "${myDir}/results/${current_timestamp}/"]  
--mebocost_dir                 Directory containing metabolite enzyme sensor gene programs [default: "${projectDir}/assets/gene_programs/metabolite_enzyme_sensor_gps/"]  

#### Analysis Options  
--species                      Set the species of the samples [default: 'human']  
--spatial_key                  Key for spatial data in the AnnData object [default: 'spatial']  
--n_neighbors                  Number of neighbors to consider in the analysis [default: 4]  
--cat_covariates_keys          Keys for categorical covariates in the AnnData object [default: 'replicate']  
--cell_type_key                Key for cell type annotations [default: 'manual_annotations']  
--sample_key                   Key for sample identification [default: 'replicate']  
--spot_size                    Size of spots in spatial data [default: 30]  

#### Model Architecture Options  
--cat_covariates_embeds_nums   Number of embeddings for categorical covariates [default: 3]  
--conv_layer_encoder           Type of convolutional layer encoder to use [default: 'gcnconv']  
--active_gp_thresh_ratio       Threshold ratio for active gene programs [default: 0.01]  

#### Training Options  
--n_epochs                     Number of epochs for training [default: 400]  
--n_epochs_all_gps             Number of epochs for all gene programs [default: 25]  
--lr                           Learning rate [default: 0.001]  
--lambda_edge_recon            Lambda value for edge reconstruction [default: 500000]  
--lambda_gene_expr_recon       Lambda value for gene expression reconstruction [default: 300]  
--lambda_l1_masked             Lambda value for L1 masked regularization [default: "0."]  
--lambda_l1_addon              Lambda value for L1 addon regularization (de novo GP regularization) [default: "30."]  
--edge_batch_size              Batch size for edge processing [default: 4096]  
--n_sampled_neighbors          Number of sampled neighbors [default: 4]  
--use_cuda_if_available        Whether to use CUDA if available [default: "True"]  

#### Clustering Options  
--latent_leiden_resolution     Resolution for Leiden clustering on latent space [default: 0.2]  

#### Visualization Options  
--enable_tissuumaps            Enable Tissuumaps visualization [default: 'true']  


### Installation   
To run the workflow you need to install:  
- [Nextflow]  (https://www.nextflow.io/docs/latest/install.html)
- [NicheCompass]  (https://nichecompass.readthedocs.io/en/latest/installation.html)

 
