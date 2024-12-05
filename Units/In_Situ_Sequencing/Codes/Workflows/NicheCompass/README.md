## NicheCompass Nextflow Pipeline

### Overview
This workflow runs NicheCompass which is a graph deep learning method designed to analyze spatial omics data by characterizing cell niches through cellular communication principles. It integrates knowledge of inter- and intracellular interaction pathways to learn an interpretable latent space of cells across multiple tissue samples, facilitating the construction and querying of spatial reference atlases.
For more details check [NicheCompass page](https://nichecompass.readthedocs.io/en/latest/index.html).  

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

### Installation  
To run the workflow you need to install:
- [Nextflow](https://www.nextflow.io/docs/latest/install.html)
- [NicheCompass](https://nichecompass.readthedocs.io/en/latest/installation.html)

 
