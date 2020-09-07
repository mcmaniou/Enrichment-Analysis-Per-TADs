# enrichment-analysis-perTADs

This project presents an algorithm for correlation of Chromosomal Locations with functional biological processes (e.g. Gene Ontology Terms, KEGG Pathways, Transcription Factors) taking into account the TADs, that they correspond to. 

## Installing

You can clone the project in your local machine for develpoment and testing purposes using git:

```
git clone https://github.com/mcmaniou/enrichment-analysis-perTADs
```

### Prerequisites

In order to get the project up and running you need to make sure the following packages are installed in your local machine:

```
install.packages(c("tidyverse","ggplot2","data.table","dplyr","tidyr","ggseqlogo","seqinr","httr","jsonlite","xml2","enrichR","stats","purrr","igraph","ggraph","hrbrthemes","extrafont","gridExtra","ggpubr"))

devtools::install_github("nikopech/saveImageHigh")
```

If you install the "extrafont" package make sure to download the fonts with the following command. It only needs to be done once.

```
font_import() 
```

And from Bioconductor:

```
BiocManager::install(c("KEGGREST","pathview","PWMEnrich","PWMEnrich.Hsapiens.background"))
```

## Running the project 

The project consists of one main script:

- enrichmentAnalysis.R

and three additional scripts, that contain the functions called by the main script:

- goPathwayEnrich.R
- motifEnrich.R
- visualization.R

In order to run the project, using the sample input provided in the Datasets folder, run the following command:

```
source("enrichmentAnalysis.R")
```

The output is stored in the Outputs folder.

In the main script, the user can change the following parameters: 
1. The method used for adjustment of the p-values of the enrichment analysis ```p.adjust.method```
(default value is ```"fdr"```, accepted values are ``` "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"```)
2. The threshold of significant p-values ```cut.off```
(default value is 0.05)
3. The minimum number of genes in over-represented terms ```min.genes```
(default value is ```3```)
4. The OS of the local machine ```system```
(default is ```"win"```)

## Data

The raw sample input data used in this project were provided by the Institute of Applied Biosciences (INAB), 
Centre for Research and Technology Hellas (CERTH) and then processed with the open-source tool TnterTADs (https://github.com/nikopech/InterTADs). 

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details


