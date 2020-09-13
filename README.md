# enrichment-analysis-perTADs

This project presents an algorithm for correlation of Chromosomal Locations with functional biological processes (e.g. Gene Ontology Terms, KEGG Pathways, Transcription Factors), taking into account the TADs that they correspond to. 

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

In order to run the project, using the sample input provided in the [Datasets](Datasets) folder, run the following command:

```
source("enrichmentAnalysis.R")
```

If you want to change the input data, the path of the file must be provided in the main script.

The output is stored in the [Outputs](Outputs) folder.

The main script has the following inputs: 
1. ```dbs```: vector with a list of the Enrichr libraries used for the enrichment analysis (default is the vector ```c("GO_Molecular_Function_2018", "GO_Biological_Process_2018", "KEGG_2019_Human")```, for the complete list of the Enrichr libraries type the command ```listEnrichrDbs()```) 
2. ```genes.cover```: vector with the number of genes covered by the Enrichr libraries, which are in the ```dbs``` (taking into account the chosen libraries, default is ```c(11459,14433,7802)```) 
3. ```dir_name``` : the name of the input data folder (default value is ```Datasets```)
4. ```output_folder``` : the name of the outputs folder (default value is ```Outputs```)
5. ```filepath```: the filepath of the input dataset
6. ```p.adjust.method```: the method used for adjustment of the p-values of the enrichment analysis (default value is ```"fdr"```, accepted values are ``` "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"```)
7. ```cut.off```: the threshold of significant p-values (default value is ```0.05```)
8. ```min.genes```: the minimum number of genes in over-represented terms (default value is ```3```)
9. ```system```: the OS of the local machine (default value is ```"win"```)

## Data

The raw sample input data used in this project were provided by the Institute of Applied Biosciences (INAB), 
Centre for Research and Technology Hellas (CERTH) and then processed with the open-source tool [InterTADs](https://github.com/nikopech/InterTADs). 

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.


