# iSensors <img src="https://img.shields.io/badge/R-package-blue" alt="R package" height="24">

**iSensors** is a package for analyzing signaling activity scores of gene panels based on single-cell RNA-seq data stored in Seurat objects. The package allows you to compute signaling activity scores for specified gene panels and store the results as new assays.

## Installation

You can install the package from GitHub using `devtools`:

```R
# Make sure you have devtools installed
install.packages("devtools")

# Install the iSensor package from main
devtools::install_github("MironovaLab/iSensors")
# or from other branch
devtools::install_github("MironovaLab/iSensors", ref = "devel")
```
In case dealing with a mistake like:
```
Failed to install 'unknown package' from GitHub
```
You need to go to GitHub Settings -> Developer Settings and generate your token.
Run
```
usethis::edit_r_environ()
```
and copy your token to the GITHUB_PAT= your token/

## Quick Start

Here's an example of how to load the package, check its version, explore documentation, and run a basic analysis:

```R
# Load the package
library(iSensors)

# Check installed version
packageVersion("iSensors")

# Access general help page
help("iSensors")

# Load test data (Seurat object)
testData <- readRDS("testData/testSeurData.rds")

# Load test panel with meta-panels
testPanel <- LoadSensors(setName = 'testPanelSet', species = 'AT', hormone = 'aux', customPanels = TRUE,
                          randomInfo = list('n' = 3, 'sizes' = c(100, 200, 300), majortrend = TRUE),
                          metaPanels = list(
                            'meta1' = list('srcPanels' = c("AT_aux_cis_DR5_ARF1", "AT_aux_cistrans_DR5_ARF5_2_up"), rule = mean),
                            'meta2' = list('srcPanels' = c("AT_aux_cis_DR5_ARF1", "AT_aux_cistrans_DR5_ARF5_2_up"), rule = prod))
                          )

# Calculate signaling scores for selected panels
result <- CalcSensors(testData,
                      seurLayer = "data",
                      panelSet = testPanel,
                      signals = c("mean_normed", "median"))
```

iSensorsTransPanelCreate
=========================

Introduction
--------------------------------------------------

The **iSensorsTransPanelCreate.R** function generates an object for a trans-type gene panel compatible with iSensors R package.

The function receives a list of gene IDs as input (either as a vector or as a txt file), and the the corresponding trivial gene names in txt format (optional).

The function generates a list of three items. The first is a list **genes** containing gene IDs. The second is a data frame **gene_metadata** containing gene IDs and the corresponding short and full gene names. The third is a data frame **panel_metadata** containing the information about the species, gene panel type, gene panel description, creation date. The function writes the panel to an object in the current environment and saves it as an rda file in the current working directory.

Requirements
------------

### R packages

- stringr
- magrittr
- dplyr
- Biostrings

Installation of this packages in **R** is carried out using the commands

```
install.packages('stringr')
install.packages('magrittr')
install.packages('dplyr')

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Biostrings")
```

Usage
------

#### Parameters of function

| Parameter | Description |
| ----------- | ----------- |
| panel_name | The name of the panel to create. The function will create an object with this name in the current environment and save the rda file with this name in the current working directory. |
| gene_list | List of genes. Can be supplied as a vector `c('AT1G01010', 'AT1G01030', 'AT1G01040')`, or as a .txt file.   |
| species | Name of species for which the trans panel is created.|
| trivial_names_file | List of trivial gene names as a .txt file. (optional) |
| panel_description    | Description of the panel in natural language.|

#### Gene list TXT-file format

 ```
AT1G01010
AT1G01030
AT1G01040
```
#### Trivial names list TXT-file format

Tab should be used as a delimiter

| gene | gene_name | gene_trivial_name |
| ----------- | ----------- | ----------- |
|AT1G01010 | NAC001 | NAC domain containing protein 1 |
|AT1G01020 | ARV1 | none |
|AT1G01030 | NGA3 | NGATHA3 |

#### iSensor panel object

This is an example of panel with two genes 

- *genes*

`  "AT1G01020" "AT1G01060" `

- *gene_metadata*

| Gene | Gene name | Gene full name |
| ---- | ---------- | --------------|
| AT1G01020 |ARV1|  none| 
AT1G01060| LHY | LATE ELONGATED HYPOCOTYL |

- *panel_metadata*

| Species | Panel type | Panel description | Date Created |
| ---- | ---------- | --------------|----|
| Arabidopsis thaliana | trans | This is an example of panel with two genes | 2025-06-20 |

  


#### Usage examples

- *Arabidopsis thaliana*, input as vector.

` iSensorsTransPanelCreate(panel_name = 'trans_panel', gene_list = c('AT1G01010', 'AT1G01030', 'AT1G01040'), species = 'Arabidopsis thaliana', panel_description = 'This is an example of trans panel') `

- *Arabidopsis thaliana*, input as txt-file.

` iSensorsTransPanelCreate(panel_name = 'trans_panel' ,gene_list = 'gene_list.txt', species = 'Arabidopsis thaliana', panel_description = 'This is an example of trans panel') `

- *Arabidopsis thaliana*, input as txt-file, with trivial names list.

` iSensorsTransPanelCreate(panel_name = 'trans_panel', gene_list = 'gene_list.txt', species = 'Arabidopsis thaliana', trivial_names_file = 'Arabidopsis_trivial_names_example.txt', panel_description = 'This is an example of trans panel') `

iSensorsCisTransPanelCreate
=========================

Introduction
--------------------------------------------------

The **iSensorsCisTransPanelCreate.R** function performs the recognition of binding sites in promoters using positional weight matrices and generates an object for a cis-trans type panel, if the recognition is limited to differentially expressed genes, and cis-trans type gene panel if not. Format for cis/cis-trans panels is compatible with iSensors R package.

The function receives (1) a set of promoters in FASTA format, (2) a positional probability matrix in Homer format, (3) a list of RNA-seq experiments with logFC and adjusted p-values for each gene (optional) and (4) the the corresponding trivial gene names in txt format (optional).

The function generates a list of three items. The first is a list **genes** containing gene IDs. The second is a data frame **gene_metadata** containing gene IDs, recognized sites, coordinates of each site in the genome, location (forward or reverse strand) and distance of each site relative to the TSS, and the corresponding short and full gene names. The third is a data frame **panel_metadata** containing the information about the species, gene panel type, gene panel description, creation date. 
The function writes the panel to an object in the current environment and saves it as an rda file in the current working directory.

Requirements
------------

### R packages

- universalmotif
- Biostrings
- stringr
- magrittr
- dplyr
- purrr

Installation of this packages in **R** is carried out using the commands

```
install.packages('stringr')
install.packages('magrittr')
install.packages('dplyr')
install.packages('purrr')

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("Biostrings")
BiocManager::install("universalmotif")
```

Usage
------

#### Parameters of function

| Parameter | Description |
| ----------- | ----------- |
| panel_name | The name of the panel to create. The function will create an object with this name in the current environment and save the rda file with this name in the current working directory. |
| species | Name of species for which the trans panel is created.|
| promoters_set | Path to a FASTA file containing promoter sequences. For the header structure, see *Promoters set format*|
| ppm | Path to a positional probability matrix file in Homer format. For the description of format, see *PPM format* |
| deg_list | Path to a file containing logFC and adjusted p-value information obtained from a set of transcriptomic experiments that determine the differential expression status of a gene. For the description of format, see *DEG list format* |
| min_dataset_number | The minimum number of transcriptomics experiments in which a gene is differentially expressed to be considered differentially expressed in the cis-trans panel. For example, if `min_dataset_number = 1`, then a gene must be differentially expressed in at least one experiment to be considered differentially expressed in the cis-trans panel. |
| trivial_names_file | List of trivial gene names as a .txt file. (optional) |
| panel_type | Type of panel. Could be `cis`, `UP` and `DOWN`. |
|  transcriptomes_info | Information about the used transcriptomes, written by the user in free form in natural language.|

#### Promoters set format

The function should receive a set of promoters as a **FASTA** file. Each individual promoter must be a separate sequence in the file. The header *must* contain the following information, separated by the `-` sign:

- Gene ID;
- Gene orientation: 1 if gene is located on forward strand, 0 if gene is located on reverse strand;
- Chromosome;
- Promoter absolute start coordinate;
- Promoter absolute end coordinate.

##### Example of promoter in promoters set

 ```
>AT1G01010-1-1-2130-3630
AACATTTCAAACCACTTGTTCTCTTTTATGTTTTGGTAAGAGCTATCTTCTAAATTTATAATACGCATAAATTCAAAAGTAAAAGAAAATTTTGGTCATGAATGTTGTTTAAGTCATTTGGAGATACGAAATCAAATCTCCTTGTAGATTTTGTTTTTAGAATGTCGTTCCTTTTTCATCATCTTAGCTATATCTACAGCTATATATCCTATCTTTAAACCTATATTATTTTTTCCTCTCTTCACCAAAGCCATGTTTTTTAGTTGTGGCGAAAAATAAGAAATCCATACATCAACATATCGCTTTCGTTACCTTAAATTTTGGCTTGTTATGAAGGCATGTCATAACGTTTCTAGTCACAACTCACAAGCATACCAACGACCATGATAAATCCAAAAAGTAGAAACAATCTATTATCTAAACCCCCAAAAGACAAAAGAAAAAAGTAGAAAGAAAAGGTAGGCAGAGATATAATGCTGGTTTTATTTGTTTGTTAAAAGATATTGCTATTTCTGCCAATATTAAAACTTCACTTAGGAAGACTTGAACCTACCACACGTTAGTGACTAATGAGAGCCACTAGATAATTGCATGCATCCCACACTAGTACTAATTTTCTAGGGATATTAGAGTTTTCTAATCACCTACTTCCTACTATGTGTATGTTATCTACTGGCGTGGATGCTTTTAAAGATGTTACGTTATTATTTTGTTCGGTTTGGAAAACGGCTCAATCGTTATGAGTTCGTAAGACACATACATTGTTCCATGATAAAATGCAACCCCACGAACCATTTGCGACAAGCAAAACAACATGGTCAAAATTAAAAGCTAACAATTAGCCAGCGATTCAAAAAGTCAACCTTCTAGATGGATTTAACAACATATCGATAGGATTCAAGATTAAAAATAAGCACACTCTTATTAATGTTAAAAAACGAATGAGATGAAAATATTTGGCGTGTTCACACACATAATCTAGAAGACAGATTCGAGTTGCTCTCCTTTGTTTTGCTTTGGGAGGGACCCATTATTACCGCCCAGCAGCTTCCCAGCCTTCCTTTATAAGGCTTAATTTATATTTATTTAAATTTTATATGTTCTTCTATTATAATACTAAAAGGGGAATACAAATTTCTACAGAGGATGATATTCAATCCACGGTTCACCCAAACCGATTTTATAAAATTTATTATTAAATCTTTTTTAATTGTTAAATTGGTTTAAATCTGAACTCTGTTTACTTACATTGATTAAAATTCTAAACCATCATAAGTAAAAAATAATATGATTAAGACTAATAAATCTTAATAGTTAATACTACTCGGTTTACTACATGAAATTTCATACCATCAATTGTTTTAATAATCTTTAAAATTGTTAGGACCGGTAAAACCATACCAATTAAACCGGAGATCCATATTAATTTAATTAAGAAAATAAAAATAAAAGGAATAAATTGTCTTATTTAAACGCTGACTTCACTGTCTTCCTCCCTCC
```
#### PPM format

The function receives a positional probability matrix as input and transforms it into a positional weight matrix. The matrix format is Homer. Its structure is:
- a header starting with a greater than sign
- four columns separated by tabs, containing the probabilities of finding four letters of the genetic alphabet in each position. The number of rows is equal to the length of the matrix (excluding the header).

##### Example of PPM

 ```
>ARF1 - MA0942.1      
0.312 0.207 0.22  0.262
0.008 0.088 0.004 0.901
0.002 0.002 0.993 0.002
0.007 0.001 0.003 0.989
0.001 0.992 0.002 0.004
0.002 0.002 0.993 0.004
0.004 0.001 0.993 0.003
0.001 0.276 0.002 0.722
```

#### DEG list format

Table contains gene ID's and log2FC and p-adj values for each RNA-seq experiment separated by tabulation. The number of experiments can be any. 

#### DEG list example

 ```
AT1G01200 0.798468892 0.999963365 0.262868659 0.858116364
AT1G01210 0.24716958  0.999963365 0.299284788 0.795936156
AT1G01220 -0.130354152  0.999963365 0.243543128 0.811947466
AT1G01225 -0.398548781  0.999963365 0.080885801 0.980423396
AT1G01230 0.094237385 0.999963365 -0.401992625  0.372281038
 ```




#### iSensor panel object

This is an example of panel with two genes 

- *genes*

`  "AT1G04730" "AT1G05055" `

- *gene_metadata*

| Gene | Chromosome | Site start | Site end | Strand | Site | To TSS | Gene name | Gene full name | 
| ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- | ---- |
| AT1G04730 | 1 | 1331105 | 1331112 | 1 | TATCGGAA | 28 | CTF18 | CHROMOSOME TRANSMISSION FIDELITY 18 |
| AT1G05055 | 1 | 1451800 | 1451807 | 1 | TGTCGTGA | 922 | GTF2H2 | general transcription factor II H2 |

*Strand* takes the value `1` if the site is on the direct chain and `0` if on the reverse chain.

- *panel_metadata*

| Species | Promoter length | Motif model name | Panel type | Transcriptomes experiment info | Date Created |
| ---- | ---- |----|----|----|----|
| Arabidopsis thaliana | 1500 | ARF1 - MA0942.1 | cis-trans | Auxin 1h, auxin 4h | 2025-06-20 |

  


#### Usage examples

- *Arabidopsis thaliana*, cis panel.

` iSensorsCisTransPanelCreate(panel_name = 'cis_panel', 
      species = 'Arabidopsis thaliana', 
      promoters_set = 'At_TAIR10_promoters.fas',
        ppm = 'ARF1.txt',
          panel_type = 'cis', 
      trivial_names_file = 'at_trivial.txt')`

- *Arabidopsis thaliana*, UP cis-trans panel.

` iSensorsCisTransPanelCreate(panel_name = 'cis-trans_panel', 
                      species = 'Arabidopsis thaliana',
      promoters_set = 'At_TAIR10_promoters.fas',
                      ppm = 'ARF1.txt',
                      deg_list = 'auxin_degs.txt',
                      min_dataset_number = 1,
                      panel_type = 'UP',
                      transcriptomes_info = 'Auxin 1h, auxin 4h') `

