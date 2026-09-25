# iSensors <img src="man/figures/logo.png" align="right" width="140" alt="iSensors logo">

<img src="https://img.shields.io/badge/R-package-blue" alt="R package" height="20">

**[Tutorial](https://mironovalab.github.io/isensors-tutorial.html)**: a step-by-step
analysis with example data, from loading panels to plotting scores.

**iSensors** scores hormone signalling activity in single cells. Each *sensor* is a
panel of genes that report on one part of a signalling pathway, for example the
auxin response factors, auxin transport or cytokinin biosynthesis. For every cell,
iSensors averages the expression of the panel genes and stores the scores as new
assays of your Seurat object, ready for `FeaturePlot()`, `VlnPlot()` or
comparisons between cell types.

The package ships 676 ready-made panels for auxin and cytokinin in *Arabidopsis
thaliana*, tomato and 98 other plant species, and lets you build your own.

- [Installation](#installation)
- [Quick start](#quick-start)
- [How the scores are calculated](#how-the-scores-are-calculated)
- [Available panels](#available-panels)
- [Random and meta panels](#random-and-meta-panels)
- [Your own panels](#your-own-panels)
- [Tutorial](#tutorial)

## Installation

iSensors needs two Bioconductor packages, installed once:

```r
install.packages("BiocManager")
BiocManager::install(c("Biostrings", "universalmotif"))
```

Then install iSensors from GitHub:

```r
install.packages("devtools")
devtools::install_github("MironovaLab/iSensors")
```

To try a development version, name its branch, e.g.
`install_github("MironovaLab/iSensors", ref = "iSensors-dev")`.

If GitHub refuses the download (`Failed to install 'unknown package' from GitHub`),
create a personal access token under GitHub *Settings → Developer settings*, open
your R environment file with `usethis::edit_r_environ()` and add the line
`GITHUB_PAT=<your token>`.

## Quick start

```r
library(iSensors)
library(Seurat)

# Your log-normalised Seurat object (the tutorial provides example data)
seurat_obj <- readRDS("my_seurat_object.rds")

# Arabidopsis auxin panels
panelSet <- LoadSensors(setName = "ArabidopsisAuxin", species = "ATH", hormone = "aux")

# One score per panel and cell, added as the assay "iSensors_mean"
seurat_obj <- CalcSensors(seurat_obj, panelSet = panelSet, signals = "mean")

DefaultAssay(seurat_obj) <- "iSensors_mean"
FeaturePlot(seurat_obj, features = "ATH-aux-trans-ARF")
```

Each panel becomes one feature of the new assay, named like its panel file.

## How the scores are calculated

`CalcSensors()` reads the `data` layer of the default assay (change it with
`seurLayer`). For the `"mean"` signal the score of panel *P* in cell *j* is the
mean expression of the panel genes:

$$S_{Pj} = \frac{1}{|P|} \sum_{g \in P} x_{gj}$$

With Seurat's standard log-normalisation (`NormalizeData()`, scale factor 10,000),
where $x_{gj} = \ln\left(1 + 10^4 \, C_{gj} / \sum_k C_{kj}\right)$ for the raw
counts $C$, this is

$$S_{Pj} = \frac{1}{|P|} \sum_{g \in P} \ln\left(1 + \frac{C_{gj}}{\sum_k C_{kj}} \cdot 10^4\right)$$

Genes with zero variance across all cells are left out of *P*. A panel needs at
least 3 such genes to be scored; panels with fewer detected genes are skipped with
a warning that names them.

| `signals` | Score per cell | Assay |
|---|---|---|
| `"mean"` (default) | mean of the panel genes | `iSensors_mean` |
| `"median"` | median of the panel genes, ignoring zeros | `iSensors_median` |

The `"mean_normed"` and `"median_normed"` signals of earlier versions were removed
in 1.3.0.

`CalcSensors()` also accepts a genes × cells matrix instead of a Seurat object and
then returns an `iSensors` object with the scores in `$signals`.

## Available panels

Panel files are named `species-hormone-type-name`, for example
`ATH-aux-trans-ARF`. `LoadSensors()` filters on the first three parts:

| Filter | Values |
|---|---|
| `species` | species code, e.g. `ATH` (*Arabidopsis thaliana*), `SLY` (tomato), `OSA` (rice), `ZMA` (maize) |
| `hormone` | `aux` (auxin), `cyt` (cytokinin) |
| `type` | `trans`, `cis` or `reg` |

The three panel types:

- **trans**: genes of one part of the pathway, e.g. all ARF transcription factors
  or all auxin transporters.
- **cis**: genes whose promoters contain the binding site of a response factor,
  e.g. DR5 or IR8 sites recognised by ARF1.
- **reg**: cis genes that are also up- or down-regulated by the hormone in
  transcriptome experiments (suffix `-up` or `-down`).

| Species | Hormone | trans | cis | reg |
|---|---|---:|---:|---:|
| *Arabidopsis thaliana* (`ATH`) | auxin | 8 | 17 | 32 |
| *Arabidopsis thaliana* (`ATH`) | cytokinin | 4 | 10 | 20 |
| *Solanum lycopersicum* (`SLY`) | auxin | 6 | 18 | – |
| 98 other species | auxin | 6 each* | – | – |

\*ARF, IAA, PAT, receptors, synthesis and transport, where the gene family is
annotated for the species.

Arabidopsis trans panels: auxin `A-ARF`, `ARF`, `ConjugationDeconjugation`, `IAA`,
`PolarAuxinTransport`, `Receptors`, `Synthesis`, `Transport`; cytokinin `A-ARR`,
`B-ARR`, `Receptors`, `Synthesis`.

```r
ListSensorPanels()                          # all panel files
InspectSensorPanel("ATH-aux-trans-ARF.rda") # genes and metadata of one panel
```

<details>
<summary>All species codes</summary>

AAG *Anthoceros agrestis* · AAR *Aethionema arabicum* · ACERTR *Acer truncatum* ·
ACH *Actinidia chinensis* · ALY *Arabidopsis lyrata* · AMA *Avicennia marina* ·
AMHYB *Amaranthus hybridus* · AOX *Aquilegia oxysepala* · ARHY *Arachis hypogaea* ·
ATH *Arabidopsis thaliana* · ATR *Amborella trichopoda* · BCA *Brassica carinata* ·
BNA *Brassica napus* · BOL *Brassica oleracea* · BRA *Brassica rapa* ·
BVU *Beta vulgaris* · CAMSI *Camellia sinensis* · CAN *Capsicum annuum* ·
CANSAT *Cannabis sativa* · CAR *Cicer arietinum* · CAV *Corylus avellana* ·
CBR *Chara braunii* · CCAN *Coffea canephora* · CCL *Citrus clementina* ·
CDE *Ceratophyllum demersum* · CFA *Carpinus fangiana* · CHI *Cardamine hirsuta* ·
CIL *Carya illinoinensis* · CLA *Citrullus lanatus* · CME *Cucumis melo* ·
COL *Corchorus olitorius* · CPA *Carica papaya* · CQU *Chenopodium quinoa* ·
CRE *Chlamydomonas reinhardtii* · CRU *Capsella rubella* · CSA *Cucumis sativus* ·
DCA *Daucus carota* · DIN *Davidia involucrata* · DZI *Durio zibethinus* ·
ECA *Erigeron canadensis* · EGR *Eucalyptus grandis* · EGUT *Erythranthe guttata* ·
ESA *Eutrema salsugineum* · FAN *Fragaria × ananassa* · FVE *Fragaria vesca* ·
GHI *Gossypium hirsutum* · GMA *Glycine max* · GRA *Gossypium raimondii* ·
HAN *Helianthus annuus* · HMA *Hydrangea macrophylla* · LAL *Lupinus albus* ·
LJA *Lotus japonicus* · LONJA *Lonicera japonica* · LSA *Lactuca sativa* ·
MBI *Magnolia biondii* · MCO *Micromonas commoda* · MDO *Malus domestica* ·
MES *Manihot esculenta* · MPO *Marchantia polymorpha* · MTR *Medicago truncatula* ·
NNU *Nelumbo nucifera* · NTA *Nicotiana tabacum* · OEU *Olea europaea* ·
OSA *Oryza sativa* · PAX *Petunia axillaris* · PCO *Prasinoderma coloniale* ·
PGR *Punica granatum* · PPA *Physcomitrium patens* · PPE *Prunus persica* ·
PSA *Pisum sativum* · PSO *Papaver somniferum* · PTR *Populus trichocarpa* ·
PVU *Phaseolus vulgaris* · QLO *Quercus lobata* · RCH *Rosa chinensis* ·
RSI *Rhododendron simsii* · SAS *Striga asiatica* · SBO *Salvia bowleyana* ·
SBR *Salix brachista* · SCI *Simmondsia chinensis* · SED *Sechium edule* ·
SGI *Sequoiadendron giganteum* · SHI *Sapria himalayana* · SLY *Solanum lycopersicum* ·
SMO *Selaginella moellendorffii* · SPA *Schrenkiella parvula* · SPE *Solanum pennellii* ·
STU *Solanum tuberosum* · SUN *Selenicereus undatus* · TAR *Trochodendron aralioides* ·
TCA *Theobroma cacao* · THA *Tarenaya hassleriana* · TPR *Trifolium pratense* ·
TWI *Tripterygium wilfordii* · UGI *Utricularia gibba* · VMA *Vaccinium macrocarpon* ·
VMU *Vigna mungo* · VPL *Vanilla planifolia* · VVI *Vitis vinifera* · ZMA *Zea mays*

</details>

## Random and meta panels

By default `LoadSensors()` adds control panels: two panels of 200 and 500 randomly
chosen genes (`random1`, `random2`) and `majortrend`, the mean of all genes. Set
your own with `randomInfo`, or turn them off with `random = FALSE`:

```r
panelSet <- LoadSensors(setName = "ArabidopsisAuxin", species = "ATH", hormone = "aux",
                        randomInfo = list(n = 3, sizes = c(100, 200, 300), majortrend = TRUE))
```

Random panels are drawn when the scores are calculated; call `set.seed()` before
`CalcSensors()` for reproducible controls.

A **meta panel** combines the scores of existing panels with a function of your
choice:

```r
panelSet <- LoadSensors(
  setName = "ArabidopsisAuxin", species = "ATH", hormone = "aux",
  metaPanels = list(
    DR5_ARF1_ARF5 = list(srcPanels = c("ATH-aux-cis-DR5-ARF1", "ATH-aux-cis-DR5-ARF5-1"),
                         rule = mean)
  )
)
```

## Your own panels

Two functions create panels and save them as `iSensors/<panel_name>.rda` in the
working directory (the folder is created if needed). They also return the panel
invisibly, so `panel <- iSensorsTransPanelCreate(...)` keeps a copy in your
session. Load saved panels together with the default panels with
`LoadSensors(..., customPanels = TRUE)`.

A panel must contain at least 3 genes; both functions refuse to create a smaller
one, and `LoadSensors()` refuses to load one.

### From a gene list: `iSensorsTransPanelCreate()`

```r
iSensorsTransPanelCreate(panel_name = "ATH-aux-trans-myPanel",
                         gene_list = c("AT1G01010", "AT1G01030", "AT1G01040"),
                         species = "Arabidopsis thaliana",
                         panel_description = "Three genes of interest")
```

| Argument | Description |
|---|---|
| `panel_name` | Name of the panel and its file. Follow the `species-hormone-type-name` scheme to make the filters work. |
| `gene_list` | Gene IDs as a vector or a `.txt` file with one ID per line. |
| `species` | Species name, stored in the panel metadata. |
| `trivial_names_file` | Optional tab-separated file with gene ID, short name and full name. |
| `panel_description` | Free-text description. |

<details>
<summary>File formats</summary>

Gene list, one ID per line:

```
AT1G01010
AT1G01030
AT1G01040
```

Trivial names, tab-separated:

```
AT1G01010	NAC001	NAC domain containing protein 1
AT1G01020	ARV1	none
AT1G01030	NGA3	NGATHA3
```

</details>

### From a binding-site motif: `iSensorsCisTransPanelCreate()`

Scans promoters with a position weight matrix and keeps genes with a binding site
(`panel_type = "cis"`), optionally only those up- or down-regulated in
transcriptome experiments (`panel_type = "UP"` or `"DOWN"`).

```r
# cis panel
iSensorsCisTransPanelCreate(panel_name = "ATH-aux-cis-myMotif",
                            species = "Arabidopsis thaliana",
                            promoters_set = "At_TAIR10_promoters.fas",
                            ppm = "ARF1.txt",
                            panel_type = "cis",
                            trivial_names_file = "at_trivial.txt")

# genes with the site that are up-regulated in at least one experiment
iSensorsCisTransPanelCreate(panel_name = "ATH-aux-reg-myMotif-up",
                            species = "Arabidopsis thaliana",
                            promoters_set = "At_TAIR10_promoters.fas",
                            ppm = "ARF1.txt",
                            deg_list = "auxin_degs.txt",
                            min_dataset_number = 1,
                            panel_type = "UP",
                            transcriptomes_info = "Auxin 1h, auxin 4h")
```

| Argument | Description |
|---|---|
| `panel_name`, `species` | As above. |
| `promoters_set` | FASTA file of promoter sequences. |
| `ppm` | Position probability matrix in Homer format. |
| `panel_type` | `"cis"`, `"UP"` or `"DOWN"`. |
| `deg_list` | Table of log2 fold changes and adjusted p-values per experiment (for `"UP"`/`"DOWN"`). |
| `min_dataset_number` | Number of experiments in which a gene must be differentially expressed (adjusted p < 0.05). |
| `trivial_names_file` | Optional, as above. |
| `transcriptomes_info` | Free-text description of the experiments. |

<details>
<summary>File formats</summary>

**Promoters (FASTA).** One sequence per promoter. The header holds, separated by
`-`: gene ID, strand (`1` forward, `0` reverse), chromosome, promoter start and
end coordinates.

```
>AT1G01010-1-1-2130-3630
AACATTTCAAACCACTTGTTCTCTTTTATGTTTTGGTAAGAGCTATCTTC...
```

**Position probability matrix (Homer).** A header line starting with `>`, then one
row per position with the tab-separated probabilities of A, C, G and T:

```
>ARF1 - MA0942.1
0.312	0.207	0.22	0.262
0.008	0.088	0.004	0.901
0.002	0.002	0.993	0.002
0.007	0.001	0.003	0.989
0.001	0.992	0.002	0.004
0.002	0.002	0.993	0.004
0.004	0.001	0.993	0.003
0.001	0.276	0.002	0.722
```

**Differentially expressed genes.** Tab-separated: gene ID, then log2 fold change
and adjusted p-value for each experiment (any number of experiments).

```
AT1G01200	0.798	0.999	0.263	0.858
AT1G01210	0.247	0.999	0.299	0.796
```

</details>

### Panel structure

A panel is a list of class `GenePanel` with:

- `genes`: the gene IDs;
- `genes_metadata`: a table with gene ID, short and full name (cis panels add
  chromosome, site position, strand, site sequence and distance to the TSS);
- `panel_metadata`: species, panel type, description and creation date (cis panels
  add promoter length, motif name and transcriptome information).

## Tutorial

The [iSensors tutorial](https://mironovalab.github.io/isensors-tutorial.html) walks
through a full analysis of an example dataset, from installing the package and
loading panels to calculating and plotting scores. Its scripts are in the
[iSensors-supplementary](https://github.com/MironovaLab/iSensors-supplementary/tree/main/Tutorial-iSensors)
repository.

## Authors

Maxim Rybakov, Elena Zemlyanskaya, Vladislav Dolgikh and Victoria Mironova
([MironovaLab](https://github.com/MironovaLab)). Released under the MIT licence.
