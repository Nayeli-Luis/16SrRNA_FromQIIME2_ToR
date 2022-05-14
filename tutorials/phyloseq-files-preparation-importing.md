Files preparation and importing to `phyloseq`
================
Luis-Vargas, Maira N.

email: <nayeli.luis@ciencias.unam.mx>

## What is `phyloseq`?

`phyloseq` is a R package which brings tools in order to analyse
microbial communities from many types of data. The Github of `phyloseq`
([click here](https://joey711.github.io/phyloseq/), has many tutorials
for make several analysis and data visualization. But, in my opinion, is
unclear how the data must be imported, in this tutorial I hope to show
an easy way to import data for using `phyloseq`.

## First of all

We need four files, three of these files are insde `.qza` directories
from QIIME2, the last one is your \`sample-metadata.csv’ file. The next
table shows the directories where the files are, you just have to unzip
them.

From these file, we will create others, which are the ones that will be
import into `phyloseq`. These are: `feature-table.tsv` and
`otu_tax_matrix.csv`.

## The preparation

### `feature-table.tsv`

To get the `feature-table.tsv` file we need to transform the `biom`
file. If you want to know what a biom file is, yo can go to the [Biom
page](https://biom-format.org). In this case we have to transform the
`biom` file in a `.tsv` file. For that you’ll active your QIIME2 and
copy-paste the next comand.

`biom convert -i feature-table.biom -o feature-table.tsv --to-tsv`

The `.tsv` file looks like this:

``` r
library(tidyverse)
library(magrittr)

asv_table <- read.table("../datasets/feature-table.tsv", sep = "\t")
head(asv_table)
```

    ##                                 V1  V2 V3 V4    V5    V6   V7 V8   V9  V10  V11
    ## 1 7316c3832ec38cd0495800fdb90ee8b6   0  0  0 10218  8097 3945  0 3752 7435 6179
    ## 2 45500d433aca766bd01414417f075e13   0  0  0 28819 88824 1058  0  110  268  214
    ## 3 0c4044221ce447bccf82e751ba1d0d8c   0  0  0 15193 29690 2193 34  378 1055  860
    ## 4 2d2e697b25bebf9a20b02c01c72c9bb2   0  0  0  8815  7727  704 96  436  646  694
    ## 5 b2ff8f86352843d0d37e34277c937bc4   0  0  0    98   268  834  0    0    0    0
    ## 6 a33fdfe70b0d55aedfa9361dc98b79e1 140  0  0 12161 17431 2061  0  140  449  299
    ##    V12  V13  V14  V15  V16  V17   V18  V19  V20  V21  V22   V23  V24  V25   V26
    ## 1 8332 5416  368 1307  222  284 89584 7641 5701 8752 5582  4282 3848 4251  7463
    ## 2  269  184 1432  850  392    0  8200  557  101    0    0     0    0    0     0
    ## 3 1159  678 3690 2425 1190  197 17653 1467  307  251  198   836  262  414  1357
    ## 4  716  597 1647 1338  996   88  6741  583 1646 2334 4139 11153 4866 6056 39285
    ## 5    0    0    0    0    0   82    90    0    0    0    0     0    0    0     0
    ## 6  370  221  568  587  254 1497 18600 1228  143  196  279     0    0    0     0
    ##     V27  V28 V29 V30 V31  V32  V33  V34   V35  V36  V37   V38   V39   V40 V41
    ## 1  2119 1758   0   0   0  614  725  809  2734 2586  924  2575  2689  1866   0
    ## 2     0    0   0   0   0 1955 1067    0  9261 5890 2307   625   618   457   0
    ## 3   382  252   0   0 111 5463 6094 7506 17343 9174 4597   121   151   150  65
    ## 4 12309 9948   0   0   0  106  169  135   841  547  222   254   282   242   0
    ## 5     0    0   0   0   5    4    0    2  1193  623  276 28332 30638 28400   0
    ## 6     0    0   0 120  60 5202 5466 6477  6158 2681  959   131   187    88   0
    ##     V42  V43  V44  V45  V46 V47 V48 V49
    ## 1  2851 1060 1772 1061 1518   0   0 171
    ## 2   835  353    0    0    0   0   0   0
    ## 3   166  145  499  270  516   0  96 187
    ## 4   313   24 1209  620  515   0   0  62
    ## 5 22611 4576    0    0    0   0 176 464
    ## 6   110  161    0    0    0   0  49 130

The first column is the FeatureID (ASV) and the next ones are the
samples. So, the next step is put column names and we can do this
importing the `sample-metadata.csv` file, passing the names of the
samples as columnames.

``` r
metadata <- read.csv("../datasets/sample-metadata.csv", row.names = 1)
sample.names <- rownames(metadata)

colnames(asv_table) <- c("Feature.ID", sample.names)

head(asv_table)
```

    ##                         Feature.ID  A1 A2 A3    A4    A5   A6 A7   A8   A9  A10
    ## 1 7316c3832ec38cd0495800fdb90ee8b6   0  0  0 10218  8097 3945  0 3752 7435 6179
    ## 2 45500d433aca766bd01414417f075e13   0  0  0 28819 88824 1058  0  110  268  214
    ## 3 0c4044221ce447bccf82e751ba1d0d8c   0  0  0 15193 29690 2193 34  378 1055  860
    ## 4 2d2e697b25bebf9a20b02c01c72c9bb2   0  0  0  8815  7727  704 96  436  646  694
    ## 5 b2ff8f86352843d0d37e34277c937bc4   0  0  0    98   268  834  0    0    0    0
    ## 6 a33fdfe70b0d55aedfa9361dc98b79e1 140  0  0 12161 17431 2061  0  140  449  299
    ##    A11  A12  A13  A14  A15  A16   A17  A18  A19  A20  A21   A22  A23  A24    B1
    ## 1 8332 5416  368 1307  222  284 89584 7641 5701 8752 5582  4282 3848 4251  7463
    ## 2  269  184 1432  850  392    0  8200  557  101    0    0     0    0    0     0
    ## 3 1159  678 3690 2425 1190  197 17653 1467  307  251  198   836  262  414  1357
    ## 4  716  597 1647 1338  996   88  6741  583 1646 2334 4139 11153 4866 6056 39285
    ## 5    0    0    0    0    0   82    90    0    0    0    0     0    0    0     0
    ## 6  370  221  568  587  254 1497 18600 1228  143  196  279     0    0    0     0
    ##      B2   B3 B4  B5  B6   B7   B8   B9   B10  B11  B12   B13   B14   B15 B16
    ## 1  2119 1758  0   0   0  614  725  809  2734 2586  924  2575  2689  1866   0
    ## 2     0    0  0   0   0 1955 1067    0  9261 5890 2307   625   618   457   0
    ## 3   382  252  0   0 111 5463 6094 7506 17343 9174 4597   121   151   150  65
    ## 4 12309 9948  0   0   0  106  169  135   841  547  222   254   282   242   0
    ## 5     0    0  0   0   5    4    0    2  1193  623  276 28332 30638 28400   0
    ## 6     0    0  0 120  60 5202 5466 6477  6158 2681  959   131   187    88   0
    ##     B17  B18  B19  B20  B21 B22 B23 B24
    ## 1  2851 1060 1772 1061 1518   0   0 171
    ## 2   835  353    0    0    0   0   0   0
    ## 3   166  145  499  270  516   0  96 187
    ## 4   313   24 1209  620  515   0   0  62
    ## 5 22611 4576    0    0    0   0 176 464
    ## 6   110  161    0    0    0   0  49 130

Our first file (`feature-table.tsv`) is ready!

### `otu_tax_matrix.csv`.

This file has the OTUIDs, the taxonomy and the samples, and the taxonomy
should be separated. I mean, each taxonomic hierarchy must be in a
column. To get this file we need the `taxonomy.tsv` file and the
`feature-table.tsv`.

First, the `taxonomy.tsv` has three columns: `Feature.ID`, `Taxon` and
`Confidence`.

``` r
tax <- read.table("../datasets/taxonomy.tsv", sep = '\t', header = T)
colnames(tax)
```

    ## [1] "Feature.ID" "Taxon"      "Confidence"

The first step is delete the column `Confidence`. Next, we’ll rename the
column `Feature.ID` in `OTUID`, because we will merge the datasets by
this column. And finally, we’ll separate the column `Taxon` in different
columns.

``` r
tax %<>% 
  select(-Confidence) %>%
  separate(Taxon, c("Kingdom", "Phylum"), sep = "; p__") %>%
  separate(Phylum, c("Phylum", "Class"), sep = "; c__") %>%
  separate(Class, c("Class", "Order"), sep = "; o__") %>%
  separate(Order, c("Order", "Family"), sep = "; f__") %>%
  separate(Family, c("Family", "Genus"), sep = "; g__") %>%
  separate(Genus, c("Genus", "Specie"), sep = "; s__") %>%
  mutate_all(na_if,"")

tax$Kingdom <- sub (pattern = "d__", replacement = "", 
                       tax$Kingdom)

head(tax)
```

    ##                         Feature.ID  Kingdom           Phylum
    ## 1 7316c3832ec38cd0495800fdb90ee8b6 Bacteria   Proteobacteria
    ## 2 45500d433aca766bd01414417f075e13 Bacteria Actinobacteriota
    ## 3 0c4044221ce447bccf82e751ba1d0d8c Bacteria Actinobacteriota
    ## 4 2d2e697b25bebf9a20b02c01c72c9bb2 Bacteria   Proteobacteria
    ## 5 b2ff8f86352843d0d37e34277c937bc4 Bacteria   Proteobacteria
    ## 6 a33fdfe70b0d55aedfa9361dc98b79e1 Bacteria Actinobacteriota
    ##                 Class             Order             Family      Genus
    ## 1 Gammaproteobacteria         Ga0077536          Ga0077536  Ga0077536
    ## 2      Actinobacteria Pseudonocardiales Pseudonocardiaceae Crossiella
    ## 3      Actinobacteria Pseudonocardiales Pseudonocardiaceae Crossiella
    ## 4 Gammaproteobacteria   Nitrosococcales   Nitrosococcaceae    wb1-P19
    ## 5 Gammaproteobacteria   Nitrosococcales   Nitrosococcaceae    wb1-P19
    ## 6      Actinobacteria Pseudonocardiales Pseudonocardiaceae Crossiella
    ##                 Specie
    ## 1 uncultured_bacterium
    ## 2                 <NA>
    ## 3 uncultured_bacterium
    ## 4 uncultured_bacterium
    ## 5 uncultured_bacterium
    ## 6                 <NA>

And, we’ll merge the datasets. In order to have all de taxonomy with the
number of seqs by ASV and by sample.

``` r
asv_matrix <- merge(x = asv_table, y = tax, by = "Feature.ID")
head(asv_matrix)
```

    ##                         Feature.ID A1 A2 A3 A4 A5 A6 A7 A8 A9 A10 A11 A12 A13
    ## 1 00036a009561f213d4ce34976272a74e  0  0  0  0 14  0  0  0  0   0   0   0   0
    ## 2 000492f90810f8c40399573f9b78fbf4  0  0  0 28 26  0  0  0  0   0   0   0   0
    ## 3 0006d142c7d8d12f6a25328b73d993b3  0  0  0  8  0  0  0  0  0   0   0   0   0
    ## 4 000c2c548d916008a072ecb085317efb  0  0  0  0 23  0  0  0  0   0   0   0   0
    ## 5 000fd226d21abb5d0c38bfb05f7e9fa5  0  0  0  0  0  0  0  0  0   0   0   0   0
    ## 6 00130d1108179cf86d8ef16c2858f706  0  0  0  0  0  0  0  0  0   0   0   0   0
    ##   A14 A15 A16 A17 A18 A19 A20 A21 A22 A23 A24 B1 B2 B3 B4 B5 B6 B7 B8 B9 B10
    ## 1   0   0   0  13   0   0   0   0   0   0   0  0  0  0  0  0  0  0  0  0   0
    ## 2   0   0   0   0   0   0   0   0   0   0   0  0  0  0  0  0  0  0  0  0   0
    ## 3   0   0   0   0   0   0   0   0   0   0   0  0  0  0  0  0  0  0  0  0   0
    ## 4   0   0   0   0   0   0   0   0   0   0   0  0  0  0  0  0  0  0  0  0   0
    ## 5   0   0   0   0   0   0   0   0   0   0   0  0  0  0  0  0  9  0  0  0   0
    ## 6   0   0   0   0   0   0   0   0   0   0   0  0  0  0  0  0  0  0  0  0   0
    ##   B11 B12 B13 B14 B15 B16 B17 B18 B19 B20 B21 B22 B23 B24  Kingdom
    ## 1   0   0   0   0   0   0   0   0   0   0   0   0   0   0 Bacteria
    ## 2   0   0   0   0   0   0   0   0   0   0   0   0   0   0 Bacteria
    ## 3   0   0   0   0   0   0   0   0   0   0   0   0   0   0 Bacteria
    ## 4   0   0   0   0   0   0   0   0   0   0   0   0   0   0 Bacteria
    ## 5   0   0   0   0   0   0   0   0   0   0   0   0   0   0 Bacteria
    ## 6   0   0   0   0   0   0   0  18   0   0   0   0   0   0 Bacteria
    ##             Phylum               Class                 Order
    ## 1 Desulfobacterota          uncultured            uncultured
    ## 2            NB1-j               NB1-j                 NB1-j
    ## 3      Myxococcota         bacteriap25           bacteriap25
    ## 4   Proteobacteria Alphaproteobacteria   Paracaedibacterales
    ## 5   Proteobacteria Alphaproteobacteria      Rhodospirillales
    ## 6       Firmicutes             Bacilli Thermoactinomycetales
    ##                   Family               Genus                      Specie
    ## 1             uncultured          uncultured                  metagenome
    ## 2                  NB1-j               NB1-j             uncultured_soil
    ## 3            bacteriap25         bacteriap25            uncultured_delta
    ## 4   Paracaedibacteraceae Candidatus_Captivus        uncultured_bacterium
    ## 5             uncultured          uncultured uncultured_Rhodospirillales
    ## 6 Thermoactinomycetaceae          Planifilum        uncultured_bacterium

## Importing datasets to `phyloseq`

Finally, we can import our data to `phyloseq`:

``` r
library(phyloseq)

# Column Feature.ID as rownames in the asv_table
rownames(asv_table) <- asv_table$Feature.ID
asv_table %<>% select(-Feature.ID)

# Column Feature.ID as rownames in the df taxonomy
rownames(tax) <- tax$Feature.ID
tax %<>%
  select(-Feature.ID) %>%
   as.matrix(.) # Convert to a matrix


ASV <- otu_table(asv_table, taxa_are_rows = TRUE)
TAX <- tax_table(tax)
metadata <- sample_data(metadata)
phy_tree <- read_tree("../datasets/tree.nwk")

phyloseq_object <- phyloseq(ASV, TAX, metadata, phy_tree)
phyloseq_object
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 16536 taxa and 48 samples ]
    ## sample_data() Sample Data:       [ 48 samples by 2 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 16536 taxa by 7 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 16536 tips and 16511 internal nodes ]

I recommend you to save this physeq object as a RDS

``` r
ps <- saveRDS(phyloseq_object, "../datasets/phyloseq_object.RDS")
```
