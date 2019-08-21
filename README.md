## GUST - Genes Under Selection in Tumors
GUST predicts oncogenes (OGs), tumor suppressor gens (TSGs) and passenger genes (PGs) using somatic SNVs detected in a collection of tumors. It can classify single genes from targetted sequencing or multiple genes from whole-exome sequencing. 

## Install. 
````
# Be patient. It may take some time.
 devtools::install_github('liliulab/gust')
````

## Usage
```` 
library(gust)

# to make predictions. 
# A sample input can be downloaded from the "examples" folder. 
 gust(input.file.name='./examples/TCGA.ACC.mutect.somatic.maf.gz', output.folder='./examples/', output.prefix='TCGA.ACC');
# Compare outputs from the above command with files in the "examples" folder.
 
# to plot distribution
 m=plot.this.one('TP53', './examples/', 'TCGA.ACC') 
````

The input file needs to be in the VCF format. Required fields include Tumor_Sample_Barcode, Chromosome, Start_Position, dbSNP_RS, Reference_Allele, Tumor_Seq_Allele2, FILTER, One_Consequence, Hugo_Symbol, Gene, Feature, ENSP, HGVSc, HGVSp_Short, Amino_acids, Codons, ENSP, RefSeq, and Entrez_Gene_Id. We recommend using SnpEff to annotate somatic variants, which will generate a VCF file containing the required fields.
Additional usages can be found via using "?gust" after installing the R gust package.

## Precomputed results
GUST classifications of genes in 33 cancer types using the TCGA whole-exome sequencing data are available at https://www.liliulab.com/gust.

## Reference
Manuscript is available at BioRxiv and currently under review.

Chandrashekar P, Ahmadinejad N, Wang J, Sekulic A, Egan JB, Asmann YW, Mayley C, Liu L (2019) Contextual Classifications of Cancer Driver Genes. bioRxiv 715508; doi: https://doi.org/10.1101/715508 . 

## Contributors
Algorithm of GUST was developed by Li Liu. Database was developed by Pramod Chandrashekar. Please contact liliu at asu.edu for any questions or suggestions.
