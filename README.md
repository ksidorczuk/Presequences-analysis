# Presequences-analysis

## Scripts to run analyses:

### _targets.R
This targets pipepline manages all the datasets used in the project: sequence and annotation files, extraction of targeting signals or transmembrane domains, annotation processing, homology reduction, etc.

### 01-get_sequence_numbers.R
This script summarizes information extracted from raw datasets using number of extracted presequences/regions, number of those containing nonstandard amino acids, shorter than 10 and suitable for further analyses (only standard amino acids, length >= 10). 

### 02-run_predictions.R
This script uses datasets obtained in a targets pipeline and runs prediction of antimicrobial properties for peptides from each dataset. 


### 03-1-run_ngram_analyses_amp.R 
This script runs all analyses of motifs for presequences and AMPs. Each type of considered n-gram is analysed separately and report based on **ngram_analysis_amp.Rmd** is generated.

### 03-2-run_ngram_analyses_enc.R 
This script runs all analyses of motifs for presequences using amino acid encodings. Each type of considered n-gram is analysed separately and report based on **ngram_analysis_encoding.Rmd** is generated.

### 03-3-run_ngram_analyses_taxonomy.R 
This script runs all analyses of motifs for presequences divided into main taxonomic groups. Each type of considered n-gram is analysed separately and report based on **ngram_analysis_taxonomy.Rmd** is generated.

### 03-4-run_ngram_analyses_counts.R 
This script runs all analyses of motifs occurring more than once in a dataset. Each group of considered n-gram (e.g. bigrams with gaps, bigrams without gaps) is analysed separately and report based on **ngram_analysis_counts.Rmd** is generated.

### 03-5-run_ngram_analyses_amp_enc.R 
This script runs all analyses of motifs for AMPs using amino acid encodings. Each type of considered n-gram is analysed separately and report based on **ngram_analysis_amp_encoding.Rmd** is generated.

### 03-6-run_ngram_analyses_taxonomy_SP.R 
This script runs all analyses of motifs for SPs divided into taxonomic groups using different levels. Each type of considered n-gram is analysed separately and report based on **ngram_analysis_taxonomy_SP.Rmd** is generated.

### 04-analyse_motif_results.R
This script performs analysis of manually curated results obtained from all motif analyses and generates all plots with motif frequencies.

### 05-statistical_analysis.R
This script compares motif distribution between taxonomic groups for each analysed dataset using Mann-Whitney test with Benjamini-Hochberg correction. 

### 06-organisms_motifs_correlation.R
This script was used to determine if there is a correlation between number of frequently occurring motifs and the number of organisms in given taxonomic group. It requires data files generated in **04-analyse_motif_results.R**

### 07-proteins_imported_via_ES.R
This script contains analyses of small group of proteins that possess internal targeting peptides that allow their import via ensomembrane system. 



## Reports that summarise results:

### Prediction_results.Rmd
This report summarises analyses of antimicrobial activity prediction in all studied datasets.

### clustering.Rmd
This report performs analyses of amino acid composition, physicochemical properties and n-gram composition of each dataset. Based on these features clustering of sequences is performed using PCA, t-SNE and UMAP.

### clustering_selecting_AMPs.Rmd
This report was used to extract sequences of AMPs that clustered most closely with transit peptides.

### taxonomic_representation.Rmd
This report summarises taxonomic representation in datasets downloaded from UniProt.

