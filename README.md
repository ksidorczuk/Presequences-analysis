# Presequences-analysis

### _targets.R
This targets pipepline manages all the datasets used in the project: sequence and annotation files, extraction of targeting signals or transmembrane domains, annotation processing, homology reduction, etc.

### 01-get_sequence_numbers.R
This script summarizes information extracted from raw datasets using number of extracted presequences/regions, number of those containing nonstandard amino acids, shorter than 10 and suitable for further analyses (only standard amino acids, length >= 10). 

### 02-run_predictions.R
