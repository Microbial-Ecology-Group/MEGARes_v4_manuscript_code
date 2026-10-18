#### Introduction
This file describes the analysis workflow of SNV recall and Precision, and subsequent calculation of F1 scores. Steps include ground-truth data preparation, calculation of True-Positive SNVs (TP), False Positive SNVs (FP), and False Negative SNVs (FN), and a supplemental analysis of (read) abundance distribution of TP and FP SNVs. Recall is calculated as TP/(TP+FN). Precision is calculated as TP/(TP+FP). Both recall and precision where calculated based on formulas provided in [(Andreu-Sanchez et al. 2021)]([https://doi.org/10.3389/fgene.2021.648229](https://doi.org/10.3389/fgene.2021.648229)).

#### Step 1: Data preparation
This step makes two "True Positive" SNP Keys. Each "key" is a binary matrix, wherein column names represent the each "input-fasta" used for mock data generation, row-names represent all ground-truth genes and SNVs, and the intersecting cells include 0 for absent and 1 for present. For this analysis, there were 5 input-fastas per key and 10 total inputs. 

To recap, each input was used for the creation of 3 mock-metagenomic samples of 5M, 25M, and 50M reads. Each "sample" included "R0#" to indicate the ground-truth genes and SNVs that should be present under perfect SNV/Gene capture. Two sets of ground-truth genes/SNVs were used to generate 5 input-fasta. While input-fastas generated with the same ground truth set have similar gene/SNV composition, they do not share the same profile of genes/SNVs. Similarly, the abundance distribution of these genes/SNVs was generated at random for each resulting sample. 

The .csv files listing information of ground-truth genes/SNVs present in each sample, and the SNV-keys associated with each ground-truth set, were used to construct both "True-positive" keys. Thus, they are instrumental for subsequent analysis. This process involved re-formatting some outputs for down-stream parsing, which is also included in the script. 

The binary matricies were generated with the R script "V2_generate_analysis_inputs.Rmd". 
For each input fasta (R01-R10), the csv file containing the ground truth genes was loaded and used to interitively construct each binary matrix. Thus, for each ground-truth gene set (n=2), seven input files were loaded: 
- the set ground truth SNV key (n=1)
- the set ground-truth gene key (n=1)
- each input_fasta ground-truth gene key (n=5) 

Two files were output for each set: 
- set ground-truth gene presence/absence matrix 
- set ground-truth SNV presence/absence matrix 

#### Step 2: Analyze results

The analysis of recall, precision, and abundance of TP and FP SNVs is performed with a script called "*V2_SNV_calling_performance_abundance_analysis.Rmd*". 

To refresh, SNV calling was performed on each non-host sequencing depth group (5M, 25M, 50M) independently, as well has for the samples with simulated host reads. For each depth group and for host-removed samles, SNV calling analytic matrices were generated for each of 3 minimum alignment lengths (45bp, 75bp, 100bp). Thus, this analysis loads 12 distinct analytic matrices, each with n=10 samples. Each matrix has a corresponding metadata file (n=12). 

Three additional files also need to be loaded for the final analysis: both ground-truth binary matrices, described above, and a list of gene-accessions to filter. The latter list is for specifically filtering gene-accessions that contain the "RequiresSNPConfirmation" designation. These accessions are removed from output as they were not considered for ground-truth genes/SNVs to be injected into mock data, and were not removed from input genomes used to generate mock data. Not analyzing genes with this designation was a decision made when the study was designed. However, these genes were present in the database.fasta file during SNV calling, and thus appear in these analytic matrices. Updated versions of the python script used to generate the analytic matrices include a step to filter SNVs called in these "RequiresSNPConfirmation" accessions. 

This R script includes 3 functions to prepare-inputs for parsing (prepare_snp_data), tally SNVs meeting different criteria (get_snp_counts3), and calculate performance metrics (recall, precision) for subsequent plotting and analysis. The first two functions are run 9 times, one for each sequencing depth x alignment length combination. The outputs of are then combined into a single data frame, and the final function is run to calculate recall and precision metrics. The performance metrics (Recall, Precision) are then plotted. Analysis of host-removed SNV recall and precision is run in the same manner. In order to account for sample name differences in the host filtered data, a slightly different version of the function to calculate recall and precision is used. The only difference is that it splits the sample name differently than with the samples without simulated host. 

Results of recall and precision are reported as mean (within each sequencing depth group) plus standard deviation. 

#### Step 3: Additional metrics
Host removal: 
The inclusion of mock host reads allowed for subsequent analysis of host removal. After performing host removal using AMR++, host removal statistics were analyzed. To calculate percent host removal, total host-removed reads were divided by theoretical host reads present in the sample (20 million input reads). Prior to this calculation, host read counts were divided by two to get the average of removed reads between the forward and reverse sequencing files. Host removal across the 10 host samples was reported as mean + standard deviation. 

ARG reads retained after host removal: 
Estimation of ARG read depth (non SNV) was performed by running the 5M read depth samples (that did not have host reads) and the corresponding mock host samples (based on 5M read depth samples) through the complete AMR++ pipeline. The former was performed without host removal, as there were no simulated host reads in these samples. The latter was performed with host removal. Both sets of samples contain the same simulated ARG reads, with each non-host sample having a corresponding mock-host sample allowing for direct comparison. ARG read depths were calculated for each set of samples based on the corresponding AMR analytic matrix. The ARG read depths after host removal were divided by the ARG read depths from non host samples to calculate the total percentage of ARG reads retained after host removal. The result was reported as the mean and SD of percent ARG reads retained across all samples. 
