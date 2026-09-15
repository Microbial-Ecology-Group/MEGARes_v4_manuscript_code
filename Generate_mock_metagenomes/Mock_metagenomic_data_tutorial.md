### Introduction
This a tutorial for making mock-metagenomic data using MEGARes v4 and the bacterial and archaeal genomes (n=19) from the Zymo gut mock (Cat: D6331). The result of this process is 30 mock metagenomes, with n=10 each sequenced to 5M, 25M, and 50M reads per sample, referred to throughout this tutorial as "depth groups". The methodology employed here creates mock paired-end samples with a high-degree of randomness. Within each depth group, n=5 inputs share similar antimicrobial resistance gene and ground-truth SNV composition. However, no two inputs have the same exact composition, and the relative abundance of each microbial genome in each mock-metagenome was selected at random. 

The workflow below was designed to replicate highly-random conditions, but could be modified to create sets of mock-metagenomes with greater similarity, particularly for simulating results from highly similar metagenomes, or to test hypotheses related to gene recall across metagenomes with slight deviations in abundance. Lastly, this workflow is run in steps, rather than as a complete pipeline. However, subsequent work could be done to create a pipeline. 
### Step 1: Genome pre-prep: 
The fasta files were downloaded from: [https://s3.amazonaws.com/zymo-files/BioPool/D6331.refseq.zip](https://s3.amazonaws.com/zymo-files/BioPool/D6331.refseq.zip)
For the sake of this analysis, the two fungal genomes were excluded from use in the mock data. This is, in part, due to the fragmented nature of the fungal genomes. Some of the fragments in each fasta file were short (< 150 bp), which caused errors when generating mock reads using insilico seq. Then, blastn was used to align the complete MEGARes database to each of the 19 genomes. MEGARes v4 was first converted to a binary blast database, then a loop was used to align each bacterial genome to the MEGARes database, as shown below: 

```
# blastn (2.13.0-gcc-8.2.0)
# Make a blast database of megares v4
makeblastdb -in megares_v4_database.fasta -dbtype nucl -out blast_mrv4_db

# Align each genome to database in a loop
for file in *.fasta; do blastn -query "$file" -db blast_mrv4_db -out "${file%.fasta}_results.txt" -outfmt 10 done
```

Then, a custom r script was used to filter the raw blast output in order to identify the "best supported" resistance gene at each unique location. The function "filter_alignment()" was used to filter raw alignment output as follows: 
- For alignments with the same start position, keep the aligned gene with the highest percent identity
- For alignments with the same stop position, keep the aligned gene with the highest percent identity 
- After those initial filtering steps, only keep alignments with ≥ 75% identity 
- Filter any aligned genes with the designator "RequiresSNPConfirmation" - associated with the MEGARes database. These genes were also excluded from downstream analysis. 
- Keep only aligned genes with ≥ 100bp alignment length. 
- Keep only aligned genes with a bitscore ≥ 100. 

The second function, "prepare_position_file", was used to subset the columns used later in the mock-data generation process. This kept the columns that included the fragment name (in each genome fasta file), the start position, and the stop position. This created one file per genome, called *genome_name_start_stop_pos.csv* . This position files are used in step 2 of the mock data generation process. Both R functions were run on each input genome independently. The primary goal was to identify all unique positions either containing an ARG, or with high-homology to an ARG sequence. These regions are subsequently replaced with a "ground-truth" set of ARGs later in the workflow. Regions of alignment-overlap were not considered in this initial filtering step, but are subsequently handled in later steps. Additionally, these identified regions are not used as ground truth later in analysis. The two functions used to accomplish this are included below:
```r
# R v.4.6.0
# dplyr v.1.2.1
library(dplyr)
filter_alignment <- function(input_csv){
  colnames(input_csv) <- c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
  filt_blast <- input_csv |> group_by(as.factor(qstart)) |>
    slice_max(order_by = pident, n = 1, with_ties = FALSE) |>
    ungroup()
  filt_blast <- filt_blast |> group_by(as.factor(qend)) |>
    slice_max(order_by = pident, n = 1, with_ties = FALSE) |>
    ungroup()
  filt_blast <- filt_blast[filt_blast$pident >= 75, ]
  filt_blast <- filt_blast |> filter(!str_detect(sseqid, "RequiresSNPConfirmation"))
  filt_blast <- filt_blast[filt_blast$length >= 100, ]
  filt_blast <- filt_blast[filt_blast$bitscore >= 100, ]

  return(filt_blast)
}

prepare_position_file <- function(input_filt){
  output_filt <- input_filt[,c(1,7,8)]
  colnames(output_filt) <- c("fragment","start","stop")
  return(output_filt)
}
```

This initial filtering step uncovered 1,050 alignment positions meeting the above criteria across 15/19 genomes. After accounting for alignment positions with overlap, this left 878 unique (combined) positions across the 15 genomes. Table 1 below includes the name of each bacterial genome, the number of unique "ARG" alignments in each genome after applying filter criteria (filt_pos), and the number of positions replaced in each genome after merging overlapping regions (merged_pos). 

| Genome                       | Filt_pos  | Merged_pos |
| ---------------------------- | --------- | ---------- |
| Akkermansia muciniphila      | 2         | 1          |
| Bacteroides fragilis         | 5         | 5          |
| Bifidobacterium adolescentis | 3         | 3          |
| Clostridioides difficile     | 12        | 11         |
| Clostridium perfringens      | 2         | 2          |
| Enterococcus faecalis        | 20        | 17         |
| E. coli B1109                | 171       | 142        |
| E. coli B3008                | 191       | 159        |
| E. coli B766                 | 176       | 148        |
| E. coli JM109                | 178       | 149        |
| E. coli b2207                | 196       | 162        |
| Faecalibacterium prausnitzii | 1         | 1          |
| Fusobacterium nucleatum      | 0         | NA         |
| Lactobacillus fermentum      | 0         | NA         |
| Methanobrevibacter smithii   | 0         | NA         |
| Prevotella corporis          | 1         | 1          |
| Roseburia hominis            | 0         | NA         |
| Salmonella enterica          | 93        | 75         |
| Veillonella rogosae          | 2         | 2          |
| **Total**                    | **1,053** | **878**    |

### Step 2: Prepare ground-truth ARG sequences and SNVs
**Step 2.1 - Generate random ARGs**
Used to generate a list of random ARGs from the MEGARes database. These genes will be used to generate true positive SNVs in a subsequent step. Use this syntax to do so: 
```
python3 Random_ARG_select.py database.fasta annotations.csv --seed 42 --n-per-class 40
```
For this analysis, two sets of input genes were created, titled "set 01" and "set 02", generated with seeds 300 and 301, respectively. This script selected up to 40 random gene accessions within each class of resistance, and resulted in 1072 randomly selected genes per set. 
This outputs several files: 
1) set_n_annotations.csv
2) set_n_gene_lengths.csv
3) set_n_genes.fasta

**Step 2.2 - Generate random SNVs**
This will be used to generate the list of "ground-truth" SNVs. Outputs are given in both 0 and 1 based SNP calls. Run the following script with this command, using the randomly selected set of genes produced above as the input: 
```
python3 SNP_injector_v1.py limited_database.fasta --seed 42
```
Random SNV injection was run with seed = 50 for set 1, and seed = 51 for set 2. The script randomly injects 1 to 10 SNVs in each accessions, and records both the zero- and one- based positions of each injected variant. The script outputs the following files: 
1) set_n_replacement_genes.fasta
2) set_n_zero_based_SNP_key.csv
3) set_n_one_based_SNP_key.csv

For downstream analysis, the "set_n" gene annotations and SNV keys, in tandem with the concatenated merged_csvs created in step 2, can be used to generate a "ground-truth" matrix, where columns are "sample"/replicate names, rows are either gene_accessions or SNV_accessions, and cells containing 0 = absent or 1 = present. This can be used for precision and accuracy analysis, if so desired. 

### Step 3: Prepare input files of mock-data generation 
Step 3 creates the input fasta files used for mock data generation. The steps below produce one mock input fasta per run. The steps below are run in the same directory, which should include the following files before running. 

```
generate_mock_inputs/
  ├── genome_name_1.fasta
  ├── genome_name_1_start_stop_pos.csv
  ├── genome_name_2.fasta
  ├── genome_name_2_start_stop_pos.csv
  ├── replacement_genes.fasta
  ├── Whole_genome_cleaner_v2.py
  └── run_all_genomes.sh

```

**Step 3.1 - Replace genes and pre-identified positions** This step uses run_all_genomes.sh based on Whole_genome_cleaner_v2.py. This command expected the genes used for replacement to be named "replacement_genes.fasta". Run it with this command: 
```bash
chmod +x run_all_genomes.sh
./run_all_genomes.sh
```
The output directory structure is like so: 
```
generate_mock_inputs/
 │  ├── genome_name_1.fasta
 │  ├── genome_name_1_start_stop_pos.csv
 │  ├── genome_name_2.fasta
 │  ├── genome_name_2_start_stop_pos.csv
 │  ├── replacement_genes.fasta
 │  ├── Whole_genome_cleaner_v2.py
 │  └── run_all_genomes.sh
 │ 
 └── replacement_results/
      ├── genome_name_1/
      │     ├── genome_name_1_insertions.csv
      │     ├── genome_name_1_replaced.fasta
      │     └── unused_genes_####.fasta
      ├── genome_name_2/
      │     ├── genome_name_2_insertions.csv
      │     ├── genome_name_2_replaced.fasta
      │     └── unused_genes_####.fasta
      └── unused_genomes/
	      ├── unused_genome_1_replaced.fasta
          └── unused_genome_2_replaced.fasta
```

The seed used for each of the 10 replicates is denoted in the table below: 

| Replicate | Gene Set | Seed |
| --------- | -------- | ---- |
| R01       | Set 1    | 201  |
| R02       | Set 1    | 202  |
| R03       | Set 1    | 203  |
| R04       | Set 1    | 204  |
| R05       | Set 1    | 205  |
| R06       | Set 2    | 206  |
| R07       | Set 2    | 207  |
| R08       | Set 2    | 208  |
| R09       | Set 2    | 209  |
| R10       | Set 2    | 210  |

**Step 3.2 - Combine modified output genomes into a single fasta file** This will only combine genomes that have genes replaced. The python script opens each sub_directory in the "replacement_results" directory, and looks for the string "replaced.fasta", so genomes that did not require replacement can be placed into an additional directory "unreplaced/", and names updated to "genome_name_replaced.fasta".  Then run:
```
python3 concatenate_fastas_v2.py /path/to/parent/ output_name.fasta
```

In case any genomes were printed with "newline" characters, run this script to ensure all genomes are re-linearized. Linearization is applied to sequences under each genome/fragment ID within the fasta file. 
```
python3 linearize.py file_to_fix.fasta file_to_fix_lin.fasta
```

 **Step 3.3 - Combine output position files into a single .csv file** This script requires pandas to be included in your python3 environment. Gene inputs and positions for each input fasta are concatenated into a single .csv file like so: 
 
```
python3 merge_csvs.py --root /dir/with/sub_dirs/ --output outname.csv
```

For this analysis, there should be 10 input fasta files in total. Thus, the above process for step 2 was repeated 10 times, each in a directory with a unique name. 5 of 10 inputs were each generated with the same "replacement_genes" generated in step 1, which was included in the final sample name. Replacement of "excised regions" was performed at random, with a different --set_seed used for each input. These each of the 10 fastas will be used to make 3 mock metagenomes of 3 different sequencing depths. This makes final input fasta for generating the mock metagenomes. 
##### Step 3 Directory setup 
Preparation for step 3.1 requires three input file types to be placed in the same directory: 
1) genome.fasta = the unmodified genome file for the taxa of interest, one per taxa. 
2) taxa_start_stop_pos.csv = the start and stop positions (aka "position file") of gene regions to be replaced, one position file per genome.fasta
3) replacement_genes.fasta = the "ground-truth" genes to replace the previously-identified regions in the position file

Additionally, two program files are required to be placed in this directory: 
1) Whole_genome_cleaner_v2.py
2) run_all_genomes.sh

With the proper inputs in place, modifications to run_all_genomes.sh can be made: 
- Line 5: current_gene_pool = "replacement_genes.fasta" - this can be changed to the name of the fasta file containing "ground-truth" genes you would like to add
- Line 6: output_dir="" - can be changed to whatever you would like the output results directory to be 
- Line 32 : can change the --seed argument to run with whichever seed is preferable. 

This bash script runs the "Whole_genome_cleaner_v2.py" in a loop, and applies it to n=x genomes you include in the processing directory. The outputs of each run are included in /replacement_results/genome_name/. Importantly, this script is designed to iterate through each unique position in genome specified the position file, remove the specified position, and replace the position in the genome with a "ground-truth gene". When two or more unique positions in the position file overlap, the length cumulative length of all overlapping positions is removed as one "excised region", and replaced with one ground-truth gene. Ground truth genes are selected at random from the "replacement_genes.fasta" file, without replacement. After the python script has finished working on one genome, it outputs three files in /replacement_results/genome_name/: 

1) genome_name_isertions.csv = a csv file containing pertinent information about where in the genome a position was excised, and what gene replaced it. 
2) genome_name_replaced.fasta = a fasta file containing the genome with ground-truth sequences. The genome used for mock metagenomic sample creation 
3) unused_genes_####.fasta = a fasta file of the unused replacement genes, where #### = the number of replacement genes left to be randomly selected from. This is by Whole_genome_cleaner_v2.py when iterating on the next genome. 

Step 3.2 combines the replaced genomes to one input fasta file. To prepare for step 2.2, it is helpful to have a directory called "unused_genomes/" that contain the genome.fasta files that you also want to include in the mock metagenomic dataset, but do not want to replace genes in. These will be used in the final "input.fasta" file. After running "run_all_genomes.sh", you can copy this "unused_genomes" directory into the replacement_results sub-directory. Unused genomes should include the statement "replaced" prior to ".fasta" in the name. 

#### Step 4: Generate Mock Metagenomic read files. 
This step was performed with [InSilicoSeq](https://insilicoseq.readthedocs.io/en/latest/iss/install.html), and is fairly straight-forward. InSilicoSeq can be downloaded as a conda environment, as below. 

```
conda create --name insilicoseq -c bioconda insilicoseq
conda avtivate insilicoseq
```

Each input.fasta file was used to create mock metagenomic datasets of three sequencing depths: 5M paired-end reads, 25M paired-end reads, and 50M paired-end reads. Each mock sample was also assigned a sample name "S001-S030". The code block below includes the syntax to generate samples at each sequencing depth, and was run once per input-fasta at each sequencing depth. 

```
# 5M reads 
iss generate -g input_fasta_files/Input_R##_lin.fasta --n_reads 10M --model NovaSeq --cpus 10 --seed ### --output mock_meta_10M/S###_rd010m_R## --store_mutations

# 25M reads
iss generate -g input_fasta_files/Input_R##_lin.fasta --n_reads 50M --model NovaSeq --cpus 100 --seed ### --output mock_meta_50M/S###_rd010m_R## --store_mutations

# 50M reads
iss generate -g input_fasta_files/Input_R##_lin.fasta --n_reads 100M --model NovaSeq --cpus 100 --seed ### --output mock_meta_100M/S###_rd010m_R## --store_mutations
```

--store_mutations creates a vcf file of "sequence errors" injected as part of the given InSilicoSeq error model. While InSilicoSeq can also compress mock read files as part of the ISS command, it appears to uses a gzip process. This can become problematic when simulating greater sequencing depths, as the time needed to gzip each individual file increases with file size. As an alternative, files can be output unzipped and be subsequently ziped with [pigz](https://zlib.net/pigz/). 

The table below specifies the seed used to generate each mock sample: 

| Input | SN_5M | Seed_5M | SN_25M | Seed_25M | SN_50M | Seed_50M |
| ----- | ----- | ------- | ------ | -------- | ------ | -------- |
| R01   | S001  | 123     | S011   | 133      | S021   | 143      |
| R02   | S002  | 124     | S012   | 134      | S022   | 144      |
| R03   | S003  | 125     | S013   | 135      | S023   | 145      |
| R04   | S004  | 126     | S014   | 136      | S024   | 146      |
| R05   | S005  | 127     | S015   | 137      | S025   | 147      |
| R06   | S006  | 128     | S016   | 138      | S026   | 148      |
| R07   | S007  | 129     | S017   | 139      | S027   | 149      |
| R08   | S008  | 130     | S018   | 140      | S028   | 150      |
| R09   | S009  | 131     | S019   | 141      | S029   | 151      |
| R10   | S010  | 132     | S020   | 142      | S030   | 152      |

#### Step 5: Align mock samples to MEGARes v.4
This step performs alignment and subsequent filtering using NGLess (1.5). Internally, NGLess uses bwa-mem to perform alignments, and samtools to filter alignments based on minimum criteria. The base filtering parameters (minimum alignment length = 45bp, minimum average nucleotide identity = 97%, uniquely mapped reads) were the minimum alignment input parameters necessary for metaSNV v2 to perform SNV calling.

Download [NGLess](https://ngless.readthedocs.io/en/latest/install.html) conda envrionment: 
```
conda create --name ngless -c bioconda ngless
```

Manually running NGLess requires several inputs to be created before starting the mapping/filtering process. 
- sample_filename.txt - a list of the unique file-names (1 per sample, no extension)
- NGLess_output - a directory for NGLess to output alignment and intermediate files.

Two NGLess scripts were used for this process: 
- map.ngl --> performs alignment to the complete MEGARes database, filters alignments based on specified parameters, and outputs alignment.bam files only containing those alignments meeting said parameters. 
```
# map.ngl
ngless "1.5"
import "parallel" version "1.1"
import "mocat" version "0.0"
import "samtools" version "0.0"

# Load inputs 
## Current sample 
current = run_for_all(readlines("/path/to/sample_filenames.txt"))
input = paired('path/to/forward/reads/' + current + '_R1.fastq.gz', 'path/to/reverse/reads/' + current + '_R2.fastq.gz')

fastaRefDb="path/to/megares_v4_database.fasta"

## Read pre-processing
input = preprocess(input) using |read|:
    read = substrim(read, min_quality=20)
    if len(read) < 45:
        discard

mapped = map(input, fafile= fastaRefDb, mode_all=True)
# only keep mappings with at least 97% identity and 45 bp length
mapped = select(mapped) using |mr|:
 mr = mr.filter(min_match_size=45, min_identity_pc=97, action={unmatch})

# only keep reads that mapped uniquely (discard multi-mappers)
mapped_unique = select(mapped, keep_if=[{mapped}, {unique}])
mapped_unique = samtools_sort(mapped_unique)
write(mapped_unique, ofile='NGLess_outputs/' + current + '.unique.sorted.bam')

```

- filter.ngl --> filters an alignment.bam file based on specified parameters. 
```
# filter.ngl
ngless "1.5"
import "parallel" version "1.1"
import "mocat" version "0.0"
import "samtools" version "0.0"

current = run_for_all(readlines("/path/to/sample_filenames.txt"))
input = samfile('path/to/NGLess_outputs/' + current + 'unique.sorted.bam')

# only keep mappings with at least 97% identity and 100 bp length
filtered = select(input) using |mr|:
 mr = mr.filter(min_match_size=100, min_identity_pc=97, action={unmatch})
 
# only keep reads that mapped uniquely (discard multi-mappers)
filtered_unique = select(filtered, keep_if=[{mapped}, {unique}])
filtered_unique = samtools_sort(filtered_unique)
write(filtered_unique, ofile='NGLess_outputs/' + current + '_100bp_97ani.unique.sorted.bam')
```

map.ngl was used on raw mock.fastq.gz samples, with the minimum base alignment length (45bp), for all three sequencing depths. filter.ngl was used to further filter the output bam files to only get reads with a minimum alignment length of 75bp and 100bp, respectively. Each unique combination of sequencing depth and minimum alignment length was processed in a separate directory. 

Both scripts run on a per-sample basis, and iterate to the next sample with each successive command. It can be run with a bash loop, as written below, one iteration per sample (10 in "-le 10" can be swapped with number of unique samples.). The number of iterations can be found with "ws -l sample_filenames.txt". Since the filter.ngl file is created at the start of a run, all specified parameters are included within the script, save for threads. The number of threads with which to run the command is "-j NN".

```
# Run syntax
i=1; while [ $i -le 10 ]; do ngless -j100 filter.ngl; i=$((i+1)); done
```

### Step 6: Call SNVs with metaSNV v2. 
SNVs are called using the main metaSNV python file, which outputs the files needed to construct the count matrix. Unlike with NGLess, if you prematurely create the output directory, then metaSNV will throw an error. The output directory is therefor created with the metaSNV.py command. There are other dependencies in the overall software, so this may be best managed as a conda environment. The only file that needs to be created manually is "bam_filepaths.txt", which should contain the absolute path to each bam file, one per file. This can also be done on the command line like so: 
```
## get all file paths in alignment directory
find $PWD > file_paths.txt
## remove first and last entry
sed -i '1d;$d' file.txt
```

Download [metaSNV v2](https://github.com/metasnv-tool/metaSNV/tree/master) conda environment: 
```
### create new environment 
conda create --name metaSNV -c bioconda -c conda-forge 'metasnv>=2.0.1'

### load metaSNV environment 
conda activate metaSNV
```

Run metaSNV:
```
# Run syntax
metaSNV.py [path/to/metaSNV/output_dir/] [path/to/bam_filepaths.txt] [path/to/MEGARes_db.fasta]

### Internal AMR++ defaults 
metaSNV output DIR         = test_results/SNV_analysis_output/
path to bam_filepaths.txt  = test_results/Alignment/Sam_files/bam_filepaths.txt path to MEGARes db.        = data/amr/megares_v4_database.fasta
```

### Step 7: Format SNV output for recall and precision analysis. 
This step uses a custom python script, clean_metaSNP_1.0.1.py. This script differs slightly from the second version, which is used in the actual AMR++ pipeline. The main difference is the formatting of SNV annotations, which contain an extra "accession" for easier parsing for recall and precision analysis. This script is run like so: 
```
## Clean outputs
python3 clean_metaSNP_1.0.1.py /path/to/raw/called/snp/files/ snpfile_name /path/to/files/in/run/order.txt
```

This outputs a SNV-level count matrix, and a corresponding annotation file. Only the count matrix is used for recall and precision analysis. With that, all inputs for supplemental analysis in R are generated. 
