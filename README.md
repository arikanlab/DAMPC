<h1 align="center">DAMPC: Dual Activity Microbial Peptides Catalog</h1>

If you use this repository, please cite:

Sağıroğlu Ö, Arıkan M, DAMPC: dual activity microbial peptides catalog, *bioRxiv*, 2025 [Link].

---

# Table of contents
- [Overview](#overview)
- [Step-by-Step Workflow](#step-by-step-workflow)
- [References](#references)

---

# Overview
**DAMPC** is a comprehensive catalog of microbial peptides predicted to have **dual activity**—both anticancer (ACP) and antimicrobial (AMP) properties.  

The DAMPC datasets are publicly available through the links below:

| Catalog   | Catalog link |
|-----------|-----------|
| DAMPCS - Peptides | [Link 1](https://example.com/data1) |
| DAMPC - Metadata | [Link 2](https://example.com/data2) |

DAMPC is also accessible as a searchable catalog at https://arikanlab.com/Home/DAMPC

---

# Step-by-Step Workflow
DAMPC was constructed by:

1) Collecting microbial **smORFs/short peptides** from public catalogs,  
2) Running multiple **ACP and AMP predictors** in a consensus approach, and 
3) Comparing the predictions against known **ACP/AMP databases**

---

## Requirements
- **conda**
- **git**
- **seqkit** (for sequence prep)

---

## Part 1 — Downloading and preparing smORFs

### Download Global Microbial smORFs Catalog v1.0
Go visit: https://gmsc.big-data-biology.org/downloads

From Nucleotide sequence (.fasta) section:
Download 100AA smORF catalog: GMSC10.100AA.fna.xz

#### Download DBsmORF Catalog
Go visit: http://104.154.134.205:3838/DBsmORF/

From Download section:
Download both RefSeq and HMP data. Be sure to select Protein Sequences in the Download Fields section.


### Filtering
For filtering, we used **Seqkit**

Install SeqKit via conda:
```bash
conda install -c bioconda seqkit
```
Since the catalogs contain methionine (M) as the initial amino acid and M is trimmed during the analyses, the sequence length was set to 11–51.

This command selects sequences with lengths between 11 and 51 amino acids, writes them in single lines, and saves the result to length_filtered.fasta:
```bash
seqkit seq -g -m 11 -M 51 input.fasta | seqkit seq -w 0 > length_filtered.fasta 
```
**Note:** SeqKit can also work directly with compressed files, both as input and output.

This command checks each sequence in length_filtered.fasta. If the first amino acid is M (Methionine), it is removed; otherwise, the sequence is kept as is. The output is saved to output.fasta:
```bash
awk '/^>/ {print; next} {if(substr($0,1,1)=="M") print substr($0,2); else print}' length_filtered.fasta > output.fasta
```
**Note:** Macrel already performs this trimming automatically, so for Macrel this step is not  necessary.

If you have a large FASTA file, you can split it into smaller parts with SeqKit to make the tools run more efficiently.
This command splits output.fasta into 5 equal parts and saves them into the folder parts/:
```bash
seqkit split2 -p 5 output.fasta -O parts/
```

## Part 2 — Running ACP and AMP prediction tools

### AntiCP2
#### Setup
```bash
conda create -n anticp2 python=3.7
conda activate anticp2

pip install pandas==1.3.5
pip install numpy==1.21.6
pip install scikit-learn==0.22.2
pip install pickle-mixin==1.0.2

git clone https://github.com/raghavagps/anticp2.git
cd anticp2
```

#### Usage
```bash
python3 anticp2.py -i input.fasta -o anticp2_output.csv -j 1 -t 0.5 -m 2 -w 10 -d 2
```

#### Notes

- -j : Job type {1: predict, 2: design, 3: scan} 
- -m : Model {1: ACP/AMP, 2: ACP/non-ACP} 
- -t : Threshold (0–1, default 0.5) 
- -w : Window length (scan mode only) 
- -d : Display {1: ACP only, 2: All peptides} 
- Output File: The Seq ID column contains the ">" symbol. This should be taken into account for further analyses. Positive results: AntiCP, Negative results: Non-AntiCP. File format is CSV.

---

### ACPred-BMF
#### Setup
```bash
conda create -n acpredbmf python=3.8
conda activate acpredbmf

pip install protobuf==3.19.0
pip install numpy==1.22.3
pip install pandas==1.3.4
pip install tensorflow==2.8.0
pip install scikit-learn==1.0.1
pip install keras==2.8.0

git clone https://github.com/RUC-MIALAB/ACPred-BMF.git
cd ACPred-BMF
```

#### Usage
```bash
bash main.sh alternative
```

#### Notes

- Alternative model: trained on ACPs/random peptides as positive/negative samples 
- Input: Place FASTA files of peptide sequences inside the data/ directory 
- Output File: There is no column for Seq ID. Positive results: 1, Negative results: 0. File format is CSV.

---

### Con_ACP
#### Setup
```bash
conda create -n conacp python=3.10
conda activate conacp

pip install numpy==1.22.3
pip install torch==2.7.1
pip install biopython

git clone https://github.com/bzlee-bio/con_ACP.git
cd con_ACP
```

#### Usage
```bash
python3 inf.py --input input.fasta --batch_size 20 --model_type ACP_Mixed_80 --device cpu --output output.txt
```

#### Notes

- model_type: ACP_Mixed_80 (default), ACP2_main, ACP2_alter, ACP500_ACP164, ACP500_ACP2710, LEE_Indep 
- device: cpu or gpu 
- Input: FASTA format file 
- Output File: There is no column for sequences. Positive results: ACP, Negative results: – . File format is TXT.

---

### AMP-Scanner-v2
#### Setup
```bash
git clone https://github.com/dan-veltri/amp-scanner-v2.git
cd amp-scanner-v2

conda env create -f environment_original_paper_model.yml
conda activate ascan2_orig
```

#### Usage
```bash
python amp_scanner_v2_predict_tf1.py \
  -fasta input.fasta \
  -model trained-models/OriginalPaper_081917_FULL_MODEL.h5
```

#### Notes

- Uses the original paper model provided in the repository 
- Input: FASTA file with peptide sequences 
- Output File: Positive results are shown as " AMP " and negative results are shown as " Non-AMP ".

---

### Macrel
#### Setup
```bash
conda create -n env_macrel -c bioconda macrel
conda activate env_macrel
```

#### Usage
```bash
cd /yourFile
macrel peptides -f input.fasta -o out_peptides --keep-negatives
```

#### Notes

- --keep-negatives: keeps negative predictions in the output file 
- Important: /yourFile is the directory where your input FASTA file is located 
- Output File: Positive results are shown as " True " and negative results are shown as " False ".

---

## Part 3 — Comparisons with ACP and AMP databases

Since Con_ACP result files do not contain sequences by default, with add_sequences.py, you can automatically add the corresponding sequences from a FASTA file into the result table based on matching IDs:

```bash
python add_sequences.py FASTA_FILE INPUT_FILE OUTPUT_FILE [options]

```
**Arguments** 

- FASTA_FILE – Path to the FASTA file containing sequences. 
- INPUT_FILE – Con_ACP result file (tab-delimited, contains at least an ID column). 
- OUTPUT_FILE – Path where the merged output will be written. 

**Options**

- --id-column (default: ID) – Column in the input file with sequence IDs. 
- --seq-column (default: Sequence) – Name of the new column to store sequences. 
- --sep (default: \t) – Field separator in the input/output file. Example: --sep , for CSV. 
- --chunksize (default: 100000) – Rows processed per chunk. 
- --fasta-id-split (default: " ") – Use only the part of FASTA ID before this character. Use "" to disable. 
- --encoding (default: utf-8) – File encoding. 
- --verbose – Print extra diagnostic info (chunk sizes, match counts, sample IDs). 

---

Since ACPred-BMF result files do not contain Seq IDs by default, with add_ids.py, you can automatically add the corresponding IDs from a FASTA file into the result table based on matching sequences:

```bash
python add_ids.py FASTA_FILE INPUT_FILE OUTPUT_FILE [options]

```
**Arguments** 

- FASTA_FILE – Path to the FASTA file containing sequences. 
- INPUT_FILE – ACPred_BMF result file (tab-delimited, contains at least an ID column). 
- OUTPUT_FILE – Path where the merged output will be written. 

**Options** 

- -seq-column name of sequence column (default: Sequence) 
- --id-column output ID column name (default: ID) 
- --sep field separator (, for CSV, "\t" for TSV) 
- --fasta-id-split " " split FASTA header at first space (use "" to disable) 
- --upper, --strip-ws normalize sequences before matching 
- --verbose – Print extra diagnostic info (chunk sizes, match counts, sample IDs). 

---

Use the to_fasta.py script to convert CSV, TSV, or TXT files into FASTA format by specifying which columns contain the sequence IDs and sequences:

```bash
# Header is on the 2nd row, using column names
python to_fasta.py input.tsv -o out.fasta --id-col ID --seq-col Sequence --header-row 2

# No header present, use 0-based column indices
python to_fasta.py input.txt -o out.fasta --id-col 0 --seq-col 1 --header-row 0

# Header is on the 1st row (default), delimiter auto-detected
python to_fasta.py input.csv -o out.fasta --id-col accession --seq-col peptide

# Manually specify delimiter (here: tab)
python to_fasta.py input.txt -o out.fasta --id-col 0 --seq-col 1 --header-row 0 --delimiter "\t"
```
**Explanation of key options** 
- --id-col: Column containing the sequence ID. Can be a column name (case-insensitive) if there is a header. Or a 0-based index if there is no header. 
- --seq-col: Column containing the sequence itself. Works the same way as --id-col. \
- --header-row: Tells the script where the header row is. 1 = header is in the first line (default). 2 = header is in the second line, etc. 0 = no header present (then you must use numeric column indices). 
- --delimiter: (Optional) Manually set the delimiter. Default: auto-detect (csv.Sniffer + fallback to common delimiters like comma, tab, semicolon, pipe, or whitespace). Example: --delimiter "\t"

---

Use the fasta_compare.py script to compare two FASTA files and write out the common sequences and the sequences unique to the second file:

```bash
python fasta_compare.py <first.fasta> <second.fasta>
```

---
# References
- Duan, Y., Santos-Júnior, C.D., Schmidt, T.S. et al. A catalog of small proteins from the global microbiome. Nat Commun 15, 7563 (2024). https://doi.org/10.1038/s41467-024-51894-6. Website: https://gmsc.big-data-biology.org/

- Sberro H, Fremin BJ, Zlitni S, Edfors F, Greenfield N, Snyder MP, Pavlopoulos GA, Kyrpides NC, Bhatt AS. Large-Scale Analyses of Human Microbiomes Reveal Thousands of Small, Novel Genes. Cell. 2019 Aug 22;178(5):1245-1259.e14. doi: 10.1016/j.cell.2019.07.016. Website: http://104.154.134.205:3838/DBsmORF/

- Durrant MG, Bhatt AS. Automated Prediction and Annotation of Small Open Reading Frames in Microbial Genomes. Cell Host Microbe. 2021 Jan 13;29(1):121-131.e4. doi: 10.1016/j.chom.2020.11.002. Website: http://104.154.134.205:3838/DBsmORF/

- Shen W, Sipos B, Zhao L. SeqKit2: A Swiss army knife for sequence and alignment processing. Imeta. 2024 Apr 5;3(3):e191. doi: 10.1002/imt2.191. 

- Agrawal P, Bhagat D, Mahalwal M, Sharma N, Raghava GPS. AntiCP 2.0: an updated model for predicting anticancer peptides. Brief Bioinform. 2021 May 20;22(3):bbaa153. doi: 10.1093/bib/bbaa153. GitHub: https://github.com/raghavagps/anticp2

- Han B, Zhao N, Zeng C, Mu Z, Gong X. ACPred-BMF: bidirectional LSTM with multiple feature representations for explainable anticancer peptide prediction. Sci Rep. 2022 Dec 19;12(1):21915. doi: 10.1038/s41598-022-24404-1. GitHub: https://github.com/RUC-MIALAB/ACPred-BMF

- Lee B, Shin D. Contrastive learning for enhancing feature extraction in anticancer peptides. Brief Bioinform. 2024 Mar 27;25(3):bbae220. doi: 10.1093/bib/bbae220. GitHub: https://github.com/bzlee-bio/con_ACP

- Veltri D, Kamath U, Shehu A. Deep learning improves antimicrobial peptide recognition. Bioinformatics. 2018 Aug 15;34(16):2740-2747. doi: 10.1093/bioinformatics/bty179. GitHub: https://github.com/dan-veltri/amp-scanner-v2

- Santos-Júnior CD, Pan S, Zhao XM, Coelho LP. Macrel: antimicrobial peptide screening in genomes and metagenomes. PeerJ. 2020 Dec 18;8:e10555. doi: 10.7717/peerj.10555. GitHub: https://github.com/BigDataBiology/macrel 
