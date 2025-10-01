# South African Listeria Outbreak WGS Project

## Introduction
South Africa faced a public health crisis that would become the world's largest recorded outbreak for a particular infection. Unravelling the cause of the outbreak could only be done by WGS of DNA from infected patients. The goal of this script is to:
Confirm the identity of the organism,
Determine the antimicrobial resistance (AMR) profile of these pathogens.
Detect if there might be a toxin that is accelerating the death rate
Suggest antibiotics/treatment options that could be used managing the cases

## Step 1: Creating directories and downloading data.

`mkdir Project && cd Project`

`mkdir raw_data && cd raw_data`

`wget https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/SA_Polony_100_download.sh`

`bash SA_Polony_100_download.sh`

## Step 2: Quality control and combinig qc report
I ran FastQC on all raw reads to assess quality.

`fastqc raw_data/*.fastq.gz -o qc_reports/`

I combined all qc report using multiqc

`multiqc qc_reports`

I opened the multiqc html file showing a combined report. I looked out for good quality of our data and delete those with poor quality from our raw data set. we also delete all data without both the forward and backward sequences.

## Step 3: Trimming
We trimmed adapters and low-quality reads using fastp. first we created a sh document to run through bash

`nano trim.sh`

```
#!/bin/bash
mkdir trimmed_data trimmed_report
for R1 in *_1.fastq.gz
do
    # Detect matching R2 file
    R2=${R1/_1.fastq.gz/_2.fastq.gz}
   sample=$(basename "$R1" _1.fastq.gz)

  fastp \
        -i "$R1" \
        -I "$R2" \
        -o trimmed_data/${sample}_1.trim.fastq.gz \
        -O trimmed_data/${sample}_2.trim.fastq.gz \
        -h trimmed_reports/${sample}_fastp.html \
        -j trimmed_reports/${sample}_fastp.json \
        --thread 8
done

```

save and close the bash script

`bash trim.sh`

## Step 4: Genome Assembly (SPAdes)
I assembled my genome using spades. First I created a bash script then I ran the script

`nano spades.sh`
```
#!/bin/bash
mkdir -p spades_output

for R1 in trimmed_data/*_1.trim.fastq.gz
do
    R2=${R1/_1.trim.fastq.gz/_2.trim.fastq.gz}
    sample=$(basename "$R1" _1.trim.fastq.gz)
    echo "Running SPAdes for $sample..."
    
    spades.py \
      -1 "$R1" \
      -2 "$R2" \
      -o spades_output/$sample \
      --threads 8 \
      --phred-offset 33\
      --memory 32 

done

```

`bash spades.sh`

I combined our contigs.fasta into a folder so I can run futher steps easily

`nano compile.sh`

```
#!/bin/bash
# Gather contigs.fasta from SPAdes outputs and rename with sample name

mkdir -p assemblies
for d in spades_output/*; do
    if [ -f "$d/contigs.fasta" ]; then
        sample=$(basename "$d")
        cp "$d/contigs.fasta" assemblies/${sample}_contigs.fasta
        echo "Copied $d/contigs.fasta to assemblies/${sample}_contigs.fasta"
    fi
done
```

I downloded the contigs.fasta file then copied the nucleotides into the NCBI BLAST site (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&BLAST_SPEC=&LINK_LOC=blasttab&LAST_PAGE=tblastn). This gives the closest organism and percentage of similarity between the genomes. There was strong similarity between our samples and _Listeria monocytogenes_ . The predominant strains were FDA00006667, FDA00009448, BL88/034, BL88/097 and 11-04869. One sample had a stron similarity to the _Listeria _ phage LP-030-2 suggesting a lysogenic relationship considering that this sample also had about 98% similarity to the organism.

## Step 5: Evaluating my Assembly (QUAST)
I evaluated my assembly using quast

`nano run_quast.sh`

```
#!/bin/bash
# Run QUAST on all assemblies in the "assemblies" folder

# Create output directory
mkdir -p quast_results

# Loop through assemblies and run QUAST
for fasta in assemblies/*.fasta; do
    sample=$(basename "$fasta" _contigs.fasta)
    echo "Running QUAST on $sample ..."
    quast.py "$fasta" -o quast_results/${sample}_quast
done

echo "All QUAST analyses completed! Results are in quast_results/"

```
`bash run_quast.sh`

I used a bash script to combine the quast reports
`mkdir quast_combined`
`quast.py assemblies/*.fasta -o quast_combined`

## Step 6: Antimcrobial resistance detection (ABRicate)

`nano abricate.sh`

```
#!/bin/bash

# Create output folder for abricate results
mkdir -p abricate_results

# Loop through all assemblies inside the assembly folder
for assembly in assemblies/*.fasta
do
    # Extract sample name (remove .fasta extension)
    sample=$(basename $assembly .fasta)

    echo "Running ABRicate for sample: $sample"

    # Run abricate and save output
    abricate $assembly > abricate_results/${sample}_abricate.txt
done

#Summarize all results into one table

`abricate --summary abricate_results/*.txt > abricate_results/abricate_summary.tab`
```

`bash abricate.sh`

This gives result for each sample and a summary of everything in a file. From our abricate reports summary we see the antimicrobial resistant gene present and the percentage coverage and identity. Genes with >90% are considered as being present and pronounced. This indicates that both antibiotics can not be used in treatment.We saw Fosfomycin and linomycin with genes fosx and lmo0919_fam respectively. Research has shown that the  first line therapy for _Listeria monocytogenes_ is typically ampicillin (amoxicillin) with or without gentamicin. Since there is no resistance to this gene, we can suggest this as our treatment. For patients allergic to betalactams, Trimethoprim-sulfamethoxazole (TMP-SMX) is an alternative.


## Step 7: Toxin screening.
We screen for virulence factoes with VFD in abricate

`abricate --db vfdb assemblies/*.fasta > abricate_vfdb_results.tab`

compile all results in a .txt format

`abricate --summary abricate_vfdb_results.tab > abricate_vfdb_summary.txt`

There were various toxins and virulent factors present in the samples and at different percentages (>90%). The major classical toxins among others were:

_hly_ (listeriolysin O) → pore-forming toxin, key for escaping the phagosome.

_plcA, plcB_ (phospholipases) → break down host cell membranes, helping escape & spread.

_llsA–llsY_ (Listeriolysin S operon) → produces a bacteriocin-like peptide Listeriolysin S, thought to increase virulence and help survival in the gut.

_iap/cwhA_ (invasion-associated protein, autolysin) → helps with adhesion and invasion, sometimes linked to toxin-like cell wall activity.

_mpl_ (metalloprotease) → activates PlcB toxin.

Other virulence genes found include:
Invasion & spread proteins: internalins _(inlA–inlK), actA, lap_.

Regulators: _prfA_.

Stress/adaptation factors: _clp family, bsh, oatA, pdgA, prsA2_.


# SUMMARY
✅ 60 samples were provided but after quality control only 48 samples were viable for WGS
✅ Work flow: Fastqc ➡️ Fastp ➡️ SPAdes ➡️ QUAST ➡️ ABRicate ➡️ ABRicate (VFDB)
✅ Resiistant genes foxp and lmo0919_fam were found indicating fosfomycin and linomycin resistance.
✅ Public health concerns: High virulence genes and toxin genes explain severe neonatal infections and surveilance and outbreak control methods are critical.

