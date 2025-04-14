# Beagle_phasing
## To start
#### Beagle will be used in conda for this tutorial so is necessary to install anaconda, java is also necessary 
### Step 
#### 1 - Creating VCF file with plink;
#### 2 - Creating Beagle enviroment and install package;
#### 3 - Phasing.
## 1 -  Creating VCF file with plink

## 2 - Creating beagle enviroment and install package
#### First anaconda need to be installe in the device, then we can menage to install all the necessary package and beagle.
```python
conda create -n beagle-env openjdk bcftools htslib -c bioconda -c conda-forge
```
#### This command: 
- create an enviroment called  `beale-env`;
- install `openjdk` (necessary to launch `.jar`);
- install `bcftools` and htslib necessary for working with VCF.
> Note: beagle is not avaible as a conda package, it's necessary to dowload the `.jar` manually.
```python
conda activate beagle-env
```
#### Next step is to activate the conda enviroment
#### Then we dowload Beagle using wget command. 
```python
wget https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar
```
> In this case we are using the Beagle.27Feb25.75f version
#### Once is all set we can lauch the script `launchbeagle.txt` to run beagle 
#make sure the file ar compressed and have a index (.tbi)
bgzip input.vcf
tabix -p vcf input.vcf.gz

## Phasing & launchbeagle script 
#### Need to chaek for dupliccate markers 
bcftools norm -d snps -Oz -o TPR_nodup.vcf.gz TPR.vcf.gz
tabix -p vcf TPR_nodup.vcf.gz

# 🧬 Genotype Phasing with hsphase in R

This repository provides a simple pipeline to perform haplotype phasing on genotypic data using the [`hsphase`](https://cran.r-project.org/package=hsphase) package in R.

## 📁 Contents

- `run_hsphase.R`: Main R script to perform phasing from PLINK `.raw` format
- `haplotypes_hsphase.txt`: Output file with phased haplotypes (generated after running the script)
- `example_data/`: Optional folder for example genotype input files (if you want to include a test dataset)

---

## 📦 Requirements

- R (version ≥ 4.0 recommended)
- R package: `hsphase`
- [PLINK v1.9](https://www.cog-genomics.org/plink/1.9/) or newer

---

## ⚙️ How to Run

### 1. Prepare your genotype data

Use PLINK to export your binary genotype data in `.raw` format:

```bash
plink --bfile your_data --recode A --out TPR_genoA
```
## 2. Run the phasing script

You can run the phasing script from your terminal or directly from an R environment:

```bash
Rscript run_hsphase.R
```
## 3. Output

After running the script, the following output file will be generated:

- `haplotypes_hsphase.txt`:  
  A tab-separated table containing phased haplotypes.  
  - **Rows**: individuals  
  - **Columns**: SNP markers  
  - **Format**: Each cell contains phased allele data (0/1 format)

This output is suitable for downstream haplotype-based analyses or visualization.
