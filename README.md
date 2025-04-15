# Beagle_phasing
## To start
#### Beagle will be used in conda for this tutorial, so it is necessary to install Anaconda. Java is also required.

### Steps:
1. **Creating VCF file with Plink**;
2. **Creating Beagle environment and installing the package**;
3. **Phasing**.

## 1 - Creating VCF file with Plink
#### First, I need to extract the samples for the analysis and replace the correct IDs.
#### Extracting the samples that are needed based on the breed (in this case, TPR) and create the `sample.ped` file:
```bash
grep '^TPR' firstinput.ped > sample.ped
```
#### Replacing the IDSample with the MATR and create the `updated_sample.ped` file
#### `SampleID_ID.txt`
```bash
IDsample MATR 
horse001 TPR001
horse002 TPR002
horse003 TPR003

```
```bash
awk 'NR==FNR {map[$1]=$2; next} {if ($2 in map) $2=map[$2]; print}' SampleID_ID.txt sample.ped > updated_sample.ped
```
#### Finally, we can create the cleaned `input.vcf` file (filtered for MAF, mind, and geno), which will be the final file used for Beagle:
```plink
plink --ped update_sample.ped --map map.map --geno 0.1 --mind 0.1 --maf 0.05 --chr-set 32 --recode vcf --out input
```
## 2 - Creating beagle enviroment and install package
#### First anaconda need to be installe in the device, then we can menage to install all the necessary package and beagle.
```python
conda create -n beagle-env openjdk bcftools htslib -c bioconda -c conda-forge
```
#### This command: 
- create an enviroment called  `beale-env`;
- install `openjdk` (necessary to launch `.jar`);
- install `bcftools` and htslib necessary for working with VCF.
> [!WARNING]
> Note: Beagle is not avaible as a conda package, it's necessary to dowload the `.jar` manually.
```python
conda activate beagle-env
```
#### Next step is to activate the conda enviroment
#### Then we dowload Beagle using wget command. 
```python
wget https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar
```
> [!NOTE]
> In this case we are using the Beagle.27Feb25.75f version
#### Once is all set we can lauch the script `launchbeagle.sh` to run beagle 
> [!NOTE]
> make sure the file ar compressed and have a index (.tbi)
```python
bgzip input.vcf
tabix -p vcf input.vcf.gz
```
## 3 - Phasing & launchbeagle script 
#### The `launchbeagle.sh` script will check for duplicate markers and launch Beagle to phase the `input.vcf.gz` file. Below are the explanations of the script.
#### Header extraction
```bash
bcftools view -h "$VCF_ORIG" > "$OUTDIR/header.vcf"
```
#### Excluding duplicates
```bash
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\t%FORMAT[\t%GT]\n' "$VCF_ORIG" \
  | awk '!seen[$1"\t"$2]++' > "$OUTDIR/body.vcf"
```
#### Create the final VCF file
```bash
cat "$OUTDIR/header.vcf" <(awk 'BEGIN {OFS="\t"} {print $0}' "$OUTDIR/body.vcf") \
  | bgzip > "$VCF_NODUP"
```
```bash
tabix -p vcf "$VCF_NODUP"
```
#### The cleaned VCF file is saved, and Beagle phasing can now start.
## End of `launchBeagle.sh` and the phased output file is correct.

# Analyzing trio data in VCF files for assessing Mendelian consistency and phasing accuracy
#### To get started, you'll need:
1 - Phased VCF/BCF file: This is typically the output from phasing software like Beagle, where genotypes are marked with | for phased genotypes.
2 - Pedigree file: It defines the relationships within the family (ID, father, mother). You've already shared an example format for this.
### Step 1 : Prepare files
#### The VCF file should have a similar structure:
```bash
##fileformat=VCFv4.2
##filedate=20250401
##source="beagle.28Jun21.220.jar"
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated ALT Allele Frequencies">
##INFO=<ID=DR2,Number=A,Type=Float,Description="Dosage R-Squared: estimated squared correlation between estimated REF dose [P(RA) + 2*P(RR)] and true REF dose">
##INFO=<ID=IMP,Number=0,Type=Flag,Description="Imputed marker">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DS,Number=A,Type=Float,Description="estimated ALT dose [P(RA) + 2*P(AA)]">
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  ID_18EQ10137/02        ID_18EQ10143/01        ID_19EQ30390/01        ID_19EQ31113/01        ID_21EQ01288/02        ID_21EQ01289>
1       14892   UKUL1   G       A       .       PASS    .       GT      0|0     0|0     0|1     0|0     1|1     1|0     0|0     0|1     0|0     0|0     0|0     0|0     0|1     1|0     0|1     0|0     0|0  >
1       29374   UKUL3   G       A       .       PASS    .       GT      0|0     0|1     1|0     0|0     0|0     0|1     1|0     1|0     0|1     1|1     1|1     0|1     1|0     0|1     1|0     1|1     1|0  >
1       115704  BIEC2_25        A       G       .       PASS    .       GT      0|0     1|0     0|0     0|0     0|0     0|0     0|0     0|0     1|0     0|0     0|0     1|0     0|0     0|0     0|0     0|0  >
```
#### The ped file should have a similar structure
```bash
child1  father1  mother1  
father1  0        0       
mother1  0        0        
```
