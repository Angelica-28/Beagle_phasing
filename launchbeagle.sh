#!/bin/bash

# === Adjustable option ===
VCF_ORIG="path/to/vcffile"
VCF_NODUP="path/to/input_nodup.vcf.gz"
BEAGLE_JAR="path/to/beagle.27Feb25.75f.jar"
OUTDIR="/beagle_output"
OUTPREFIX="input__phased"

WINDOW=50       # Windows cM
OVERLAP=5       # Overlap cM
NE=100         # Ne
THREADS=4

mkdir -p "$OUTDIR"

echo "🔍 Removing duplicates from $VCF_ORIG..."

# Header extraction
bcftools view -h "$VCF_ORIG" > "$OUTDIR/header.vcf"

# excluding duplicates
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\t%FORMAT[\t%GT]\n' "$VCF_ORIG" \
  | awk '!seen[$1"\t"$2]++' > "$OUTDIR/body.vcf"

# Create the final VCF file
cat "$OUTDIR/header.vcf" <(awk 'BEGIN {OFS="\t"} {print $0}' "$OUTDIR/body.vcf") \
  | bgzip > "$VCF_NODUP"
tabix -p vcf "$VCF_NODUP"

echo "✅ Clean VCF saved in: $VCF_NODUP"

echo "🚀 Starting phasing with Beagle..."

# Phasing
java -Xmx8g -jar "$BEAGLE_JAR" \
  gt="$VCF_NODUP" \
  out="$OUTDIR/$OUTPREFIX" \
  window=$WINDOW \
  overlap=$OVERLAP \
  ne=$NE \
  impute=false \
  nthreads=$THREADS

echo "✅ Phasing completed! Output in: $OUTDIR/${OUTPREFIX}.vcf.gz"



