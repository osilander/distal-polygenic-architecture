#!/bin/bash
#SBATCH --job-name=maf_multiwin
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --account=uoa04097
#SBATCH --output=logs/maf_multiwin_%A_%a.out
#SBATCH --error=logs/maf_multiwin_%A_%a.err

set -euo pipefail

module load BCFtools
module load HTSlib
module load BEDTools

mkdir -p logs per-snp-bed windows tmp

LIST="filled-vcf-list.txt"
VCF=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$LIST")

if [ -z "${VCF}" ]; then
  echo "No VCF on line ${SLURM_ARRAY_TASK_ID} of ${LIST}"
  exit 1
fi

PREFIX=$(basename "$VCF" .autosomes.gnomad_af_filled.vcf.gz)

PER_SNP="per-snp-bed/${PREFIX}.maf.per_snp.bed.gz"
ENRICHED="tmp/${PREFIX}.${SLURM_JOB_ID}.per_snp.enriched.bed"

echo "Input:  $VCF"
echo "Prefix: $PREFIX"

SAMPLE=$(bcftools query -l "$VCF" | head -1)
if [ -z "${SAMPLE}" ]; then
  echo "Could not determine sample name from VCF: $VCF"
  exit 1
fi
echo "Sample: $SAMPLE"

# Build base per SNP BED if missing
# Columns: chrom start end es_signed_maf lp maf
if [ ! -s "$PER_SNP" ]; then
  TMP_PER="tmp/${PREFIX}.${SLURM_JOB_ID}.per_snp.bed"

  bcftools view \
    -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22 \
    -m2 -M2 \
    -i 'INFO/GNOMAD_AF!="."' \
    "$VCF" \
  | bcftools query -s "$SAMPLE" -f '%CHROM\t%POS\t%INFO/GNOMAD_AF\t[%ES]\t[%LP]\n' \
  | awk -F'\t' -v OFS='\t' '
      {
        chr=$1; pos=$2; af=$3; es=$4; lp=$5
        if (af=="." || af=="" || es=="." || es=="" || lp=="." || lp=="") next

        afn = af + 0.0
        esn = es + 0.0
        lpn = lp + 0.0

        print chr, pos-1, pos, esn, lpn, afn
      }' > "$TMP_PER"

  bgzip -c "$TMP_PER" > "$PER_SNP"
  tabix -s 1 -b 2 -e 3 "$PER_SNP"
  rm -f "$TMP_PER"
fi

# Build enriched per SNP BED (plain text) once
# Columns:
# 1 chrom
# 2 start
# 3 end
# 4 es_signed_maf
# 5 lp
# 6 maf
# 7 es_abs_maf
# 8 es_randsign_maf  (random sign applied to abs)
#
# Deterministic randomisation (same for all studies) unless you change srand seed.
zcat "$PER_SNP" \
| awk -F'\t' -v OFS='\t' '
  BEGIN{srand(1)}
  {
    es=$4+0.0
    abs=es; if(abs<0) abs=-abs
    s=(rand()<0.5?-1:1)
    rs=abs*s
    print $1,$2,$3,$4,$5,$6,abs,rs
  }' > "$ENRICHED"

# Window sizes
for W in 25000 50000 100000 200000 500000 1000000 2000000 3000000 5000000 7000000 10000000; do
  WIN="grids/windows.${W}.bed"
  if [ ! -s "$WIN" ]; then
    echo "Missing windows file: $WIN"
    exit 1
  fi

  OUT="windows/${PREFIX}.w${W}.summary.tsv"
  TMP_MAP="tmp/${PREFIX}.${SLURM_JOB_ID}.w${W}.mapped.tsv"

  # One bedtools pass per window size
  # - c4 median, stdev
  # - c7 median (abs)
  # - c8 median (random sign)
  # - c5 median, max (LP)
  # - c6 median (MAF)
  # - c4 count (SNP count)
  bedtools map \
    -a "$WIN" \
    -b "$ENRICHED" \
    -c 4,4,7,8,5,5,6,4 \
    -o median,stdev,median,median,median,max,median,count \
    -null NaN > "$TMP_MAP"

  awk -F'\t' -v OFS='\t' '
    BEGIN{
      print "chrom","start","end","median_es_signed","stdev_es_signed","median_es_abs","median_es_randsign","median_lp","max_lp","median_af","count_snps"
    }
    { print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11 }' "$TMP_MAP" > "$OUT"

  rm -f "$TMP_MAP"
  echo "Wrote $OUT"
done

rm -f "$ENRICHED"

echo "Done $PREFIX"
