
# Requirements
# python3 conda install python=3.8.5
# pysam conda install -c bioconda pysam
# pandas  conda install -c anaconda pandas
# numpy conda install -c anaconda numpy
# scipy conda install -c anaconda scipy
# matplotlib  conda install -c anaconda matplotlib

# reference genome and corresponding blacklist
REFERENCE="GRCh38.fasta" or "GRCh38.fasta.gz"
BLACKLIST_GRCH38="./spectre-cnv/data/grch38_blacklist.bed"  # Optional but recommended

COVERAGE_DIR="mosdepth_dir_from_sample"
RESULTSD_DIR="results_dir_path"   # Must exisit
SNV_VCF="sample_SNPs.vcf.gz"      # Opcional

# --bin-size need to be the same used in mosdepth
python3 ./spectre.py CNVCaller \
  --bin-size 1000 \
  --coverage ${COVERAGE_DIR} \
  --output-dir ${RESULTSD_DIR} \
  --sample-id sample_name \
  --reference  ${REFERENCE} \
  --snv ${SNV_VCF} \
  --black_list ${BLACKLIST_GRCH38}

# Other optional parameters are:
# --only_chr  for a list of chromosomes to analyze, comma separated, useful for sex chromosomes
# --ploidy    for dealing with sex chromosomes

# Spectre will create a sample.mdr file with genomic regions that contain "N" and thus not mappable
# this file can be later input with --metadata parameter, otherwise it will be created every time.
# For multiple samples run, it is recommended to run the metadata command first:

python3 ./spectre.py removeNs \
  --reference  ${REFERENCE} \
  --output-dir ${RESULTSD_DIR} \
  --output-file genome.mdr \
  --bin-size 1000 \
  --save-only
