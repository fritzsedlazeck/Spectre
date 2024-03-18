
# Requirements
# python3 conda install python=3.8.5
# pysam conda install -c bioconda pysam
# pandas  conda install -c anaconda pandas
# numpy conda install -c anaconda numpy
# scipy conda install -c anaconda scipy
# matplotlib  conda install -c anaconda matplotlib

# reference genome and corresponding blacklist
REFERENCE="GRCh38.fasta" or "GRCh38.fasta.gz"
BLACKLIST_GRCH38="./data/grch38_blacklist.bed"  # Optional but recommended

COVERAGE_PATH="mosdepth_dir_from_sample.regions.bed.gz"  # Must exisit
RESULTSD_DIR="results_dir_path"   # Must exisit
SNV_VCF="sample_SNPs.vcf.gz"      # Opcional

# --bin-size need to be the same used in mosdepth
spectre CNVCaller \
    --coverage ${COVERAGE_PATH} \
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

spectre removeNs \
  --reference  ${REFERENCE} \
  --output-dir ${RESULTSD_DIR} \
  --output-file genome.mdr \
    --save-only


# Docker

# Build docker image
docker-compose -f docker-compose-spectre.yaml build

# Run Spectre with the latest docker image from the build process
# 1kb bin size and 10mb min CNV length
docker run -it --rm --name spectre \
-v /mnt/g/giab/hg002/ont-ul/paper/simulations/survivor/benchmark/1k/test_v37_spectre_sim_benchmark_30X_10MB:/input \
-v /mnt/g/giab/hg002/ont-ul/paper/simulations/survivor/benchmark/spectre_out/test_v37_spectre_sim_benchmark_30X_10MB.1.a:/output \
-v /mnt/d/Projects/spectre_data/docker_test/supplementals:/supplementals spectre:latest CNVCaller \
--coverage /input/ \
--sample-id spectre_full_hg002-GRCh37-ONT-1k \
--output-dir /output/ \
--reference /supplementals/grch37.fa \
--metadata /supplementals/grch37_N.mdr \
--blacklist ./data/grch37_blacklist_spectre.bed \
--min-cnv-len 100000

docker run -it --rm --name "spectre" \
-v "${mosdepth_path}":/input \
-v "${output_path}":/output \
-v "${reference_dir}":/supplementals spectre:0.2.5.alpha.1 CNVCaller \
--coverage /input/$mosdepth_coverage_filename \
--sample-id spectre_full-GRCh38-ONT-1k \
--output-dir /output/ \
--reference /supplementals/grch38.fa.gz \
--metadata ./data/grch38.mdr \
--blacklist ./data/grch38_blacklist_spectre.bed \
--ploidy 2 \
--min-cnv-len 100000 \
--ploidy-chr chrX:4 \
--dev >> "${log_path}"