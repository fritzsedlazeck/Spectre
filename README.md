
![Spectre](./logo.png)
# Spectre - Long read CNV caller

## Required programs (conda)

Setup a conda environment for Spectre (copy and paste the following commands)
```bash
conda create -n spectre python=3.8.5 pysam==0.21.0 numpy==1.24.3 pandas==2.0.1 matplotlib==3.7.1 scipy==1.10.1 -y
conda activate spectre
```
Alternatively, you can use pip for installing the packages stored in the requirements txt

```bash
conda create -n spectre python=3.8.5 pip -y
conda activate spectre
pip install -r requirements.txt
```
or install everything manually (check for package version in the requirements.txt file)

|Program| Conda                               |
|-------|-------------------------------------|
| python3 | conda install python=3.8.5          |
| pysam | conda install -c bioconda pysam=0.21.0     |
| pandas| conda install -c anaconda pandas==2.0.1    |
| numpy| conda install -c anaconda numpy==1.24.3     |
| scipy| conda install -c anaconda scipy==1.10.1     |
| matplotlib| conda install -c anaconda matplotlib==3.7.1 |


## How to run
Spectre need as input:
- The result of Mosdepth (directory)
- Reference genome (can be bgzip compressed)
- Window size used in Mosdepth (Make sure the binsize between Mosdepth and Spectre are matching. We suggest a binsize of 1000 base pairs.)

Optional
- VCF file containing SNV

>INFO: Make sure to include "/" at the end if you are adding directory paths.


```bash
spectre.py CNVCaller \
  --bin-size 1000 \
  --coverage mosdepth/sampleid/ \
  --sample-id sampleid \
  --output-dir sampleid_output_directory_path/ \
  --reference reference.fasta.gz \
  --snv sampleid.vcf.gz
```
### Run Spectre with multiple samples
Run Spectre with multiple samples:
>INFO: This will start the population mode automatically.

```bash
spectre.py CNVCaller \
  --bin-size 1000 \
  --coverage mosdepth/sampleid-1/ mosdepth/sampleid-1/ \
  --sample-id sampleid-1 sampleid-2 \
  --output-dir sampleid_output_directory_path/ \
  --reference reference.fasta.gz \
  --snv sampleid.vcf.gz
```

### Population mode
Run Spectre in population mode with two or more samples:
>INFO: Spectre produces an intermediate file (.spc) which contains all calculated CNVs from a given samples. They are 
> located in the output folder of given sample.

```bash
spectre.py population \
  --candidates /path/to/sample1.spc /path/to/sample2.spc \
  --sample-id output_name \
  --output-dir sampleid_output_directory_path/
```
