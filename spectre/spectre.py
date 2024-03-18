#!/usr/bin/env python3
import os
import sys
sys.dont_write_bytecode = True
import argparse
import pysam

from builtins import float

from spectre.util import logger
from spectre.spectreCNV import SpectreCNV
from spectre.spectreCNVPopulation import SpectrePopulation
from spectre.util.metadata.metadataCollector import FastaRef
from multiprocessing import Pool


class SpectreCallParam(object):
    def __init__(self):
        self.bin_size = 1  # in kb
        self.coverage_dir = ""
        self.sample_id = ""
        self.out_dir = ""
        self.reference = ""
        self.metadata = ""
        self.snv = ""
        self.n_size = 5
        self.save_only = False  # this one is not updated as we need the return to occur
        self.black_list = ""
        self.only_chr_list = ""
        self.ploidy = 2
        self.ploidy_chr_list = ""
        self.min_cnv_len = 1000000
        self.sample_coverage_overwrite = None
        self.call_from_console = False
        self.is_cancer = False
        # dev/debug + hidden params
        self.as_dev = False
        self.max_std_outlier_rm = 5
        self.mosdepth_cov_genome_chr_diff = 0.10  # 10%
        self.lower_2n_threshold = 1.5
        self.upper_2n_threshold = 2.5
        self.cov_diff_threshold = 0.80
        self.dist_proportion = 0.25
        self.candidate_final_threshold = 10000  # 10kb
        self.loh_min_snv_perkb = 10
        self.loh_min_snv_total = 100
        self.loh_min_region_size = 100000
        self.cn_neutral_arg = 0.0
        self.population_mosdepth = []
        self.threads = 1
        self.run_population_mode = False
        self.disable_max_coverage = False

    def set_params_from_args(self, user_args):
        # self.bin_size = user_args.bin_size
        self.coverage_dir = user_args.coverage_dir
        self.sample_id = user_args.sample_id
        self.out_dir = user_args.output_dir
        self.reference = user_args.reference
        self.metadata = user_args.metadata
        self.snv = user_args.snv_file
        self.snfj = user_args.snfj_file
        self.black_list = user_args.black_list_file
        self.only_chr_list = user_args.only_chr_list
        self.ploidy = user_args.ploidy
        self.ploidy_chr_list = user_args.ploidy_chr
        self.min_cnv_len = user_args.min_cnv_len
        self.sample_coverage_overwrite = user_args.sample_coverage_overwrite
        self.n_size = user_args.n_size
        self.threads = user_args.threads
        self.run_population_mode = user_args.run_population_mode
        self.disable_max_coverage = user_args.disable_max_coverage
        self.loh_min_snv_perkb = user_args.loh_min_snv_perkb
        self.loh_min_snv_total = user_args.loh_min_snv_total
        self.loh_min_region_size = user_args.loh_min_region_size
        self.is_cancer = user_args.is_cancer
        self.loh_only = user_args.loh_only
        self.cn_neutral_arg = user_args.cn_neutral
        # dev/debug + hidden params
        self.as_dev = user_args.as_dev
        self.dev_max_std_outlier_rm = user_args.max_std_outlier_rm  # 5
        self.dev_mosdepth_cov_genome_chr_diff = user_args.mosdepth_cov_genome_chr_diff  # 0.10  # 10%
        self.dev_lower_2n_threshold = user_args.lower_2n_threshold  # 1.5
        self.dev_upper_2n_threshold = user_args.upper_2n_threshold  # 2.5
        self.dev_cov_diff_threshold = user_args.cov_diff_threshold  # 0.80
        self.dev_dist_proportion = user_args.dist_proportion  # 0.25
        # self.dev_candidate_final_threshold = user_args.candidate_final_threshold  # 100000  # 100kb


class SpectreMetadataParam(object):
    def __init__(self):
        self.reference = ""
        self.bin_size = 1
        self.metadata = ""
        self.out_dir = ""
        self.n_size = 5
        self.save_only = False
        self.as_dev = False
        self.black_list = ""
        self.call_from_console = False

    def set_params_from_args(self, user_args, metadata_from_console=False):
        self.reference = user_args.reference
        # self.bin_size = user_args.bin_size
        self.metadata = user_args.metadata
        self.out_dir = user_args.output_dir
        self.n_size = user_args.n_size
        self.save_only = user_args.save_only
        self.black_list = user_args.black_list_file
        self.as_dev = user_args.as_dev
        self.call_from_console = metadata_from_console

    def set_params_from_spectre(self, user_spectre_args):
        self.reference = user_spectre_args.reference
        # self.bin_size = user_spectre_args.bin_size
        self.metadata = user_spectre_args.metadata
        self.out_dir = user_spectre_args.out_dir
        self.n_size = user_spectre_args.n_size
        self.save_only = user_spectre_args.save_only
        self.black_list = user_spectre_args.black_list
        self.as_dev = user_spectre_args.as_dev


class SpectrePopulationMode(object):
    def __init__(self):
        self.candidates = []
        self.sample_id = ""
        self.out_dir = ""
        self.reference = ""
        self.as_dev = False

    def set_params_from_args(self, user_args):
        self.candidates = user_args.population_candidates
        self.sample_id = user_args.sample_id
        self.out_dir = user_args.output_dir
        self.reference = user_args.reference
        self.as_dev = user_args.as_dev


def outside_spectre_worker(si: dict):
    worker = SpectreCNV(spectre_args=si["spectre_args"], mdr_file_path=si["mdr_file_path"],
                        coverage_filepath=si["coverage_filepath"], sample_id=si["sample_id"],
                        genome_info=si["genome_info"], debug_dir=si["debug_dir"])
    worker.cnv_call()
    return worker.cnv_analysis.intermediate_candidates_file_location


class Spectre:
    def __init__(self, as_dev=False):
        # version
        self.version = "0.2.0"
        # get a custom logger & set the logging level
        self.logger = logger.setup_log(__name__, as_dev)

        self.debug_dir = ""
        # for spectre cnv caller
        self.spectre_args = SpectreCallParam()
        self.sample_dir_list = []
        # for metadata/removeNs
        self.metadata_args = SpectreMetadataParam()
        # metadata from reference genome (Ns)
        self.mdr_file_path = dict()
        # metadata from reference genome for VCF
        self.genome = dict()
        # for benchmark
        self.benchmark = {
            "1": 1
        }
        # for population mode
        self.population_args = SpectrePopulationMode()

    def display_version(self):
        self.logger.info(f'Spectre version: {self.version}')

    def make_genome_info(self):
        MIN_CHR_LEN = 1e6
        pysam_genome = pysam.FastaFile(self.spectre_args.reference)
        genome_info = {"chromosomes": [], "chr_lengths": {}}
        for chr_name, chr_len in zip(pysam_genome.references, pysam_genome.lengths):
            if chr_len > MIN_CHR_LEN:
                genome_info["chr_lengths"][chr_name] = chr_len
                genome_info["chromosomes"].append(chr_name)
        pysam_genome.close()
        self.logger.debug(f'genome: {genome_info["chromosomes"]}')
        self.genome = genome_info

    # TODO: we need to remove this and give the mdr metadata + blacklist in a single file
    def extact_metadata_from_reference(self) -> str:
        """
        First step of Spectre, extract metadata from reference genome if no md file was provided
        :return: path to the metadata file (optional)
        """
        # ----------- Metadata extraction from ref file  -----------
        self.metadata_args.out_dir = os.path.abspath(os.path.expanduser(self.metadata_args.out_dir))

        # Extract regions with Ns
        fasta_metadata = FastaRef(metadata_args=self.metadata_args)
        self.logger.info("Evaluate if a new .mdr file needs to be created")

        # Determine if a mdr file was given if not create one in the output directory
        mdr_file_path = self.metadata_args.metadata if not self.metadata_args.call_from_console \
            else f'{self.metadata_args.out_dir}/{self.metadata_args.metadata}'
        # metadata parameter given?
        if mdr_file_path != "":
            mdr_file_path = os.path.abspath(os.path.expanduser(mdr_file_path))
        else:
            self.logger.info("Looking for default metadata.mdr in output directory")
            mdr_file_path = os.path.abspath(os.path.expanduser(f'{self.metadata_args.out_dir}/metadata.mdr'))

        # Check if mdr file exists
        if not os.path.exists(mdr_file_path):
            self.logger.info(f'No metadata file found in')
            self.logger.info(f'Extracting metadata from {self.metadata_args.reference}')
            fasta_metadata.get_n_regions(out_file_name=mdr_file_path)
        else:
            self.logger.info(f'Using existing metadata file {self.metadata_args.metadata}')
        return mdr_file_path

    # TODO: same here we need to remove this and give the mdr metadata + blacklist in a single file
    def meta_data_extraction(self):
        # ----------- Metadata extraction from ref file  -----------
        # self.metadata_args.as_dev  # Note: this is not used
        self.metadata_args.out_dir = os.path.abspath(os.path.expanduser(self.metadata_args.out_dir))
        blacklist_data_bed = self.metadata_args.black_list  # bed format
        # Extract regions with Ns
        fasta_metadata = FastaRef(metadata_args=self.metadata_args)
        self.logger.info("Extraction of metadata is activated")
        # if "call_from_console" then the meta_data_report serves as output only,
        # otherwise as both
        default_metadata_name = f'{self.metadata_args.out_dir}/metadata.mdr'  # default
        meta_data_report = self.metadata_args.metadata if not self.metadata_args.call_from_console \
            else f'{self.metadata_args.out_dir}/{self.metadata_args.metadata}'
        # metadata parameter given?
        if meta_data_report != "":
            meta_data_report = os.path.abspath(os.path.expanduser(meta_data_report))
        else:
            self.logger.info("Looking for default metadata.mdr")
            meta_data_report = os.path.abspath(os.path.expanduser(default_metadata_name))
        # TODO metadata (.mdr) file does not exist
        # TODO get true positions of Ns
        if not os.path.exists(meta_data_report):
            self.logger.info(f'Extracting metadata from {self.metadata_args.reference}')
            metadata_result = fasta_metadata.get_n_regions(out_file_name=meta_data_report)
        else:
            # TODO get true positions based on metadata (.mdr) file, shall only be called when computing the sample
            self.logger.info(f'Extracting metadata from {meta_data_report}')
            metadata_result = fasta_metadata.extract_n_regions_from_report(meta_data_report)

        if self.metadata_args.black_list != "":
            self.logger.debug("Using blacklist")
            blacklist_results = fasta_metadata.extract_blacklisted_regions()
            metadata_result = fasta_metadata.merge_metadata(metadata_result, blacklist_results)
        # return metadata object (dict) after writing to file?
        if self.metadata_args.save_only:
            pass
        else:
            self.logger.debug("returned meta Object")
            return metadata_result

    def spectre_exe(self):
        # Parameters
        self.display_version()
        self.logger.info("Starting spectre")

        self.spectre_args.out_dir = os.path.abspath(os.path.expanduser(self.spectre_args.out_dir))
        self.spectre_args.reference = os.path.abspath(os.path.expanduser(self.spectre_args.reference))

        # TODO check if mdr file does exist if not create it otherwise delegate mdr path to all samples
        # Update metadata parameters as we have init them
        # It is done here in case we need population mode not re-write the file n times
        # same for genome
        self.metadata_args.set_params_from_spectre(self.spectre_args)
        self.mdr_file_path = self.extact_metadata_from_reference()
        self.make_genome_info()  # genome information
        # set the only chromosomes list to be used
        self.spectre_args.only_chr_list = str(self.spectre_args.only_chr_list).split(",") \
            if self.spectre_args.only_chr_list != "" else self.genome["chromosomes"]

        # set the ploidy for each chromosome
        self.genome["chr_ploidy"] = dict([(x, self.spectre_args.ploidy) for x in self.genome["chromosomes"]])
        if self.spectre_args.ploidy_chr_list != "":
            for ploid_chr in self.spectre_args.ploidy_chr_list.split(","):
                # check if key value-pair is correct
                if ":" not in ploid_chr:
                    self.logger.error(f"Wrong format for ploidy_chr parameter.")
                    self.logger.error(f"You provided {self.spectre_args.ploidy_chr_list}")
                    self.logger.error(f"Please use something like this: --ploidy-chr chr1:2,chr2:2,chrX:1,...")
                    sys.exit(1)
                chr_name, ploidy = ploid_chr.split(":")
                # check if chromosome is in reference genome if missing it will be added
                if chr_name not in self.genome["chr_ploidy"].keys():
                    self.logger.warning(f"Chromosome {chr_name} not found in reference genome. Skipping {ploid_chr}!")
                    continue
                self.genome["chr_ploidy"][chr_name] = int(ploidy)

        # Setting up directories
        if not os.path.exists(self.spectre_args.out_dir):
            os.makedirs(self.spectre_args.out_dir)
        # tmp dir
        if not os.path.exists(f'{self.spectre_args.out_dir}/tmp'):
            os.makedirs(f'{self.spectre_args.out_dir}/tmp')

        if self.spectre_args.as_dev:
            self.debug_dir = f"{self.spectre_args.out_dir}/debug"
            if not os.path.exists(self.debug_dir):
                os.makedirs(self.debug_dir)
        else:
            self.debug_dir = self.spectre_args.out_dir

        # Split for pop and single
        spectre_instructions = []
        for sample_id, coverage_dir in zip(self.spectre_args.sample_id, self.spectre_args.coverage_dir):
            coverage_dir = os.path.abspath(os.path.expanduser(coverage_dir))
            instructions = {"spectre_args": self.spectre_args, "coverage_filepath": self.spectre_args.coverage_dir[0], "sample_id": sample_id,
                            "mdr_file_path": self.mdr_file_path, "genome_info": self.genome.copy(),
                            "debug_dir": self.debug_dir, "logger": self.logger
                            }
            spectre_instructions.append(instructions.copy())

        if self.spectre_args.run_population_mode:
            # Preparing spectre instructions for multiprocess
            spectre_instructions = []
            for sample_id, coverage_dir in zip(self.spectre_args.sample_id, self.spectre_args.coverage_dir):
                coverage_dir = os.path.abspath(os.path.expanduser(coverage_dir))
                instructions = {"spectre_args": self.spectre_args, "coverage_dir": coverage_dir, "sample_id": sample_id,
                                "mdr_file_path": self.mdr_file_path, "genome_info": self.genome.copy(),
                                "debug_dir": self.debug_dir, "logger": self.logger
                                }
                spectre_instructions.append(instructions.copy())

            # Distribute Samples over cores/threads
            with Pool(processes=self.spectre_args.threads) as pool:
                results = pool.map(outside_spectre_worker, tuple(spectre_instructions))
            intermediate_file_paths = [x for x in results]

            if self.spectre_args.run_population_mode and len(intermediate_file_paths) > 1:
                self.population_exe("population_file", intermediate_file_paths, self.spectre_args.out_dir, self.spectre_args.reference, self.spectre_args.as_dev)
        else:

            spectre_main = SpectreCNV(spectre_args=self.spectre_args,
                                      mdr_file_path=self.mdr_file_path,
                                      coverage_filepath=self.spectre_args.coverage_dir[0],
                                      sample_id=self.spectre_args.sample_id[0],
                                      genome_info=self.genome.copy(), debug_dir=self.debug_dir)
            spectre_main.cnv_call()

        self.logger.info("Spectre finished")
        # sys.exit(0)

    def population_exe(self, population_sample_name="", population_intermediate_files=None, outputdir="", reference="",
                       as_dev=False):
        self.logger.info("Starting Spectre population mode")
        # Adjusting parameters
        sample_ids = self.population_args.sample_id if population_sample_name == "" else population_sample_name
        population_paths = self.population_args.candidates if not population_intermediate_files else population_intermediate_files
        output_dir = os.path.abspath(os.path.expanduser(self.population_args.out_dir)) if outputdir == "" else outputdir
        as_dev = self.population_args.as_dev if not as_dev else as_dev
        reference = self.population_args.reference if not reference else reference

        # Required to load the genome information, if any other file format than .spc is provided
        if not any(".spc" in s for s in population_paths):
            if not self.genome:
                self.make_genome_info(reference)

        self.logger.info(f"Population mode: Loaded samples {population_paths}")
        # Directory setup
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        if as_dev:
            debug_dir = f"{output_dir}/debug"
            if not os.path.exists(debug_dir):
                os.makedirs(debug_dir)

        # Population mode setup
        __spectre_population_worker = SpectrePopulation(sample_id=sample_ids,
                                                                             output_dir=output_dir,
                                                                             reciprocal_overlap=0.2,
                                                                             genome_info=self.genome)
        __spectre_population_worker.load_files(population_paths)  # Loading files for population mode
        __spectre_population_worker.cnv_call_population()  # Starting the population CNV calculations


# Arguments
def get_arguments():
    spectre_help = """
    vcf_utils <command> [<args>]
    Spectre:
        CNVCaller:
            Required
                --coverage     Path to the coverage file from Mosdepth output. Expects the following files:
                                   <prefix>.regions.bed.gz
                                   <prefix>.regions.bed.gz.csi
                               Can be one or more directories. Example:
                                    --coverage /path/md1.regions.gz /path/md2.regions.gz
                --sample-id    Sample name/ID. Can be one or more ID. Example:
                                    --sample-id id1 id2
                --output-dir   Output directory
                --reference    Reference sequence used for mapping (for N removal)
            Optional, if missing it will be created
                --metadata     Metadata file for Ns removal
            Optional
                --blacklist    Blacklist in bed format for sites that will be ignored (Default = "")
                --only-chr     Comma separated list of chromosomes to use
                --ploidy       Set the ploidy for the analysis, useful for sex chromosomes (Default = 2)
                --ploidy-chr   Comma separated list of key:value-pairs for individual chromosome ploidy control
                               (e.g. chrX:2,chrY:1) If chromosome is not specified, the default ploidy will be used.
                --snv          VCF file containing the SNV for the same sample CNV want to be called
                --snfj         Breakpoints from from Sniffle which has been converted from the SNF to the SNFJ format.
                --n-size       Length of consecutive Ns (Default = 5)
                --min-cnv-len  Minimum length of CNV (Default 100kb)
                --cancer       Set this flag if the sample is cancer (Default = False)
                --population   Runs the population mode on all provided samples
                --threads      Amount of threads (This will boost performance if multiple samples are provided)

                Coverage
                --sample-coverage-overwrite     Overwrites the calculated sample coverage, which is used to normalize
                                                the coverage. e.g. a value of 30 equals to 30X coverage.
                --disable-max-coverage          Disables the maximum coverage check. This will allow to call CNVs

                LoH (requires --snv)
                --loh-min-snv-perkb             Minimum number of SNVs per kilobase for an LoH region (default=5)
                --loh-min-snv-total             Minimum number of SNVs total for an LoH region (default=100)
                --loh-min-region-size           Minimum size of a region for a LoH region (default=100000)


        RemoveNs:
            Required
                --reference    Reference genome used for mapping
                --output-dir   Output dir
                --output-file  Output file for results
            Optional
                --blacklist    Blacklist in bed format for sites that will be ignored (Default = "")
                --n-size       Length of consecutive Ns (Default = 5)
                --save-only    Will only save the metadata file and not show the results in screen (Default = False)

        Population:
            Required
                --candidates   At least 2 candidate files (.spc or .vcf) which should be taken into consideration for the population mode.
                --sample-id    Name of the output file
                --output-dir   Output directory
            Optional
                --reference    Reference sequence (Required if VCF files are used!)
        Version:
            version    Shows current version/build
    """
    parser = argparse.ArgumentParser(
        description="Spectre CNV caller",
        usage=spectre_help
    )
    subparsers = parser.add_subparsers(help=spectre_help, dest="command")

    # ############################################################################################ #
    # Version
    version_help = "Gives the version number"
    subparser_version = subparsers.add_parser("version", help=version_help)
    subparser_version.add_argument('-0', '--0', action='store_true', required=False, dest='_0', default=False, help='')

    # ############################################################################################ #
    # CNV caller
    cnv_caller_help = "..."
    subparser_cnv_caller = subparsers.add_parser("CNVCaller", help=cnv_caller_help)
    # Required
    # subparser_cnv_caller.add_argument('-b', '--bin-size', type=int, required=True, dest='bin_size', default=500,
    #                                   help='..., default = 1kb')
    subparser_cnv_caller.add_argument('-c', '--coverage', type=str, required=True, dest='coverage_dir', default="",
                                      help='..., default = None', nargs='+')
    subparser_cnv_caller.add_argument('-s', '--sample-id', type=str, required=True, dest='sample_id', default="",
                                      help='..., default = None', nargs='+')
    subparser_cnv_caller.add_argument('-d', '--output-dir', type=str, required=True, dest='output_dir', default=".",
                                      help='..., default = None')
    subparser_cnv_caller.add_argument('-r', '--reference', type=str, required=True, dest='reference', default="",
                                      help='..., default = None')
    # Optional, if missing will be created
    subparser_cnv_caller.add_argument('-m', '--metadata', type=str, required=False, dest='metadata', default="",
                                      help='..., default = None')
    # Optional
    subparser_cnv_caller.add_argument('-v', '--snv', type=str, required=False, dest='snv_file', default="",
                                      help='...')
    subparser_cnv_caller.add_argument('-sj', '--snfj', type=str, required=False, dest='snfj_file', default="",
                                      help='...')
    subparser_cnv_caller.add_argument('-l', '--blacklist', type=str, required=False, dest='black_list_file',
                                      default="",
                                      help='...')
    subparser_cnv_caller.add_argument('-o', '--only-chr', type=str, required=False, dest='only_chr_list', default="",
                                      help='...')
    subparser_cnv_caller.add_argument('-p', '--ploidy', type=int, required=False, dest='ploidy', default=2,
                                      help='..., default = 2')
    subparser_cnv_caller.add_argument('-pc', '--ploidy-chr', type=str, required=False, dest='ploidy_chr', default="",
                                      help='..., default = ""')
    subparser_cnv_caller.add_argument('-n', '--n-size', type=int, required=False, dest='n_size', default=5,
                                      help='..., default = 5')
    subparser_cnv_caller.add_argument('-mcl', '--min-cnv-len', type=int, required=False, dest='min_cnv_len',
                                      default=1000000, help='..., default = 1000000')
    subparser_cnv_caller.add_argument('-t', '--threads', type=int, required=False, dest='threads', default=1,
                                      help='..., default = 1')
    subparser_cnv_caller.add_argument('-i', '--population', action='store_true', required=False,
                                      dest='run_population_mode', default=False, help='...s, default = False')

    subparser_cnv_caller.add_argument('-sco', '--sample-coverage-overwrite', type=float, required=False,
                                      dest='sample_coverage_overwrite',
                                      default=None, help='..., numerical values of the average coverage e.g 30=30X')
    subparser_cnv_caller.add_argument('-dmc', '--disable-max-coverage', action='store_true', required=False,
                                      dest='disable_max_coverage', default=False, help='...s, default = False')

    subparser_cnv_caller.add_argument('-a', '--cancer', action='store_true', required=False, dest='is_cancer',
                                      default=False, help='..., default = False')

    # Loh
    subparser_cnv_caller.add_argument('-lohkb', '--loh-min-snv-perkb', type=int, required=False, default=10,
                                      dest='loh_min_snv_perkb', help='..., default = 5')
    subparser_cnv_caller.add_argument('-lohsnv', '--loh-min-snv-total', type=int, required=False, default=100,
                                      dest='loh_min_snv_total', help='..., default = 50')
    subparser_cnv_caller.add_argument('-lohsize', '--loh-min-region-size', type=float, required=False, default=100000.0,
                                      dest='loh_min_region_size', help='..., default = 50000')
    subparser_cnv_caller.add_argument('-lohonly', '--loh-only', action='store_true', required=False, default=False,
                                      dest='loh_only', help='')
    subparser_cnv_caller.add_argument('-cn', '--cn-neutral', type=float, required=False, default=0.0,
                                      dest='cn_neutral', help='')

    # Dev
    subparser_cnv_caller.add_argument('-0', '--dev', action='store_true', required=False, dest='as_dev', default=False,
                                      help='dev, default = False')
    subparser_cnv_caller.add_argument('-01', '--dev-max-std-outlier-rm', type=int, required=False,
                                      dest='max_std_outlier_rm', default=5, help='..., default = 5')
    subparser_cnv_caller.add_argument('-02', '--mosdepth-cov-genome-chr-diff', type=float, required=False,
                                      dest='mosdepth_cov_genome_chr_diff', default=0.10, help='..., default = 0.10')
    subparser_cnv_caller.add_argument('-03', '--lower-2n-threshold', type=float, required=False,
                                      dest='lower_2n_threshold', default=1.5, help='..., default = 2.5')
    subparser_cnv_caller.add_argument('-04', '--upper-2n-threshold', type=float, required=False,
                                      dest='upper_2n_threshold', default=2.5, help='..., default = 2.5')
    subparser_cnv_caller.add_argument('-05', '--cov-diff-threshold', type=float, required=False,
                                      dest='cov_diff_threshold', default=0.80, help='..., default = 0.80')
    subparser_cnv_caller.add_argument('-06', '--dist-proportion', type=float, required=False, dest='dist_proportion',
                                      default=0.25, help='..., default = 0.25')
    # subparser_cnv_caller.add_argument('-07', '--candidate-final-threshold', type=int, required=False,
    #                                   dest='candidate_final_threshold', default=100000,
    #                                   help='..., default = 100000')  # 100kb

    # ############################################################################################ #
    # Metadata to remove Ns
    metadata_help = "..."
    subparser_metadata = subparsers.add_parser("RemoveNs", help=metadata_help)
    # Required
    subparser_metadata.add_argument('-r', '--reference', type=str, required=True, dest='reference', default="",
                                    help='..., default = None')
    subparser_metadata.add_argument('-d', '--output-dir', type=str, required=True, dest='output_dir', default="",
                                    help='..., default = None')
    subparser_metadata.add_argument('-f', '--output-file', type=str, required=True, dest='metadata', default="",
                                    help='..., default = None')
    # Optional
    subparser_metadata.add_argument('-n', '--n-size', type=int, required=False, dest='n_size', default=5,
                                    help='..., default = 5')
    subparser_metadata.add_argument('-s', '--save-only', action='store_true', required=False, dest='save_only',
                                    default=False, help='save_only, default = False')
    subparser_metadata.add_argument('-l', '--blacklist', type=str, required=False, dest='black_list_file', default="",
                                    help='...')
    # Dev
    subparser_metadata.add_argument('-0', '--dev', action='store_true', required=False, dest='as_dev', default=False,
                                    help='dev, default = False')

    # ############################################################################################ #
    # Population mode
    cnv_population_help = "..."
    subparser_cnv_population = subparsers.add_parser("Population", help=cnv_population_help)
    # Required
    subparser_cnv_population.add_argument('-c', '--candidates', type=str, required=True, dest='population_candidates',
                                          default="", help='..., default = None', nargs='+')
    subparser_cnv_population.add_argument('-s', '--sample-id', type=str, required=True, dest='sample_id', default="",
                                          help='..., default = None')
    subparser_cnv_population.add_argument('-d', '--output-dir', type=str, required=True, dest='output_dir', default=".",
                                          help='..., default = None')
    subparser_cnv_population.add_argument('-0', '--dev', action='store_true', required=False, dest='as_dev',
                                          default=False,
                                          help='dev, default = False')
    # Optional
    subparser_cnv_population.add_argument('-r', '--reference', type=str, required=False, dest='reference', default="",
                                          help='..., default = None')

    # ############################################################################################ #
    # Dev
    dev_help = "..."
    subparser_dev = subparsers.add_parser("dev", help=dev_help)

    # ############################################################################################ #
    args = parser.parse_args()
    return args, spectre_help


def main():
    spectre_args, spectre_help = get_arguments()
    command = spectre_args.command
    try:
        run_as_dev = spectre_args.as_dev
    except AttributeError:
        run_as_dev = False
    # no command provided
    if len(sys.argv) < 2:
        print(spectre_help)
        sys.exit(1)
    # logger init
    main_logger = logger.setup_log(__name__, run_as_dev)
    # print Spectre input arguments
    main_logger.info(f"Spectre input: {command} {' '.join(sys.argv[2:])}")
    main_logger.debug(f'Debug output is enabled') if run_as_dev else None
    spectre_run = Spectre(run_as_dev)
    min_bin_size = 100
    if command == "CNVCaller":
        # logging.error("Bin size too small") if spectre_args.bin_size < min_bin_size else ""
        # ARGS:  bin_size, coverage_file_bed, sample_id, output_dir, reference="", metadata="", ...
        spectre_run.spectre_args.set_params_from_args(spectre_args)
        spectre_run.spectre_exe()
    elif command == "RemoveNs":
        main_logger.error("Bin size too small") if spectre_args.bin_size < min_bin_size else ""
        # ARGS:  reference_fas, output_dir, meta_data_report, bin_size=500, threshold=5
        metadata_call_from_console = True
        spectre_run.metadata_args.set_params_from_args(spectre_args, metadata_call_from_console)
        spectre_run.meta_data_extraction()
    elif command == "Population":
        main_logger.error("at least two candidates are required") if len(spectre_args.population_candidates) < 2 else ""
        spectre_run.population_args.set_params_from_args(spectre_args)
        spectre_run.population_exe()
    elif command == "dev":
        spectre_run.loh()
        # spectre_run.display_version()
    elif command == "version":
        spectre_run.display_version()
    else:
        # logging.error(f'No arguments provided.\n{spectre_help}')
        pass
    sys.exit(0)


if __name__ == '__main__':
    main()

