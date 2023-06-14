#!/usr/bin/env python3
import argparse
import sys
import spectreCNV
import spectreCNVPopulation
import os
import pysam
import logging as logger
from util.metadata.metadataCollector import FastaRef
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
        self.min_cnv_len = 1000000
        self.call_from_console = False
        # dev/debug + hidden params
        self.as_dev = False
        self.max_std_outlier_rm = 5
        self.mosdepth_cov_genome_chr_diff = 0.10  # 10%
        self.lower_2n_threshold = 1.5
        self.upper_2n_threshold = 2.5
        self.cov_diff_threshold = 0.80
        self.dist_proportion = 0.25
        self.candidate_final_threshold = 100000  # 100kb
        self.population_mosdepth = []
        self.threads = 1
        self.run_population_mode = False

    def set_params_from_args(self, user_args):
        self.bin_size = user_args.bin_size
        self.coverage_dir = user_args.coverage_dir
        self.sample_id = user_args.sample_id
        self.out_dir = user_args.output_dir
        self.reference = user_args.reference
        self.metadata = user_args.metadata
        self.snv = user_args.snv_file
        self.black_list = user_args.black_list_file
        self.only_chr_list = user_args.only_chr_list
        self.ploidy = user_args.ploidy
        self.min_cnv_len = user_args.min_cnv_len
        self.n_size = user_args.n_size
        self.threads = user_args.threads
        self.run_population_mode = user_args.run_population_mode

        # dev/debug + hidden params
        self.as_dev = user_args.as_dev
        self.max_std_outlier_rm = user_args.max_std_outlier_rm  # 5
        self.mosdepth_cov_genome_chr_diff = user_args.mosdepth_cov_genome_chr_diff  # 0.10  # 10%
        self.lower_2n_threshold = user_args.lower_2n_threshold  # 1.5
        self.upper_2n_threshold = user_args.upper_2n_threshold  # 2.5
        self.cov_diff_threshold = user_args.cov_diff_threshold  # 0.80
        self.dist_proportion = user_args.dist_proportion  # 0.25
        self.candidate_final_threshold = user_args.candidate_final_threshold  # 100000  # 100kb



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
        self.bin_size = user_args.bin_size
        self.metadata = user_args.metadata
        self.out_dir = user_args.output_dir
        self.n_size = user_args.n_size
        self.save_only = user_args.save_only
        self.black_list = user_args.black_list_file
        self.as_dev = user_args.as_dev
        self.call_from_console = metadata_from_console


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
    worker = spectreCNV.SpectreCNV(coverage_dir=si["coverage_dir"], bin_size=si["bin_size"],
                                   out_dir=si["out_dir"], metadata_file_fasta=si["metadata_file_fasta"],
                                   genome_info=si["genome_info"], sample_id=si["sample_id"],
                                   snv_file_vcf=si["snv_file_vcf"], only_chr_list=si["only_chr_list"],
                                   ploidy=si["ploidy_arg"],min_cnv_len=si["min_cnv_len"], as_dev=si["as_dev"],
                                   dev_params=si["dev_params"], debug_dir=si["debug_dir"])
    worker.cnv_call()
    return worker.cnv_analysis.intermediate_candidates_file_location


class Spectre:

    def __init__(self, as_dev=False):
        # version
        self.version = "0.1-alpha"
        self.logger = logger
        self.debug_dir = ""
        # for spectre cnv caller
        self.spectre_args = SpectreCallParam()
        self.sample_dir_list = []
        # for metadata/removeNs
        self.metadata_args = SpectreMetadataParam()
        # metadata from reference genome (Ns)
        self.__mdr = dict()
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

    @staticmethod
    def make_genome_info(genome_info_pysam):
        genome_info = {"chromosomes": genome_info_pysam.references,
                       "chr_lengths": genome_info_pysam.lengths,
                       "chr_lengths_by_name": {}}
        for chr_name, chr_len in zip(genome_info_pysam.references, genome_info_pysam.lengths):
            genome_info["chr_lengths_by_name"][chr_name] = chr_len
        return genome_info

    def meta_data_extraction(self):
        # ----------- Metadata extraction from ref file  -----------
        # self.metadata_args.as_dev  # Note: this is not used
        reference_fas = self.metadata_args.reference
        bin_size = self.metadata_args.bin_size
        output_dir = os.path.abspath(os.path.expanduser(self.metadata_args.out_dir))
        threshold = self.metadata_args.n_size
        # if "call_from_console" then the meta_data_report serves as output only, otherwise as both
        meta_data_report = self.metadata_args.metadata if not self.metadata_args.call_from_console \
            else f'{output_dir}/{self.metadata_args.metadata}'
        default_metadata_name = f'{output_dir}/metadata.mdr'  # default
        save_only = self.metadata_args.save_only
        fasta_metadata = FastaRef()
        blacklist_data_bed = self.metadata_args.black_list  # bed format
        self.logger.info("Extraction of metadata is activated")

        # metadata parameter given?
        if meta_data_report != "":
            meta_data_report = os.path.abspath(os.path.expanduser(meta_data_report))
        else:
            self.logger.info("Looking for default metadata.mdr")
            meta_data_report = os.path.abspath(os.path.expanduser(default_metadata_name))
        # metadata file exists
        if not os.path.exists(meta_data_report):
            self.logger.info(f'Extracting metadata from {reference_fas}')
            metadata_result = fasta_metadata.get_n_regions(filepath=reference_fas, report_output_dir=output_dir,
                                                           out_file_name=meta_data_report, threshold=threshold,
                                                           bin_size=bin_size, save_only=save_only)
        else:
            self.logger.info(f'Extracting metadata from {meta_data_report}')
            metadata_result = fasta_metadata.extract_n_regions_from_report(meta_data_report)
        if blacklist_data_bed != "":
            self.logger.debug("Using blacklist")
            blacklist_results = fasta_metadata.extract_blacklisted_regions(blacklist_data_bed)
            metadata_result = fasta_metadata.merge_metadata(metadata_result, blacklist_results)
        # return metadata object (dict) after writing to file?
        if save_only:
            pass
        else:
            self.logger.debug("returned meta Object")
            return metadata_result

    def __set_genome_info(self, reference):
        pysam_genome = pysam.FastaFile(reference)
        self.genome = self.make_genome_info(pysam_genome)
        pysam_genome.close()

    def spectre_exe(self):
        # Parameters
        self.display_version()
        self.logger.info("Spectre enabled")
        bin_size = self.spectre_args.bin_size
        coverage_dirs = [os.path.abspath(os.path.expanduser(directory)) for directory in self.spectre_args.coverage_dir]
        sample_ids = self.spectre_args.sample_id
        output_dir = os.path.abspath(os.path.expanduser(self.spectre_args.out_dir))
        reference = os.path.abspath(os.path.expanduser(self.spectre_args.reference))
        metadata = self.spectre_args.metadata
        snv_file = self.spectre_args.snv
        black_list_bed = self.spectre_args.black_list
        only_chr_list = self.spectre_args.only_chr_list
        ploidy_arg = self.spectre_args.ploidy
        min_cnv_len = self.spectre_args.min_cnv_len
        threads = self.spectre_args.threads
        as_dev = self.spectre_args.as_dev
        run_population_mode = self.spectre_args.run_population_mode
        # get the metadata
        self.metadata_args.reference = reference
        self.metadata_args.bin_size = bin_size
        self.metadata_args.metadata = metadata
        self.metadata_args.black_list = black_list_bed
        self.metadata_args.out_dir = output_dir
        self.metadata_args.n_size = self.spectre_args.n_size
        self.metadata_args.save_only = self.spectre_args.save_only
        self.__mdr = self.meta_data_extraction()
        self.__set_genome_info(reference)  # genome information

        # Setting up directories
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        if as_dev:
            self.debug_dir = f"{output_dir}/debug"
            if not os.path.exists(self.debug_dir):
                os.makedirs(self.debug_dir)
        else:
            self.debug_dir = output_dir

        # Preparing spectre instructions for multiprocess
        spectre_instructions = []
        for sample_id, coverage_dir in zip(sample_ids, coverage_dirs):
            instructions = {"coverage_dir": coverage_dir, "bin_size": bin_size,
                            "metadata_file_fasta": self.__mdr.copy(),
                            "out_dir": output_dir, "genome_info": self.genome.copy(),
                            "sample_id": sample_id, "snv_file_vcf": snv_file, "only_chr_list": only_chr_list,
                            "ploidy_arg": ploidy_arg, "as_dev": as_dev, "dev_params": self.spectre_args,
                            "debug_dir": self.debug_dir, "min_cnv_len":min_cnv_len}
            spectre_instructions.append(instructions.copy())

        # Distribute Samples over cores/threads
        with Pool(processes=threads) as pool:
            results = pool.map(outside_spectre_worker, tuple(spectre_instructions))
        intermediate_file_paths = [x for x in results]

        if run_population_mode and len(intermediate_file_paths) > 1:
            self.population_exe("population_file", intermediate_file_paths, output_dir, reference, as_dev)

        self.logger.info("Spectre finished")
        sys.exit(0)

    def population_exe(self, population_sample_name="", population_intermediate_files=None, outputdir="",reference="", as_dev=False):
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
                self.__set_genome_info(reference)

        self.logger.info(f"Population mode: Loaded samples {population_paths}")
        # Directory setup
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        if as_dev:
            debug_dir = f"{output_dir}/debug"
            if not os.path.exists(debug_dir):
                os.makedirs(debug_dir)

        # Population mode setup
        __spectre_population_worker = spectreCNVPopulation.SpectrePopulation(sample_id=sample_ids,
                                                                             output_dir=output_dir,
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
                --bin-size     Bin/Window size (same as Mosdepth)
                --coverage     Coverage directory from Mosdepth output. Expects the following files:
                                   <prefix>.mosdepth.summary.txt
                                   <prefix>.regions.bed.gz
                                   <prefix>.regions.bed.gz.csi
                               Can be one or more directories. Example:
                                    --coverage /path/dir1/ /path/dir2/
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
                --snv          VCF file containing the SNV for the same sample CNV want to be called
                --n-size       Length of consecutive Ns (Default = 5)
                --min_cnv_len  Minimum length of CNV (Default 1mb)
                --population   Runs the population mode on all provided samples
                --threads      Amount of threads (This will boost performance if multiple samples are provided)

                
        removeNs:
            Required
                --reference    Reference genome used for mapping
                --output-dir   Output dir
                --output-file  Output file for results
                --bin-size     Bin/Window size (same as Mosdepth)
            Optional
                --blacklist    Blacklist in bed format for sites that will be ignored (Default = "")
                --n-size       Length of consecutive Ns (Default = 5)
                --save-only    Will only save the metadata file and not show the results in screen (Default = False)
                
        population:
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
    subparser_cnv_caller.add_argument('-b', '--bin-size', type=int, required=True, dest='bin_size', default=500,
                                      help='..., default = 1kb')
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
    subparser_cnv_caller.add_argument('-l', '--blacklist', type=str, required=False, dest='black_list_file',
                                      default="",
                                      help='...')
    subparser_cnv_caller.add_argument('-o', '--only-chr', type=str, required=False, dest='only_chr_list', default="",
                                      help='...')
    subparser_cnv_caller.add_argument('-p', '--ploidy', type=int, required=False, dest='ploidy', default=2,
                                      help='..., default = 2')
    subparser_cnv_caller.add_argument('-n', '--n-size', type=int, required=False, dest='n_size', default=5,
                                      help='..., default = 5')
    subparser_cnv_caller.add_argument('-mcl', '--min-cnv-len', type=int, required=False, dest='min_cnv_len', default=1000000,
                                      help='..., default = 1000000')
    subparser_cnv_caller.add_argument('-t', '--threads', type=int, required=False, dest='threads', default=1,
                                      help='..., default = 1')
    subparser_cnv_caller.add_argument('-i', '--population', action='store_true', required=False,
                                      dest='run_population_mode', default=False, help='...s, default = False')

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
    subparser_cnv_caller.add_argument('-07', '--candidate-final-threshold', type=int, required=False,
                                      dest='candidate_final_threshold', default=100000,
                                      help='..., default = 100,000')  # 100kb

    # ############################################################################################ #
    # Metadata to remove Ns
    metadata_help = "..."
    subparser_metadata = subparsers.add_parser("removeNs", help=metadata_help)
    # Required
    subparser_metadata.add_argument('-r', '--reference', type=str, required=True, dest='reference', default="",
                                    help='..., default = None')
    subparser_metadata.add_argument('-d', '--output-dir', type=str, required=True, dest='output_dir', default="",
                                    help='..., default = None')
    subparser_metadata.add_argument('-f', '--output-file', type=str, required=True, dest='metadata', default="",
                                    help='..., default = None')
    subparser_metadata.add_argument('-b', '--bin-size', type=int, required=True, dest='bin_size', default="",
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
    subparser_cnv_population = subparsers.add_parser("population", help=cnv_population_help)
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
    args = parser.parse_args()
    return args, spectre_help


def main():
    spectre_args, spectre_help = get_arguments()
    command = spectre_args.command
    try:
        run_as_dev = spectre_args.as_dev
    except AttributeError:
        run_as_dev = False
    # logger init
    logger.basicConfig(
        format='spectre::%(process)d:%(levelname)s> %(message)s', level=logger.DEBUG, force=True) if run_as_dev else \
        logger.basicConfig(format='spectre::%(levelname)s> %(message)s', level=logger.INFO, force=True)
    logger.debug(f' Debug output is enabled') if run_as_dev else None
    spectre_run = Spectre(run_as_dev)
    min_bin_size = 500
    if command == "CNVCaller":
        logger.error("Bin size too small") if spectre_args.bin_size < min_bin_size else ""
        # ARGS:  bin_size, coverage_file_bed, sample_id, output_dir, reference="", metadata="", ...
        spectre_run.spectre_args.set_params_from_args(spectre_args)
        spectre_run.spectre_exe()
    elif command == "removeNs":
        logger.error("Bin size too small") if spectre_args.bin_size < min_bin_size else ""
        # ARGS:  reference_fas, output_dir, meta_data_report, bin_size=500, threshold=5
        metadata_call_from_console = True
        spectre_run.metadata_args.set_params_from_args(spectre_args, metadata_call_from_console)
        spectre_run.meta_data_extraction()
    elif command == "population":
        logger.error("at least two candidates are required") if len(spectre_args.population_candidates) < 2 else ""
        spectre_run.population_args.set_params_from_args(spectre_args)
        spectre_run.population_exe()
    elif command == "version":
        spectre_run.display_version()
    else:
        logger.error(spectre_help)


if __name__ == '__main__':
    main()
