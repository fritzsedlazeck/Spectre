import os
import sys
import pysam
from spectre.analysis.analysis import CNVAnalysis
from spectre.analysis.snv_analysis import SNVAnalysis
from spectre.util import logger


class SpectreCNV:

    def __init__(self, spectre_args, coverage_filepath, sample_id, mdr_file_path, genome_info, debug_dir=""):
        self.logger = logger.setup_log(__name__, spectre_args.as_dev)
        # from params
        self.genome_info = genome_info
        # self.metadata_reference = metadata_file_fasta no need to be saved
        self.sample_id = sample_id if sample_id != "" else "sample-"  # TODO: random string
        self.spectre_args = spectre_args
        # Other needed params
        self.mosdepth_data = None
        # output bed dir/file
        self.out_bed = os.path.join(os.path.join(self.spectre_args.out_dir, f'{self.sample_id}_cnv.bed.gz'))
        # output VCF dir/file
        self.out_vcf = os.path.join(os.path.join(self.spectre_args.out_dir, f'{self.sample_id}_cnv.vcf.gz'))
        # mosdepth coverage file
        self.mosdepth_coverage_file = self.set_coverage_files(coverage_filepath)
        # population
        self.population = {}
        # Analysis params
        # SNV
        self.snv_analysis = None  # init, used later if appropiate
        # CNV
        self.cnv_analysis = CNVAnalysis(coverage_file=self.mosdepth_coverage_file, outbed=self.out_bed,
                                        outvcf=self.out_vcf, sample_id=self.sample_id,
                                        spectre_args=self.spectre_args, metadata_ref=mdr_file_path,
                                        genome_info=genome_info, debug_dir=debug_dir)

    def get_coverage_dir_file(self, coverage_filepath):
        # for each_dir in os.listdir(coverage_dir):

        if coverage_filepath.endswith(".regions.bed.gz"):
            if not os.path.isfile(f"{coverage_filepath}.csi"):
                try:
                    self.logger.warning(f'Coverage file found but index is missing trying to index {coverage_filepath}')
                    pysam.tabix_index(coverage_filepath, preset="bed", force=True,csi=True)
                except AttributeError:
                    self.logger.error(f"Could not create (.csi) index file for: {coverage_filepath}")
                    self.logger.error("Please index the coverage file manually (e.g. Tabix)")
                    sys.exit(1)


    def set_coverage_files(self, coverage_filepath):
        self.logger.info(
            f'Spectre calculating for: {str(coverage_filepath)}')  # and bin size: {self.spectre_args.bin_size}')
        self.get_coverage_dir_file(coverage_filepath)
        return coverage_filepath

    # MAIN analysis pipeline incluide CNV and SNV analysis (CN neutral and LoH)
    def cnv_call(self):
        # SNV analysis first, if SNV data exisist use it to get the CN neutral
        if self.spectre_args.snv != "":
            self.logger.info("Analysing CN neutral state from SNV data")
            self.logger.info("Cancer mode on") if self.spectre_args.is_cancer else None
            self.snv_analysis = SNVAnalysis(snv_file=self.spectre_args.snv, genome=self.genome_info,
                                            args=self.spectre_args)
            if self.spectre_args.cn_neutral_arg != 0.0 and self.spectre_args.loh_only:
                self.snv_analysis.cn_neutral_coverage_med = self.spectre_args.cn_neutral_arg
                self.snv_analysis.cn_neutral_coverage_mean = self.spectre_args.cn_neutral_arg
                self.snv_analysis.cn_neutral_coverage_stdev = 1
            else:
                self.snv_analysis.snv_copy_number_state()
            # "self.snv_analysis.cn_neutral_coverage_med" has the estimated coverage for 2x (ploidy x) by median
            # "self.snv_analysis.cn_neutral_coverage_mean" has the estimated coverage for 2x (ploidy x) by mean
            # "self.snv_analysis.cn_neutral_coverage_stdev" has the estimated stdev for the coverage for 2x (ploidy x)
            self.cnv_analysis.snv_derived_cn_neutral = {"med": self.snv_analysis.cn_neutral_coverage_med,
                                                        "avg": self.snv_analysis.cn_neutral_coverage_mean,
                                                        "std": self.snv_analysis.cn_neutral_coverage_stdev
                                                        }
        else:
            if self.spectre_args.is_cancer:
                self.logger.warning("Cancer mode, SNV data is highly recommended. CN neutral will be estimated from "
                                    "the provided coverage")

        # check if skiping the CNV calling (loh-only)
        if not self.spectre_args.loh_only:
            # Data normalization
            self.logger.info("Data normalization and outlier removal (right tail)")
            self.cnv_analysis.data_normalization()
            # Coverage analysis
            self.logger.info(f"CNV calling - Coverage for sample: {self.sample_id}")
            # get raw CNV for original sample
            self.cnv_analysis.call_cnv_coverage()
            # CNV metrics
            self.cnv_analysis.get_cnv_metrics()
            # refine cnvs
            self.cnv_analysis.refine_cnv_calls()
            # CNV metrics
            # self.logger.warning("Disabled CNV metrics")
            self.logger.info("Calculate CNV metrics")
            self.cnv_analysis.get_cnv_metrics(refined_cnvs=True)

        # SNV analysis
        if self.spectre_args.snv != "":
            self.logger.info("CNV candidates by SNV")
            snv_file_basename_no_dot = "_".join(os.path.basename(self.spectre_args.snv).split('.'))
            self.snv_analysis.snv_af_bed = f'{self.spectre_args.out_dir}/tmp/{snv_file_basename_no_dot}.bed'
            # Compare candidates to AF/SNV
            self.snv_analysis.call_cnv_af_region(self.cnv_analysis.cnv_calls_list)
            self.cnv_analysis.cnv_calls_list_af_filtered = self.snv_analysis.cnv_calls_list_af
            # LoH
            self.logger.info("Starting LOH")
            self.snv_analysis.loh(self.cnv_analysis.existing_cnv_ids)
            self.logger.info("Saving final LOH calls")
            self.cnv_analysis.snv_loh = self.snv_analysis.loh_pass_only()

        # Check for breakpoints
        if self.spectre_args.snfj != "":
            self.logger.info("Checking breakpoints from the SNFJ file")
            self.cnv_analysis.check_snfj_breakpoints(self.spectre_args.snfj)

        # Make output files
        self.logger.info("Merge CNV and LoH candidates")
        self.cnv_analysis.cnv_loh_candidate_merge()
        self.logger.info("Final candidates are written to .spc file")
        self.cnv_analysis.write_intermediate_candidates()
        self.logger.info("Results are writen to bed file")
        self.cnv_analysis.cnv_result_bed()
        self.logger.info("Results are writen to VCF file")
        self.cnv_analysis.cnv_result_vcf()
        self.logger.info("Result plot in progress")
        self.cnv_analysis.cnv_plot()
        # End
        self.logger.info(f"Output dir: {self.spectre_args.out_dir}")
        self.logger.info(f"Done with sample {self.sample_id}")
        # sys.exit(0)
