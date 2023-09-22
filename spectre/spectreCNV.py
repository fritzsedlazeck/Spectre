import os
import logging as logger
import spectre.util.mosdepthReader
from spectre.analysis.analysis import CNVAnalysis


class SpectreCNV:

    def __init__(self, coverage_dir, bin_size, out_dir, metadata_file_fasta, genome_info, sample_id="",
                 snv_file_vcf="", only_chr_list="", ploidy=2, min_cnv_len=1000000, as_dev=False, dev_params=None,
                 debug_dir=""):
        self.as_dev = as_dev
        # logger
        logger.basicConfig(level=logger.DEBUG) if as_dev else logger.basicConfig(level=logger.INFO)
        self.logger = logger
        self.genome_info = genome_info
        self.metadata_reference = metadata_file_fasta
        self.sample_id = sample_id if sample_id != "" else "sample-"  # TODO: random string
        self.snv_file_vcf = snv_file_vcf
        self.only_chr_list = only_chr_list
        self.min_cnv_len = min_cnv_len
        self.ploidy = ploidy
        self.snv_file_bed_af = ""
        self.mosdepth_data = None
        self.out_bed = "out.bed"
        self.out_vcf = "out.vcf"
        self.bin_size = bin_size
        self.out_dir = out_dir
        self.__get_config(coverage_dir)  # fill 'self.operation_dict' working variable
        # population
        self.population = {}
        # self.__get_population_config(population)
        # dev
        self.dev_params = dev_params  # all params from spectre_args obj (SpectreCallParam)
        self.debug_dir = debug_dir
        # self.cnv_analysis = None  # TODO init by using CNVAnalysis

        self.cnv_analysis = CNVAnalysis(coverage_file=self.mosdepth_data.coverage_file,
                                        coverage_mosdepth_data=self.mosdepth_data.mosdepth_summary_data,
                                        output_directory=self.out_dir, outbed=self.out_bed, outvcf=self.out_vcf,
                                        bin_size=self.bin_size, genome_info=self.genome_info, sample_id=self.sample_id,
                                        metadata_ref=self.metadata_reference, snv_file=self.snv_file_vcf,
                                        only_chr_list=self.only_chr_list, ploidy=self.ploidy, min_cnv_len=min_cnv_len,
                                        as_dev=self.as_dev, dev_params=self.dev_params, debug_dir=self.debug_dir)

    def coverage_dir_files(self, coverage_dir):
        coverage_dir = os.path.abspath(os.path.expanduser(coverage_dir))
        coverage_file = ""
        coverage_summary_file = ""

        for each_dir in os.listdir(coverage_dir):
            if "mosdepth.summary.txt" in each_dir:
                coverage_summary_file = f'{coverage_dir}/{each_dir}'
            elif ".regions.bed.gz" in each_dir and ("csi" not in each_dir and "tbi" not in each_dir):
                coverage_file = f'{coverage_dir}/{each_dir}'
            else:
                pass
        if coverage_file != "" and coverage_summary_file != "":
            return [coverage_file, coverage_summary_file]
        else:
            self.logger.error(f'coverage file or summary not found, directory {coverage_dir} has the following files:\n'
                              f'  {os.listdir(coverage_dir)}')
        return ["", ""]

    def __get_config(self, coverage_dir):
        self.logger.info(f'Spectre calculating for: {str(coverage_dir)} and bin size: {self.bin_size}')
        # input coverage
        [coverage_file, cov_summary_file] = self.coverage_dir_files(coverage_dir)
        self.mosdepth_data = spectre.util.mosdepthReader.MosdepthReader(coverage_file, cov_summary_file)
        self.mosdepth_data.summary_data()
        # output bed dir/file
        self.out_bed = os.path.join(os.path.join(self.out_dir, f'{self.sample_id}_cnv.bed'))
        # output VCF dir/file
        self.out_vcf = os.path.join(os.path.join(self.out_dir, f'{self.sample_id}_cnv.vcf'))

    def __get_config_dict(self, coverage_dir) -> dict:
        # get basic mosdepth coverage for the provided population sample
        result = {}
        self.logger.info(f"Spectre calculating population sample for: {str(coverage_dir)} and bin size {self.bin_size}")
        [coverage_file, cov_summary_file] = self.coverage_dir_files(coverage_dir)
        result["mosdepth_data"] = spectre.util.mosdepthReader.MosdepthReader(coverage_file, cov_summary_file)
        result["mosdepth_data"].summary_data()  # calculate summary
        return result

    def __get_population_config(self, population_coverage_dirs):
        # get config for the provided population samples
        for population_coverage_dir in population_coverage_dirs:
            population_sample = os.path.join(population_coverage_dir).split("/")[-2]
            self.logger.info(f"{population_sample}\t{population_coverage_dir}")
            self.population[population_sample] = self.__get_config_dict(population_coverage_dir)

    def cnv_call(self):
        # Data normalization
        self.logger.info("Data normalization and outlier removal (right tail)")
        self.cnv_analysis.data_normalization()

        # Coverage analysis
        self.logger.info(f"CNV calling - Coverage for sample: {self.sample_id}")
        # get raw CNV for original sample
        self.cnv_analysis.call_cnv_coverage(write_csv=self.as_dev)
        self.cnv_analysis.get_cnv_metrics()
        # self.cnv_analysis.write_intermediate_candidates("raw")
        # refine cnvs
        self.cnv_analysis.refine_cnv_calls(self.as_dev)  # set to self.as_dev

        # SNV analysis
        if self.snv_file_vcf != "":
            self.logger.info("CNV candidates by SNV")
            snv_file_basename_no_dot = "_".join(os.path.basename(self.snv_file_vcf).split('.'))
            self.snv_file_bed_af = f'{self.out_dir}/{snv_file_basename_no_dot}.bed'
            self.cnv_analysis.convert_vcf_to_tabular(self.snv_file_bed_af)
            self.cnv_analysis.call_cnv_af_region()

        # CNV metrics
        # self.logger.warning("Disabled CNV metrics")
        self.logger.info("Calculate CNV metrics")
        self.cnv_analysis.get_cnv_metrics(refined_cnvs=True)

        # Make output files
        self.logger.info("Final candidates are written to spc file")
        self.cnv_analysis.write_intermediate_candidates()
        self.logger.info("Results are writen to bed file")
        self.cnv_analysis.cnv_result_bed()
        self.logger.info("Results are writen to VCF file")
        self.cnv_analysis.cnv_result_vcf()
        self.logger.info("Result plot in progress")
        self.cnv_analysis.cnv_plot()

        # End
        self.logger.info(f"Output dir: {self.out_dir}")
        self.logger.info(f"Done with sample {self.sample_id}")
        # sys.exit(0)
