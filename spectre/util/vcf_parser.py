import pandas as pd
import re
from spectre.util import logger
import os
import pysam
import gzip
from spectre.classes.cnv_candidate import CNVCandidate


class VCFParser(object):
    def __init__(self, only_chr_list=None, min_chromosome_len=1e6, as_dev=False):
        self.as_dev = as_dev
        self.logger = logger.setup_log(__name__, self.as_dev)
        self.only_chr_list = only_chr_list
        self.min_chromosome_len = min_chromosome_len

    @staticmethod
    def __get_absolute_path(input_path):
        return os.path.abspath(os.path.expanduser(input_path))

    def __create_tabix_file(self, input_path):
        self.logger.info(f"No index file found - Creating new one at {input_path}.tbi")
        pysam.tabix_index(input_path, preset="bed", force=True)

    def __create_vcf_tabix_file(self, input_path):
        self.logger.info(f"No index file found - Creating new one at {input_path}.tbi")
        pysam.tabix_index(input_path, preset="vcf", force=True, keep_original=True)

    def vcf_has_index(self, path) -> str:
        absolute_path = self.__get_absolute_path(path)

        # check if index file exists
        if not os.path.isfile(absolute_path + ".tbi"):
            self.logger.info(f"No index file found for {absolute_path}")
            self.logger.info(f"Creating new one at {absolute_path}.tbi")
            if absolute_path.endswith(".vcf"):
                self.logger.warning(f"Compressed VCF file is required for indexing")
                self.logger.warning(f"Spectre will do that for you ... this will take a while")
            self.__create_vcf_tabix_file(path)
            # was not compressed in the original version
            if absolute_path.endswith(".vcf"):
                self.logger.info(f"Compressing VCF file to {absolute_path}.gz")
                return absolute_path + ".gz"
        return absolute_path


    def vcf_pysam(self, path):
        # setup
        vcf_path = self.__get_absolute_path(path)
        vcf_tabix_path = vcf_path + ".tbi"

        # checking for index file
        if not os.path.isfile(vcf_tabix_path):
            self.__create_tabix_file(vcf_path)

        vcf_file = pysam.VariantFile(vcf_path)  # loading vcf file
        records = vcf_file.fetch()
        for record in records:
            x = record.samples.values()
            print()

        for uid, x in enumerate(vcf_file.fetch()):
            print(x.format.keys())
            af = x.format.values()[4].number

    def vcf_to_dataframe(self, path):
        self.logger.debug("Converting vcf to dataframe")
        vcf_path = self.__get_absolute_path(path)

        # pre allocation
        vcf_entries = list()
        vcf_file = open(vcf_path, "r") if "gz" not in vcf_path else gzip.open(vcf_path, "rt")
        chrom_, start_, quality_, filter_, info_, format_, sample_ = "", "", "", "", "", "", ""
        # how AF is described by tool:
        #   clair3   -> AF in FORMAT/SAMPLE
        #   longshot -> AF not present, AC contains the number of reads for REF,ALT: AF=ALT/(REF+ALT)
        #   p-m-deep -> AF not present, VAF in FORMAT/SAMPLE
        #  if in info --- af = float(list(filter(lambda x: "VAF" in x, info_.split(";")))[0].split("=")[1])
        check_chr_len = False
        if self.only_chr_list is None:
            self.only_chr_list = []
            check_chr_len = True
        for line in vcf_file:
            line = line.rstrip("\n")
            # skip header
            if line.startswith("#"):
                if check_chr_len:
                    if line.__contains__("contig"):
                        # example:  ##contig=<ID=chr1_KI270706v1_random,length=175055>
                        try:
                            chr_id_len = re.search('.*<ID=(.*),length=(.*)>', line)
                            [chr_id, chr_len] = [chr_id_len.group(1), int(chr_id_len.group(2))]
                            if chr_len >= self.min_chromosome_len:
                                self.only_chr_list.append(chr_id)
                        except AttributeError:
                            self.logger.debug(line)
            else:
                [chrom_, start_, _, _, _, quality_, filter_, info_, format_, sample_] = line.split("\t")
                if chrom_ in self.only_chr_list:
                    # searching in INFO for AF
                    if format_.__contains__("AF"):
                        # clair3 or pepper-margin-deepvariant
                        try:
                            # clair3
                            format_.split(":").index("AF")
                            af_index = format_.split(":").index("AF")
                            af = float(sample_.split(":")[af_index])
                            # self.logger.debug("like clair3")
                        except ValueError:
                            if format_.__contains__("VAF"):
                                # pepper-margin-deepvariant
                                af_index = format_.split(":").index("VAF")
                                af = float(sample_.split(":")[af_index])
                                # self.logger.debug("like pepper-margin-deepvariant")
                            else:
                                pass
                    elif info_.__contains__("AC"):
                        # longshot
                        ac = list(filter(lambda x: "AC" in x, info_.split(";")))[0].split("=")[1]
                        [ref_count, alt_count] = ac.split(",")
                        af = float(alt_count) / (float(ref_count) + float(alt_count))
                        # self.logger.debug("like longshot")
                    else:
                        # TODO use: DR/DV for calculating
                        af = "NA"
                    vcf_entries.append([chrom_, int(start_), None, af])
        vcf_file.close()
        return pd.DataFrame(data=vcf_entries, columns=["chrom_", "start_", "end_", "af_"])

    def get_mosdepth_chromosome_borders(self, mosdepth_file: str = ""):
        in_path = self.__get_absolute_path(mosdepth_file)
        file = pd.read_csv(in_path, sep="\t", header=None, names=["chrom_", "start_", "end_", "af_"])
        return file

    # def dataframe_to_tabular_file(self, df_snv: pd.DataFrame = None, mosdepth_input: str = "", out_path: str = ""):
    #     self.logger.debug("Writing dataframe to tabular")
    #     df_mosdepth = self.get_mosdepth_chromosome_borders(mosdepth_input)
    #     df_mosdepth_grouped = df_mosdepth.groupby("chrom_")
    #     df_snps_grouped = df_snv.groupby("chrom_")
    #     df_final = pd.DataFrame()
    #     for mosdepth_chrom_key in df_mosdepth_grouped.indices.keys():
    #         if mosdepth_chrom_key in df_snps_grouped.indices.keys():
    #             df_mosdepth_chrom = df_mosdepth.loc[df_mosdepth_grouped.indices[mosdepth_chrom_key]]
    #             bins = list(df_mosdepth_chrom.start_)
    #             labels = list(df_mosdepth_chrom.start_)[:-1]  # must be len(bins)-1
    #
    #             df_snps_chrom = df_snv.iloc[df_snps_grouped.indices[mosdepth_chrom_key]]
    #             df_snps_chrom["startbin_"] = pd.cut(x=df_snps_chrom.start_, bins=bins, labels=labels,
    #                                                 include_lowest=False)
    #             df_snps_chrom["startbin_af_"] = df_snps_chrom.groupby("startbin_")["af_"].transform('mean')
    #             df_snps_chrom.sort_values(by=["startbin_"], inplace=True)
    #             df_snps_chrom.drop_duplicates("startbin_", inplace=True)
    #
    #             df_snps_chrom.drop(["start_", "end_", "af_"], axis=1, inplace=True)  # clean
    #
    #             df_merged = pd.merge(df_mosdepth_chrom, df_snps_chrom,
    #                                  left_on=["chrom_", "start_"],
    #                                  right_on=["chrom_", "startbin_"],
    #                                  how="left")
    #
    #             df_merged["startbin_"] = df_merged["start_"] + 1
    #             df_merged["startend_"] = df_merged["end_"]
    #             df_merged["startbin_af_"].fillna(value=0, inplace=True)
    #
    #             df_final = pd.concat([df_final, df_merged[["chrom_", "startbin_", "startend_", "startbin_af_"]]],
    #                                  ignore_index=True)
    #             # df_final = pd.concat([df_final,df_snps_chrom],ignore_index=True)
    #
    #     df_final.to_csv(f'{out_path}', sep="\t", index=False, header=None)
    #     return df_final


class VCFtoCandidate(object):
    def __init__(self):
        pass

    def vcf_to_candidates(self, vcf_path):
        df = self.vcf_ot_dataframe(vcf_path)
        cnv_candidate_list = self.dataframe_to_candidates(df)
        return cnv_candidate_list

    def vcf_ot_dataframe(self, vcf_path: str = ''):
        # Loading VCF file
        vcf_file = open(vcf_path, "r") if "gz" not in vcf_path else gzip.open(vcf_path, "rt")
        lines = [line.strip() for line in vcf_file.readlines()]
        vcf_file.close()
        # Creating dataframe
        df = pd.DataFrame(lines, columns=['input'])
        df = df.input.str.split('\t', expand=True)
        df = df[~df.iloc[:, 0].str.contains('##')]  # header definitions starting with ##
        df.iloc[0, 0] = df.iloc[0, 0][1:]  # first line is header line of vcf entries removing the # in the first col
        df.columns = df.iloc[0]
        df = df[1:]  # removing first row which is used for the column names
        return df

    def dataframe_to_candidates(self, df: pd.DataFrame):
        cnv_candidate_list = {}
        for cnv in df.itertuples():
            # Parse basics
            candidate = CNVCandidate()
            candidate.chromosome = cnv.CHROM
            candidate.start = int(cnv.POS)
            candidate.id = cnv.ID
            info_dict = {i.split('=')[0]: i.split('=')[1] for i in cnv.INFO.split(';')}
            candidate.end = int(info_dict['END'])
            candidate.cn_status = int(info_dict['CN'])
            candidate.type = info_dict['SVTYPE']
            candidate.size = int(info_dict['SVLEN'])

            population_mode = len(df.columns[9:]) < 2  # TRUE more than 2 samples are present in the VCF
            # Form
            format_set = [i for i in cnv.FORMAT.split(':')]
            for sample_id, sample in zip(list(df.columns[9:]), list(cnv[9 + 1:])):  # +1 due to the index column
                sample_gq = 0
                if not population_mode:  # without special flags it is not possible
                    candidate.sample_origin = sample_id
                sample_cells = sample.split(':')
                if "GQ" in format_set:
                    gq_idx = format_set.index('GQ')  # get gene quality score
                    sample_gq = int(sample_cells[gq_idx])
                    candidate.statistics['z-score'] = {}
                    candidate.statistics['z-score']['sample_score'] = sample_gq
                if "GT" in format_set:
                    gt_idx = format_set.index('GT')
                    candidate.gt = sample_cells[gt_idx]

                # Could be left black as only details for the known variant are known and loaded from VCF
                if "ID" in format_set:
                    id_idx = format_set.index('ID')
                    support_samples = sample_cells[id_idx].split(',')
                    # Get all supporting cnvs from a given sample
                    for support_sample_id in support_samples:
                        # add only not NULL supports
                        if support_sample_id != 'NULL':
                            if sample_id not in candidate.support_cnv_calls.keys():
                                candidate.support_cnv_calls[sample_id] = set()
                            support_candidate = CNVCandidate(sample_origin=sample_id, bin_size=candidate.bin_size)
                            support_candidate.id = support_sample_id
                            support_candidate.statistics['z-score'] = {}
                            support_candidate.statistics['z-score']['sample_score'] = sample_gq
                            candidate.support_cnv_calls[sample_id].add(support_candidate)

            if cnv.CHROM not in cnv_candidate_list.keys():
                cnv_candidate_list[cnv.CHROM] = []
            cnv_candidate_list[cnv.CHROM].append(candidate)
        return cnv_candidate_list
