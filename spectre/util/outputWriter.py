import gzip
import string
import random
import numpy as np
import json
from collections import Counter
import pysam
import os
from spectre.util import logger
import pysam.bcftools as bcftools


class BedOutput(object):
    def __init__(self, output_file):
        self.output_bed = output_file

    def make_bed(self, chromosome_list, cnv_candidate_list):

        if "gz" in self.output_bed:
            tmp_out = self.output_bed[:-3]  # remove the .gz
        file_handler = open(self.output_bed, "w") if "gz" not in self.output_bed else open(tmp_out, "w")
        for each_chromosome in chromosome_list:
            # sort by start position
            sort_list = sorted(cnv_candidate_list[each_chromosome], key=lambda x: x.start)
            # for each_candidate in cnv_candidate_list[each_chromosome]:
            for each_candidate in sort_list:
                result_list = [each_candidate.chromosome, str(each_candidate.start), str(each_candidate.end),
                               each_candidate.type, str(each_candidate.size),
                               str(round(each_candidate.median_cov_norm, 2))]
                bed_line = "\t".join(result_list)
                file_handler.write(f'{bed_line}\n')
        file_handler.close()

        # compress and clean if needed
        if "gz" in self.output_bed:
            pysam.tabix_compress(filename_in=tmp_out, filename_out=self.output_bed, force=True)
            pysam.tabix_index(filename=self.output_bed, force=True, preset="bed")
            if os.path.isfile(self.output_bed) and os.stat(self.output_bed).st_size != 0:
                os.remove(tmp_out)


class VCFLine(object):
    def __init__(self):
        self.CHROM = ""
        self.POS = 0
        self.ID = "."
        self.REF = "N"
        self.ALT = "."  # DEL/DUP
        self.QUAL = "."
        self.FILTER = "."
        self.INFO = "."
        self.FORMAT = "GT:HO:GQ:DP"
        self.format_data = []
        self.sample_format_data = {}
        self.supp_vec = {}

    def format_vcf_line(self):
        sample_format_list = []
        if len(self.format_data) > 0:
            sample_format_list = [":".join(self.format_data)]

        for key, value in self.sample_format_data.items():
            sample_format_list.append(":".join([str(s) for s in value]))

        return "\t".join([self.CHROM, str(self.POS), self.ID, self.REF, self.ALT, self.QUAL, self.FILTER,
                          self.INFO, self.FORMAT] + sample_format_list)


class VCFOutput(object):
    def __init__(self, output_file, genome_info):
        self.supp_vec = {}
        self.output_vcf = output_file
        self.population_sample_ids = []
        self.genome_info = genome_info
        # example {'chromosomes': ['chr6'], 'chr_lengths': {'chr6': 170805979}}

    @staticmethod
    def id_generator(size=8, chars=string.ascii_uppercase + string.digits):
        return ''.join(random.choice(chars) for _ in range(size))

    def make_vcf_contigs(self):
        vcf_contigs = []
        for contig_name in self.genome_info["chromosomes"]:
            contig_length = self.genome_info["chr_lengths"][contig_name]
            vcf_contigs.append(f'##contig=<ID={contig_name},length={contig_length}>')
        return "\n".join(vcf_contigs)

    def make_vcf_header(self):
        population_mode = False
        if self.population_sample_ids:
            if len(self.population_sample_ids) > 1:
                population_mode = True

        contigs = self.make_vcf_contigs()
        vcf_header = ['##fileformat=VCFv4.2', '##FILTER=<ID=PASS,Description="All filters passed">',
                      '##source=Spectre', contigs,
                      '##ALT=<ID=DEL,Description="Deletion">',
                      '##ALT=<ID=DUP,Description="Duplications">',
                      '##ALT=<ID=LOH,Description="Loss of Heterozygosity,">',
                      '##FILTER=<ID=UNRESOLVED,Description="An insertion that is longer than the '
                      'read and thus we cannot predict the full size.">'
                      ]
        # Add Info section
        vcf_header += ['##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">',
                       '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Length of the SV">',
                       '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of copy number variant">',
                       '##INFO=<ID=CN,Number=1,Type=Integer,Description="Estimated copy number status">',
                       '##INFO=<ID=SVSUPPORT,Number=1,Type=String,Description="SV support">']
        if population_mode:
            vcf_header += ['##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="Support vector">']

        # Add Format section
        vcf_header += ['##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                       '##FORMAT=<ID=HO,Number=2,Type=Float,Description="Homozygosity proportion">',
                       '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">',
                       '##FORMAT=<ID=DP,Number=2,Type=Float,Description="Read depth">',
                       '##FORMAT=<ID=ID,Number=1,Type=String,Description="Population ID of supporting CNV calls">'
                       ]

        if self.population_sample_ids:
            if population_mode:
                s = f"##Spectre_population_samples={','.join(self.population_sample_ids)}"
            else:
                s = f"##Spectre_sample={','.join(self.population_sample_ids)}"

            vcf_header.append(s)
        return "\n".join(vcf_header)

    @staticmethod
    def make_vcf_sample_header(sample_id: list):
        return "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + sample_id)

    @staticmethod
    def set_svtype(candidate_type):
        sv_types = {"DEL": "<DEL>", "DUP": "<DUP>", "LOH": "<LOH>"}
        if candidate_type in sv_types:
            return sv_types[candidate_type]
        else:
            return "."

    def vcf_result(self, chromosome_list, cnv_candidate_list):
        vcf_lines = []
        for each_chromosome in chromosome_list:
            for each_candidate in cnv_candidate_list[each_chromosome]:
                vcf_line = VCFLine()
                vcf_line.CHROM = each_candidate.chromosome
                vcf_line.POS = each_candidate.start
                vcf_line.ID = each_candidate.id
                vcf_line.ALT = self.set_svtype(each_candidate.type)
                vcf_line.FILTER = each_candidate.filter
                vcf_line.QUAL = str(each_candidate.quality)
                vcf_line.INFO = f"END={each_candidate.end};SVLEN={each_candidate.size};SVTYPE={each_candidate.type};" \
                                f"CN={each_candidate.cn_status};"
                vcf_line.INFO += f"SVSUPPORT={'TRUE' if each_candidate.sv_support else 'FALSE'};"

                # checking if any other CNV through merging supported the given CNV
                vcf_line.supp_vec = self.supp_vec.copy()
                if not each_candidate.support_cnv_calls:
                    vcf_line.format_data = [each_candidate.gt, f'{round(each_candidate.het_score, 2)}',
                                            f"{int(each_candidate.statistics['z-score']['sample_score'])}",
                                            f'{round(each_candidate.median_raw_cov,2)}']
                else:
                    vcf_line.format_data = []
                    vcf_line.FORMAT += ":ID"  # ADD ID tag in format as everything that is following are IDs
                    vcf_line.supp_vec = dict(
                        [(str(sample_key), 0) for sample_key in each_candidate.support_cnv_calls.keys()])
                    for sample_key, sample_value in each_candidate.support_cnv_calls.items():
                        if sample_value:
                            ids = []
                            scores = []
                            gts = []
                            dps = [] # raw read depth
                            # sample_id =""
                            for candidate in list(sample_value):
                                ids.append(candidate.id)
                                scores.append(candidate.statistics['z-score']['sample_score'])
                                gts.append(candidate.gt)
                                # check if median_raw_cov exists
                                if candidate.median_raw_cov:
                                    dps.append(candidate.median_raw_cov)
                            gt = Counter(gts).most_common(1)[0][0]
                            score = np.mean(scores)
                            dp = np.mean(dps)
                            vcf_line.sample_format_data[sample_key] = [gt, "0.0", int(score),round(dp,2), ",".join(ids)]
                            vcf_line.supp_vec[sample_key] = 1
                        else:
                            vcf_line.sample_format_data[sample_key] = ["0/0", "0.0", 0, 0.0,"NULL"]
                    # add support vector only if population mode is active
                    vcf_line.INFO += ";SUPP_VEC=" + "".join([str(s) for s in vcf_line.supp_vec.values()])
                vcf_lines.append(vcf_line.format_vcf_line())
        return "\n".join(vcf_lines)

    def make_vcf(self, chromosome_list, cnv_candidate_list, sample_id, population_sample_ids=None):
        # converting population sample ids from set to list
        if not population_sample_ids:
            population_sample_ids = [sample_id]
        self.population_sample_ids = list(population_sample_ids)
        self.supp_vec = dict([(i, 0) for i in self.population_sample_ids])

        # we need to update to pysam.BGZFile(out_file, "w")
        if "gz" in self.output_vcf:
            tmp_out = self.output_vcf[:-3]  # remove the .gz
        file_handler = open(self.output_vcf, "w") if "gz" not in self.output_vcf else open(tmp_out, "w")
        vcf_header = self.make_vcf_header()
        vcf_sample_header = self.make_vcf_sample_header(self.population_sample_ids)
        vcf_lines = self.vcf_result(chromosome_list, cnv_candidate_list)
        file_handler.write(f'{vcf_header}\n')
        file_handler.write(f'{vcf_sample_header}\n')
        file_handler.write(f'{vcf_lines}\n')
        file_handler.close()
        if "gz" in self.output_vcf and len(vcf_lines) > 0:
            # sort vcf file
            bcftools.sort(f"-o", tmp_out + ".sort.tmp", tmp_out, catch_stdout=False)

            # compress and index vcf file
            pysam.tabix_compress(filename_in=tmp_out + ".sort.tmp", filename_out=self.output_vcf, force=True)
            pysam.tabix_index(filename=self.output_vcf, force=True, preset="vcf")

        # clean up
        if os.path.isfile(self.output_vcf) and os.stat(self.output_vcf).st_size != 0:
            os.remove(tmp_out) # remove the original unsorted vcf
            # sorting file
            if os.path.exists(tmp_out + ".sort.tmp"):
                os.remove(tmp_out + ".sort.tmp")

class IntermediateFile(object):
    def __init__(self, output_dir: str, as_dev=False):
        self.output_dir = output_dir
        self.logger = logger.setup_log(__name__, as_dev)

    @staticmethod
    def convert_genome_info_to_dictionary(genome_info: dict):
        tmp_genome_inf = dict([(key, []) for key in genome_info.keys()])
        return ""
        # for info in genome_info;

    @staticmethod
    def convert_candidates_to_dictionary(candidates: dict):
        tmp_candidates_dict = dict([(key, []) for key in candidates.keys()])
        for key, candidates in candidates.items():
            tmp_candidates = []
            tmp_candidates_dict[key] = tmp_candidates
            for candidate in candidates:
                tmp_dict = dict(candidate.__dict__)
                tmp_dict = {k: v for k, v in tmp_dict.items() if v}
                for key, value in tmp_dict.items():
                    if isinstance(value, set):
                        tmp_dict[key] = list(value)
                #tmp_dict.pop('logger')
                tmp_candidates.append(tmp_dict)
        return tmp_candidates_dict

    def write_intermediate_file(self, output_object, filename: str) -> str:
        output_path = f"{self.output_dir}/{filename}.spc.gz"
        with gzip.open(output_path, 'wt', encoding='UTF-8') as out_file:
            try:
                json.dump(output_object, out_file, indent="\t")
            except TypeError:
                self.logger.error("write_intermediate_file failed, json.dump error")
        return output_path
