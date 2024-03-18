import gzip
import numpy as np


# class MosdepthSummary(object):
#     def __init__(self):
#         self.chromosomes = []
#         self.chr_mean_coverage = {}
#         self.genome_mean_coverage = None
#         self.genome_bases = None
#
#     def add_chromosome(self, new_chr):
#         self.chromosomes.append(new_chr)
#
#     def add_coverage(self, new_chr, new_cov):
#         self.chr_mean_coverage[new_chr] = float(new_cov)
#
#     def add_genome_summary(self, genome_mean_coverage, genome_bases):
#         self.genome_mean_coverage = float(genome_mean_coverage)
#         self.genome_bases = int(genome_bases)
#
#     def update_genome_mean(self):
#         self.genome_mean_coverage = np.mean([value for index, value in self.chr_mean_coverage.items()])


class MosdepthReader(object):
    def __init__(self, coverage_file):
        self.coverage_file = coverage_file
        #self.summary_file = summary_file
        # self.mosdepth_summary_data = MosdepthSummary()
        # self.min_chr_lengh = 1e6  # 1mb to remove alts

    def get_bin_size(self)->int:
        """
        Determines the binsize of the mosdepth coverage regions file
        :return:
        """
        line_cnt = 0 # need only first two lines to determine bin size
        pos_list = []
        coverage_file_handler = open(self.coverage_file, "r") if "gz" not in self.coverage_file \
            else gzip.open(self.coverage_file, "rt")
        for line in coverage_file_handler:
            # skip header
            if line_cnt > 1:
                # chrom  start  end  depth
                break
            else:
                [_, start, _, _] = line.rstrip("\n").split("\t")
                pos_list.append(int(start))
            line_cnt += 1
        coverage_file_handler.close()
        return int(abs(pos_list[0] - pos_list[1]))

    # def summary_data(self):
    #     found_total_flag = False
    #     summary_file_handler = open(self.summary_file, "r") if "gz" not in self.summary_file \
    #         else gzip.open(self.summary_file, "rt")
    #     for line in summary_file_handler:
    #         # skip header
    #         if "chrom" not in line:
    #             if "_region" not in line:
    #                 # chrom  length  bases  mean  min  max
    #                 [chromosome, chr_length, bases, mean_cov, _, _] = line.rstrip("\n").split("\t")
    #                 if float(chr_length) > self.min_chr_lengh:
    #                     if chromosome == "total":
    #                         self.mosdepth_summary_data.add_genome_summary(mean_cov, bases)
    #                         found_total_flag = True
    #                     else:
    #                         self.mosdepth_summary_data.add_chromosome(chromosome)
    #                         self.mosdepth_summary_data.add_coverage(chromosome, mean_cov)
    #     if not found_total_flag:
    #         self.mosdepth_summary_data.update_genome_mean()
    #     summary_file_handler.close()
