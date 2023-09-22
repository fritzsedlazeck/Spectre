import mmap
import os
import gzip
import logging as logger


class FastaRef:

    def __init__(self, as_dev=False):
        # logger
        logger.basicConfig(level=logger.DEBUG) if as_dev else logger.basicConfig(level=logger.INFO)
        self.logger = logger
        # ... TODO
        self.__nList = list()
        self.__nFrom = 0
        self.__nTo = 0
        self.__gCnt = 0
        self.__cCnt = 0
        self.__aCnt = 0
        self.__tCnt = 0
        # report
        self.static_report = ""

    def get_n_regions(self, filepath, report_output_dir, out_file_name, threshold=5, bin_size=500,
                      save_only=False) -> dict:
        ref_dict = dict()
        i = 0
        j = 0
        g_cnt = 0
        c_cnt = 0
        a_cnt = 0
        t_cnt = 0
        chromosome_cnt = 0
        chromosome = ""
        # assure abs paths
        path = os.path.abspath(os.path.expanduser(filepath))
        report_output_dir = os.path.abspath(os.path.expanduser(report_output_dir))
        file_handler_fasta = open(path, "r") if "gz" not in path else gzip.open(path, "rt")
        for line in file_handler_fasta:
            term = line.rstrip("\n")  # cut away the newline from line (\n)
            # check if line is not a comment or header
            if not (term[:1] == '>' or term[:1] == ';'):
                for t in term:
                    c = t.upper()
                    if c == "N":
                        j += 1
                    else:
                        if (j - i) >= threshold:
                            start = i + 1
                            end = j - 1
                            # check if values are not on mod binSize = 0 position
                            if (start % bin_size) != 0:
                                start = start - (start % bin_size)
                            if (end % bin_size) != 0:
                                end = end + (bin_size - (end % bin_size))
                            ref_dict[str(chromosome)].append((start, end))
                        i = j
                        j += 1
                        # grab bp for gc statistic
                        if c == "G":
                            g_cnt += 1
                        elif c == "C":
                            c_cnt += 1
                        elif c == "A":
                            a_cnt += 1
                        elif c == "T":
                            t_cnt += 1
                    chromosome_cnt += 1
            else:
                # allocate space in dict for current chromosome
                if chromosome != "":
                    # pre-print before new chrom is select
                    self.logger.info(f'length: {str(chromosome_cnt)}\t{ref_dict[str(chromosome)]}')
                chromosome = str(term[1:]).split(" ")[0]
                ref_dict[str(chromosome)] = list()
                self.logger.info(f'Reading {term}')
                i = 0
                j = 0
                chromosome_cnt = 0

        self.__nList = ref_dict
        self.__nTo = j
        self.__cCnt = c_cnt
        self.__gCnt = g_cnt
        self.__aCnt = a_cnt
        self.__tCnt = t_cnt
        self.__get_statistic_report(report_output_dir, out_file_name)
        if not save_only:
            return ref_dict
        return {}

    def __get_statistic_report(self, output_dir, filename):
        # static_report is a long string of text
        self.static_report = "#####\tStatistic report\n"
        # gc coverage
        self.logger.info("Calculating bp statistics")
        self.static_report += "#####\tBasepair statistic\n"
        self.static_report += "BPDEF\tA\tT\tC\tG\tGC%\tTotal bp\n"
        self.static_report += "BPSTA\t" + str(self.__aCnt) + "\t" + str(self.__tCnt) + "\t" + str(self.__cCnt) + \
                              "\t" + str(self.__gCnt) + "\t" + str((self.__gCnt + self.__cCnt) / self.__nTo * 100) + \
                              "\t" + str(self.__nTo) + "\n"

        # positions of n in sequence
        self.logger.info("Calculating N positions")
        self.static_report += "#####\tN-Sequence positions in given reference file\n"
        self.static_report += "NSDEF\tFrom\tTo\n"
        for chromosome in self.__nList:
            for start, end in self.__nList[chromosome]:
                self.static_report += "NSPOS\t" + str(chromosome) + "\t" + str(start) + "\t" + str(end) + "\n"

        # write statistic report to file
        self.logger.info("Writing report")
        self.__write_report(output_dir, filename)

    @classmethod
    def extract_n_regions_from_report(cls, input_file: str) -> dict:
        """
        Loads all N regions that previously have been extracted from the reference genome and stored in a .mdr file
        :param input_file: /path/to/<file_name>.mdr
        :return: dict with all N regions according to the chromosome
        """
        result = dict()
        # assure abs paths
        path = os.path.abspath(os.path.expanduser(input_file))
        with open(path, "r") as fp:
            # map file into memory
            mm = mmap.mmap(fp.fileno(), 0, prot=mmap.PROT_READ)
            for line in iter(mm.readline, b""):
                # decode line and cut away last char from line (\n)
                term = line.decode("utf-8")[:-1]
                if term[:5].__eq__("NSPOS"):
                    cells = term.strip().split("\t")
                    if str(cells[1]) not in result:  # [1] = chromosome name
                        result[str(cells[1])] = []
                    result[str(cells[1])].append((cells[2], cells[3]))  # [2] start, [3] end
        return result

    @classmethod
    def extract_blacklisted_regions(cls, input_file: str) -> dict:
        result = dict()
        path = os.path.abspath(os.path.expanduser(input_file))
        with open(path, "r") as fp:
            mm = mmap.mmap(fp.fileno(), 0, prot=mmap.PROT_READ)
            for line in iter(mm.readline, b""):
                # decode line and cut away last char from line (\n)
                term = line.decode("utf-8")[:-1]
                chrom, start, end = term.split("\t")[:3]
                if str(chrom) not in result:
                    result[str(chrom)] = []
                result[str(chrom)].append((start,end))
        return result

    @classmethod
    def merge_metadata(cls, m1:dict, m2:dict)->dict:
        """
        Merges two dictionaries with the structure str:list
        :param m1: dict m1 with structure str:list
        :param m2: dict m2 with structure str:list
        :return: merged dict with structure str:list
        """
        result = m1.copy()
        for key2 in m2.keys():
            # check if key in results
            if key2 in result:
                # check if values are in it
                for value2 in m2[key2]:
                    # add if value2 not in result
                    if value2 not in result[key2]:
                        result[key2].append(value2)
            else:
                # if not in results -> add whole d2[key2] into it
                result[key2] = m2[key2]
        return result

    def __write_report(self, output_dir, filename):
        self.logger.info(output_dir)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        with open(os.path.join(output_dir, filename), "w") as file:
            file.writelines(self.static_report)
