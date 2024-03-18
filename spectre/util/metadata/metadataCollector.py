import mmap
import os
import gzip
from spectre.util import logger


class FastaRef:
    def __init__(self, metadata_args=None, as_dev=False):
        # logger
        self.logger = logger.setup_log(__name__, as_dev)
        # ... TODO
        self.__nList = {}
        self.__nFrom = 0
        self.__nTo = 0
        self.__gCnt = 0
        self.__cCnt = 0
        self.__aCnt = 0
        self.__tCnt = 0
        # args
        self.metadata_args = metadata_args
        self.output_file = ""
        # report
        self.static_report = ""

    def get_n_regions(self, out_file_name) -> None:
        """
        Extracting exact positions of N regions from reference genome
        :param out_file_name:
        :return:
        """
        ref_dict = dict()
        i = 0
        j = 0
        g_cnt = 0
        c_cnt = 0
        a_cnt = 0
        t_cnt = 0
        chromosome_cnt = 0
        chromosome = ""
        # bin_size = self.metadata_args.bin_size
        self.output_file = out_file_name
        # assure abs paths
        reference = os.path.abspath(os.path.expanduser(self.metadata_args.reference))
        # DEL report_output_dir = os.path.abspath(os.path.expanduser(report_output_dir))
        self.metadata_args.out_dir
        file_handler_fasta = open(reference, "r") if "gz" not in reference else gzip.open(reference, "rt")
        for raw_line in file_handler_fasta:
            line = raw_line.rstrip("\n")  # cut away the newline from line (\n)
            # check if line is not a comment or header
            if not (line[:1] == '>' or line[:1] == ';'):
                for l in line:
                    char = l.upper()
                    if char == "N":
                        j += 1
                    else:
                        # No not in N region anymore
                        if (j - i) >= self.metadata_args.n_size:
                            start = i + 1
                            end = j - 1
                            # check if values are not on mod binSize = 0 position
                            # if (start % bin_size) != 0:
                            #     start = start - (start % bin_size)
                            # if (end % bin_size) != 0:
                            #     end = end + (bin_size - (end % bin_size))
                            ref_dict[str(chromosome)].append((start, end))
                            self.logger.info(f"Found No {len(ref_dict[str(chromosome)])} N region: {chromosome}:{start}-{end}")
                        i = j
                        j += 1
                        # grab bp for gc statistic
                        if char == "G":
                            g_cnt += 1
                        elif char == "C":
                            c_cnt += 1
                        elif char == "A":
                            a_cnt += 1
                        elif char == "T":
                            t_cnt += 1
                    chromosome_cnt += 1
            else:
                # allocate space in dict for current chromosome
                if chromosome != "":
                    # pre-print before new chrom is select
                    self.logger.info(f'length: {str(chromosome_cnt)}\t{ref_dict[str(chromosome)]}')
                chromosome = str(line[1:]).split(" ")[0]
                ref_dict[str(chromosome)] = list()
                self.logger.info(f'Reading {line}')
                i = 0
                j = 0
                chromosome_cnt = 0

        self.__nList = ref_dict
        self.__nTo = j
        self.__cCnt = c_cnt
        self.__gCnt = g_cnt
        self.__aCnt = a_cnt
        self.__tCnt = t_cnt
        self.__get_statistic_report()
        # if not self.metadata_args.save_only:
        #    return ref_dict
        # return {}

    def extact_metadata_from_mdr(self, input_file: str, bin_size: int) -> dict:
        """
        Adjust the true N region positions according to the bin size.
        :param input_file: MDR file
        :param bin_size: bin size determined by looking at the positional gap in the mosdepth data.
        :return: Adjusted medata N regions.
        """
        # read in the MDR file
        with open(input_file, "r") as file_handler:
            for line in file_handler:
                # skip header
                if line.startswith("NSPOS"):
                    # split line
                    split_line = line.split("\t")
                    chromosome = split_line[1]
                    start = int(split_line[2])
                    end = int(split_line[3])
                    # check if values are not on mod binSize = 0 position
                    if (start % bin_size) != 0:
                        start = start - (start % bin_size)
                    if (end % bin_size) != 0:
                        end = end + (bin_size - (end % bin_size))
                    # add to dict
                    if self.__nList.get(chromosome) is None:
                        self.__nList[chromosome] = list()
                    self.__nList[chromosome].append((start, end))
        return self.__nList

    def load_metadata(self, mdr_file_path: str = "", blacklist_file_path: str="", bin_size: int=1000) -> dict:
        """
        Load metadata for each sample separate from mdr and blacklist.
        Positions of N regions in the MDR are adjusted in such way that the start position is the next lower position
        which aligns with the bin size. The end position is the next higher position which aligns with the bin size.
        :param mdr_file_path:
        :param blacklist_file_path:
        :param bin_size:
        :return:
        """
        # TODO call from spectre_exe
        #  1) extract from mdr
        mdr_content = self.extact_metadata_from_mdr(mdr_file_path, bin_size)
        #  2) extract from blacklist
        blacklist_content = self.extract_blacklisted_regions(blacklist_file_path)
        #  3) merge
        metadata_result = self.merge_metadata(mdr_content, blacklist_content)
        return metadata_result

    def __get_statistic_report(self):
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
        self.__write_report()

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

    def extract_blacklisted_regions(self, blacklist_file: str = "") -> dict:
        result = dict()
        if blacklist_file.__eq__(""):
            path = os.path.abspath(os.path.expanduser(self.metadata_args.black_list))
        else:
            path = os.path.abspath(os.path.expanduser(blacklist_file))
        with open(path, "r") as fp:
            mm = mmap.mmap(fp.fileno(), 0, prot=mmap.PROT_READ)
            for line in iter(mm.readline, b""):
                # decode line and cut away last char from line (\n)
                term = line.decode("utf-8")[:-1]
                chrom, start, end = term.split("\t")[:3]
                if str(chrom) not in result:
                    result[str(chrom)] = []
                result[str(chrom)].append((start, end))
        return result

    @classmethod
    def merge_metadata(cls, m1: dict, m2: dict) -> dict:
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

    def __write_report(self):
        self.logger.info(f"Writing MDR to: {os.path.join(self.metadata_args.out_dir, self.output_file)}")
        if not os.path.exists(self.metadata_args.out_dir):
            os.makedirs(self.metadata_args.out_dir)
        with open(os.path.join(self.metadata_args.out_dir, self.output_file), "w") as file:
            file.writelines(self.static_report)
