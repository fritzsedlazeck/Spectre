import os
import gzip


class OSUtil:
    @classmethod
    def get_absolute_path(cls, relative_user_path) -> str:
        """
        Extends a relative path to an absolute user path
        :param relative_user_path: relative path
        :return: absolute user path
        """
        return os.path.abspath(os.path.expanduser(relative_user_path))

    @classmethod
    def get_lines_by_chromosome(cls, filepath):
        filepath = cls.get_absolute_path(filepath)
        lines_by_chromosome = {}
        file_handler_read = gzip.open(filepath, "rt") if "gz" in filepath else open(filepath, "r")
        current_chromosome = ""
        for line in file_handler_read:
            chromosome = line.rstrip("\n").split("\t")[0]
            if current_chromosome == "":
                current_chromosome = chromosome
                lines_by_chromosome[chromosome] = 1
            elif current_chromosome != chromosome:
                current_chromosome = chromosome
                lines_by_chromosome[chromosome] = 1
            else:
                lines_by_chromosome[chromosome] += 1
        file_handler_read.close()
        return lines_by_chromosome

    @classmethod
    def get_lines_of_file(cls, filepath):
        filepath = cls.get_absolute_path(filepath)
        file_handler_read = open(filepath, "r") if "gz" not in filepath else gzip.open(filepath, "rt")
        cnt = 0
        for line in file_handler_read:
            cnt += 1

        file_handler_read.close()
        return cnt