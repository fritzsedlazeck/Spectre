import json
import gzip
import os
import pandas as pd
from spectre.util import logger


class Breakpoints(object):
    def __init__(self, breakpoints_file, as_dev=False):
        self.as_dev = as_dev
        self.logger = logger.setup_log(__name__, self.as_dev)
        self.breakpoints_file = breakpoints_file
        self.breakpoints = self.load_breakpoints()

    def load_breakpoints(self):
        if os.path.exists(self.breakpoints_file):
            if self.breakpoints_file.endswith(".gz"):
                with gzip.open(self.breakpoints_file, 'rt') as f:
                    return json.load(f)
            else:
                with open(self.breakpoints_file, 'r') as f:
                    return json.load(f)
        else:
            self.logger.error(f"Breakpoints file not found: {self.breakpoints_file}")
            return None

    def correlate_cnvs_with_breakpoints(self, cnv_candidates: dict, bin_size=1000):
        # get matching chromsomes
        cnv_candidates_chroms = list(cnv_candidates.keys())
        breakpoint_chroms = list(self.breakpoints.keys())
        # get matching chromosomes
        matching_chroms = list(set(cnv_candidates_chroms) & set(breakpoint_chroms))
        if len(matching_chroms) < 1:
            self.logger.warning(f"No matching chromosome names found between the Mosdepth data and the SNFJ file.")
            self.logger.warning(f"Check if the chromosome names are consistent between the files.")
            self.logger.warning("Skipping lookup of CNVs with breakpoints!")
        else:
            self.logger.info(f"Matching chromosome names found: {matching_chroms}")
            # filter breakpoints by matching chromosomes
            # Filter only for DEL and DUP data
            matching_breakpoints = {chrom: {bp_type: self.breakpoints[chrom][bp_type] for bp_type in ['DEL', 'DUP']} for
                                    chrom in matching_chroms}
            # Combine the entries for all available matching_breakpoints.keys according to their bp_type
            df_del_breakpoints = pd.DataFrame()
            df_dup_breakpoints = pd.DataFrame()
            for contig, value in matching_breakpoints.items():
                del_breakpoints = value["DEL"]
                dup_breakpoints = value["DUP"]
                df_del_tmp = pd.DataFrame(del_breakpoints)
                df_dup_tmp = pd.DataFrame(dup_breakpoints)
                # append df_del_tmp top df_del_breakpoints
                df_del_breakpoints = pd.concat([df_del_breakpoints, df_del_tmp])
                df_dup_breakpoints = pd.concat([df_dup_breakpoints, df_dup_tmp])

            def __row_constrains(row, c_start, c_end, offset_value=1000):
                return (row['pos'] >= c_start - offset_value) and (row['pos'] <= c_start + offset_value) and \
                    (row['end'] >= c_end - offset_value) and (row['end'] <= c_end + offset_value)

            offset = 5 * bin_size
            # Look up supporting breakpoint entries at the end and beginning of and candidate and count their encounters
            for chrom, candidates in cnv_candidates.items():
                for candidate in candidates:
                    # filter the breakpoints based on svlen and svtype
                    df_tmp_breakpoints = pd.DataFrame()
                    if candidate.type == "DEL":
                        # Get the breakpoints that are within the range of the candidate size
                        df_tmp_breakpoints = df_del_breakpoints[(df_del_breakpoints['svlen'].abs() >= candidate.size - offset) & (df_del_breakpoints['svlen'].abs() <= candidate.size + offset)]
                    else:
                        df_tmp_breakpoints = df_dup_breakpoints[(df_dup_breakpoints['svlen'].abs() >= candidate.size - offset) & (df_dup_breakpoints['svlen'].abs() <= candidate.size + offset)]
                    # Calculate the number of breakpoints that are within a offset range of start and end of the candidate
                    matching_rows_len = len(df_tmp_breakpoints[
                        (df_tmp_breakpoints['pos'] >= candidate.start - offset) &
                        (df_tmp_breakpoints['pos'] <= candidate.start + offset) &
                        (df_tmp_breakpoints['end'] >= candidate.end - offset) &
                        (df_tmp_breakpoints['end'] <= candidate.end + offset )])
                    candidate.sv_support = False if matching_rows_len == 0 else True
