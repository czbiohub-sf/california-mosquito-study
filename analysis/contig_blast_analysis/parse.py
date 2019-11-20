#!/usr/bin/env python3
# script modified from https://github.com/czbiohub/pyblastc/blob/master/scripts/parse.py

import sys
import time
from collections import defaultdict

from param import ranked_blast_output_schema, blast_outfmt6_schema
import param
from util import tsv_rows, parse_headerless_table, tsprint

# goal: find the subset of HSPs that accumulate the total_bitscores best

def interval(a, b):
    return (min(a, b), max(a, b))

def intervals_overlap(p, q, fraction=param.MIN_OVERLAP):
    p_len = p[1] - p[0]
    q_len = q[1] - q[0]
    min_len = min(p_len, q_len)
    o_len = max(0.0,  min(p[1], q[1]) - max(p[0], q[0]))
    return (o_len / min_len) > fraction

def query_interval(row):
    # decode (start, end, strand) for query
    return interval(row["qstart"], row["qend"])

def subject_interval(row):
    # decode (start, end, strand) for query
    return interval(row["sstart"], row["send"])

def hsp_overlap(hsp_1, hsp_2):
    # let's worry about subject sequence overlap later
    return intervals_overlap(query_interval(hsp_1), query_interval(hsp_2))


def find(needle, haystack):
    # return True if needle overlaps any hay
    return any(hsp_overlap(needle, hay) for hay in haystack)

class Optimizer:

    def __init__(self, hsps):
        # List of HSPs from the same query to the same subject sequence,
        # ordered by decreasing bitscore.
        self.hsps = hsps
        self.optimal_set = None
        self.score = None

    def solve(self):
        # Find a subset of disjoint HSPs with maximum sum of bitscores
        # Initial implementation:  Super greedy.
        optimal_set = [self.hsps[0]]
        for next_hsp in self.hsps[1:]:
            if not find(next_hsp, optimal_set):
                optimal_set.append(next_hsp)
        self.score = sum(hsp["bitscore"] for hsp in optimal_set)
        self.optimal_set = optimal_set

    def solution_row(self):
        r = dict(self.optimal_set[0])
        r["hsplen"] = sum(hsp["hsplen"] for hsp in self.optimal_set)
        r["pident"] = sum(hsp["pident"] * hsp["hsplen"] for hsp in self.optimal_set) / r["hsplen"]
        r["bitscore"] = sum(hsp["bitscore"] for hsp in self.optimal_set)
        # these are new
        if "qlen" not in r:
            if '~' in r["qseqid"]:
                r["qlen"] = int(r["qseqid"].split("~")[1].split("_")[3])
            else:
                r["qlen"] = int(r["qseqid"].split("_")[3])
        r["qcov"] = r["hsplen"] / r["qlen"]
        r["hsp_count"] = len(self.optimal_set)
        return r


def parse_blastn_outfmt6(blast_output_fmt6):  # orig_m8
    """ process the output format 6 blast results """

    # 1 per-HSP filters
    # 2 group HSPs by (query id, subject sequence id)
    # 3 aggregate within each group to produce a score for the group
    # 4 emit top scoring group per query_id

    t_start = time.time()

    table_iterator = parse_headerless_table(tsv_rows(blast_output_fmt6), blast_outfmt6_schema)

    # Ouput column indices
    HSPs = defaultdict(list)

    for row in table_iterator:

        # local HSP sequence similarity filter
        if row["pident"] < param.MIN_PIDENT:
            continue

        # add HSP to group
        group_id = (row["qseqid"], row["sseqid"])

        HSPs[group_id].append(row)

    winners = dict()
    for group_id, hsps in HSPs.items():
        o = Optimizer(hsps)
        o.solve()
        query_id, _ = group_id
        sr = o.solution_row()
        print("\t".join(str(sr[k]) for k in ranked_blast_output_schema.keys()))

    # Choose the best alignments:
    ## option 1: prioritize coverage

    ## option 2: prioritize precision

    t_end = time.time()
    tsprint(f"{blast_output_fmt6}: Run time {t_end - t_start} seconds.")
    return "it worked"


def main():
    assert len(sys.argv) > 1
    blast_output_fmt6 = sys.argv[1]
    parse_blastn_outfmt6(blast_output_fmt6)

if __name__ == "__main__":
    main()
