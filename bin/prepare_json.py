#!/usr/bin/env python3

import sys
import json


def read_one_line(path):
    with open(path) as f:
        content = f.read().strip()
    return content

def parse_two_lines_metrics(path, keys: int, values: int):
    h = open(path)
    lines = h.readlines()
    header = lines[keys].strip().split("\t")
    values = lines[values].strip().split("\t")
    h.close()
    return {k:v for k,v in zip(header, values)}


def main(dup_metrics, wgs_metrics, algn_metrics, contamination_metrics):
    dup = parse_two_lines_metrics(dup_metrics, 6, 7)
    wgs = parse_two_lines_metrics(wgs_metrics, 6, 7)
    contam = parse_two_lines_metrics(contamination_metrics, 0, 1)
    algn = parse_two_lines_metrics(algn_metrics, 6, 9)

    with open("summary.json", "w") as f:
        json.dump(dict(**dup, **wgs, **contam, **algn), f, indent=4)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
