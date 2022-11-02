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
    header = lines[keys].strip().replace('"', '').split("\t")
    values = lines[values].strip().replace('"', '').split("\t")
    h.close()
    return {k:v for k,v in zip(header, values)}


def main(contamination_metrics):
    contam = parse_two_lines_metrics(contamination_metrics, 0, 1)

    with open("summary.json", "w") as f:
        json.dump(contam, f, indent=4)

if __name__ == "__main__":
    main(sys.argv[1])
