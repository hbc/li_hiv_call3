#!/usr/bin/env python
"""Prepare control library files: convert CSV/Excel input to FASTA.

Usage:
  prepare_controls.py <input_csv>
"""
import csv
import sys

def main(in_file):
    with open(in_file) as in_handle:
        reader = csv.reader(in_handle)
        cur_library = None
        cur_items = []
        for parts in reader:
            if "Start" in parts:
                if cur_items:
                    write_library(cur_library, cur_items)
                cur_library = parts[0]
                cur_items = []
            elif cur_library and parts[0]:
                name, percent, start, end, size, seq = parts
                cur_name = "%s_%s-%s_%s" % (name, start, end, percent)
                cur_items.append((cur_name, seq))
        if cur_items:
            write_library(cur_library, cur_items)

def write_library(libname, items):
    out_file = "%s-control.fa" % libname.replace(" ", "_")
    with open(out_file, "w") as out_handle:
        for name, seq in items:
            out_handle.write(">%s\n%s\n" % (name, seq))

if __name__ == "__main__":
    main(*sys.argv[1:])
