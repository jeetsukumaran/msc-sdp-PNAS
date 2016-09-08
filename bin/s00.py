#! /usr/bin/env python

import sys

def log(msg):
    sys.stderr.write("- {}\n".format(msg))

def open_output_file_for_csv_writer(filepath, append=False):
    if filepath is None or filepath == "-":
        out = sys.stdout
    elif sys.version_info >= (3,0,0):
        out = open(filepath,
                "a" if append else "w",
                newline='')
    else:
        out = open(filepath, "ab" if append else "wb")
    return out

