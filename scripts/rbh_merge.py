#!/usr/bin/env python

"""
This script is for merging two RBH files into one new file.

Input: two RBH files
 - The path to two RBH files
 - The output path where the the new RBH file will be saved
"""

import argparse
import os
import sys

# import rbh_tools
scriptpath = os.path.dirname(os.path.realpath(__file__))
source_path = os.path.join(scriptpath, "../source")
sys.path.insert(1, source_path)
import rbh_tools

def parse_args():
    """
    - Two RBH files, allow multiple entries, allow expansions
    - One output file
    - A flag to specify whether the rbh files are database rbhs
    """
    parser = argparse.ArgumentParser(description='Merge two RBH files into one new file.')
    parser.add_argument("-r", '--rbh_files', nargs='+', help='The path to two or more RBH files')
    parser.add_argument("-o", '--output_file', help='The output path where the the new RBH file will be saved')
    parser.add_argument("-d", '--database', action='store_true', help='Flag to specify whether the rbh files are database rbhs')
    args = parser.parse_args()
    # make sure that the input files exist
    for file in args.rbh_files:
        if not os.path.exists(file):
            parser.error(f"The file {file} does not exist.")
    # make sure that the output file ends with .rbh
    if not args.output_file.endswith('.rbh'):
        parser.error("The output file must end with .rbh")
    # Make sure that there are only two rbh files as input, otherwise raise an error
    if len(args.rbh_files) != 2:
        parser.error("Only two RBH files are allowed as input")
    return parser.parse_args()

def main():
    args = parse_args()
    filepath_1 = args.rbh_files[0]
    filepath_2 = args.rbh_files[1]
    output_file = args.output_file

    # if it is a database, use the special function for database rbh files
    if args.database:
        df = rbh_tools.combine_rbh_db(filepath_1, filepath_2)
        df.to_csv(output_file, sep='\t', index=False)
    else:
        # merge the two rbh files
        df = rbh_tools.combine_rbh(filepath_1, filepath_2)
        # output the file
        df.to_csv(output_file, sep='\t', index=False)

if __name__ == '__main__':
    main()
