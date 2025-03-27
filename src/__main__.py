#!/usr/bin/env python

"""
Command line interface to m6alinker
"""

import argparse
from src.module import process_files, calculate_transcript_features, convert_to_genome_coordinates, rolling_progress


def parse_command_line():
    "parses args for the m6alinker funtion"

    # init parser and add arguments
    parser = argparse.ArgumentParser()

    # add input file args
    parser.add_argument(
        "-I", "--input_file",
        required=True,
        help="Path to the m6anet output file for annotation.")

    # add gtf file args
    parser.add_argument(
        "-G", "--gtf_file",
        required=True,
        help="Path to the GTF file for reference.")

    # add output prefix args
    parser.add_argument(
        "-O", "--output_prefix",
        help="Prefix for the output file. (Defalt: annotated_output)",
        default="annotated_output")

    # parse args
    args = parser.parse_args()

    '''
    # check that user only entered one action arg
    if sum([args.next, args.last, args.info]) > 1:
        raise SystemExit(
            "only one of 'next', 'last' or 'info' at a time.")
    '''
    return args


def main():
    "run main function on parsed args"

    # get arguments from command line as a dict-like object
    args = parse_command_line()

    # Pass arguments to process function
    process_files(args.input_file, args.gtf_file, args.output_prefix)


if __name__ == "__main__":
    main()
