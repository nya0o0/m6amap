#!/usr/bin/env python

"""
Command line interface to m6amap
"""

import argparse
import subprocess
import os
import sys
from src.module import process_files, calculate_transcript_features, convert_to_genome_coordinates, rolling_progress


def run_m6anet(eventalign_path, output_dir, n_processes=4, num_iterations=1000):

    try:
        import m6anet
    except ImportError:
        raise ImportError("m6anet is not installed. Please install it with 'conda install m6anet'")
    
    # Step 1: Dataprep
    dataprep_command = [
        "m6anet", "dataprep",
        "--eventalign", eventalign_path,
        "--out_dir", output_dir,
        "--n_processes", str(n_processes)
    ]

    print("Running m6anet dataprep...")
    result = subprocess.run(dataprep_command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError("m6anet dataprep failed.")
    else:
        print("m6anet dataprep finished.")

    # Step 2: Inference
    inference_command = [
        "m6anet", "inference",
        "--input_dir", output_dir,
        "--out_dir", output_dir,
        "--n_processes", str(n_processes),
        "--num_iterations", str(num_iterations)
    ]

    print("Running m6anet inference...")
    result = subprocess.run(inference_command, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError("m6anet inference failed.")
    else:
        print("m6anet inference finished.")
        print(f"Output files save to {output_dir}/.")


def parse_command_line():
    "parses args for the m6alinker funtion"

    # init parser and add arguments
    parser = argparse.ArgumentParser(description="Run m6anet and annotate m6anet output according to GTF file.")

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
    
    # add run m6anet args
    parser.add_argument(
        "--run_m6anet", action="store_true",
        help="Run m6anet inference pipeline."
    )

    parser.add_argument(
        "--eventalign", help="Path to eventalign.txt"
    )

    parser.add_argument(
        "--out_dir", help="Output directory for m6anet", default="m6anet_output"
    )

    parser.add_argument(
        "--n_proc", type=int, default=4, help="Number of processes for m6anet"
    )

    parser.add_argument(
        "--iterations", type=int, default=1000, help="Number of iterations for m6anet"
    )

    # parse args
    args = parser.parse_args()

    return args


def main():
    "run main function on parsed args"

    # get arguments from command line as a dict-like object
    args = parse_command_line()

    # Pass arguments to process function
    if args.run_m6anet:
        run_m6anet(args.eventalign, args.out_dir, args.n_proc, args.iterations)

    if args.input_file and args.gtf_file:
        process_files(args.input_file, args.gtf_file, args.output_prefix)


if __name__ == "__main__":
    main()
