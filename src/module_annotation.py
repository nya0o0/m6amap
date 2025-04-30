#!/usr/bin/env python

import pandas as pd
import numpy as np
import gffutils
from collections import defaultdict
import warnings
import subprocess
import streamlit as st

def process_files(input_file, gtf_file, output_prefix):
    """
    Annotate m6anet results using a GTF file and display progress in Streamlit.

    Args:
        input_file (str): Path to the input CSV file containing modification sites.
        gtf_file (str): Path to the GTF file.
        output_prefix (str): Prefix for the output annotated file.
    """

    with st.status(" Starting annotation process...", expanded=True) as status:
        st.write(f"📂 Loading input site file: `{input_file}`...")
        modification_sites = pd.read_csv(input_file)
        modification_sites['transcript_id'] = modification_sites['transcript_id'].str.split('.').str[0]
        

        st.write(f"🧬 Preparing GTF database from `{gtf_file}`...")
        db_file = gtf_file + ".db"
        try: # Check whether the database already exists
            db = gffutils.FeatureDB(db_file, keep_order=True)
            st.write(":D GTF database loaded.")
        except:
            st.write("🛠️ Creating GTF database...")
            db = gffutils.create_db(gtf_file, dbfn=db_file, force=True, keep_order=True, disable_infer_genes=True, disable_infer_transcripts=True)
            db = gffutils.FeatureDB(db_file, keep_order=True)
            st.write(":D GTF database loaded.")


        # Extract transcript IDs from the modification sites
        st.write("🔍 Extracting transcript data from GTF database...")
        transcript_ids_to_keep = set(modification_sites["transcript_id"])
        
        # Filter features from the GTF database to include only relevant transcript IDs
        transcript_data = []
        filtered_features = []
        for feature in db.all_features():
            if "transcript_id" in feature.attributes:
                transcript_id = feature.attributes["transcript_id"][0].split(".")[0]
                if transcript_id in transcript_ids_to_keep:
                    filtered_features.append(feature)  # Keep only relevant features
                    if feature.featuretype == "transcript":
                        transcript_data.append({
                            "transcript_id": feature.id.split(".")[0],
                            "gene_name": feature.attributes.get("gene_name", [None])[0],
                            "transcript_type": feature.attributes.get("transcript_type", [None])[0],
                        })
        
        # Turn the transcrip data to DataFrame and drop dupilcated data
        transcript_data = pd.DataFrame(transcript_data).drop_duplicates("transcript_id")

        
        # Group filtered features by transcript_id
        st.write("🧬 Computing transcript features...")
        grouped_features = defaultdict(list)
        for feature in filtered_features:
            transcript_id = feature.attributes["transcript_id"][0].split(".")[0]  # Remove version
            grouped_features[transcript_id].append(feature)

        # Get the transcripts features
        tx_features, exons = calculate_transcript_features(grouped_features, transcript_data)
        
        # Apply the convert_to_genome_coordinates function to each row of the modification_sites
        st.write("📍 Mapping transcript coordinates to genome...")

        modification_sites["annotation"] = modification_sites.apply(
            lambda row: convert_to_genome_coordinates(row["transcript_id"], row["transcript_position"], exons, tx_features), axis=1
        )
        # Expand dictionary columns
        annotation_df = modification_sites["annotation"].apply(pd.Series)
        modification_sites = pd.concat([modification_sites, annotation_df], axis=1).drop(columns=["annotation"])

        # Merge with transcript biotype information
        annotated_modification_sites = modification_sites.merge(transcript_data, on="transcript_id", how="left")
        
        # Sort the annotated result according to the transcript ids
        annotated_modification_sites = annotated_modification_sites.sort_values(by='transcript_id')

        # Save output
        st.write("💾 Saving annotated results...")

        output_file = f"{output_prefix}_{input_file.split('/')[-1]}"
        annotated_modification_sites.to_csv(output_file, index=False)

        st.write(f"✅ Annotation completed! File saved as `{output_file}`:D")
        status.update(label="✅ Annotation successful", state="complete")

def calculate_transcript_features(grouped_features, transcript_data):
    '''
    Compute transcript features: transcript length, CDS length, and UTR lengths.

    Args:
        grouped_features (dict): Dictionary of transcript IDs mapped to their genomic features.
        transcript_data (DataFrame): DataFrame containing transcript data (transcript_ids, gene_names, transcript_types).

    Returns:
        DataFrame: Updated transcript metadata with computed feature lengths.
    '''

    transcript_features = []  # Store the results for all transcripts
    exons = defaultdict(list)  # Store exon features dictionary

    for tx_id, features in grouped_features.items():
        cds_length = 0
        utr5_length = 0
        utr3_length = 0
        tx_length = 0
        start_codon_pos = None
        stop_codon_pos = None

        # Separate features by type and store exon features in a dictionary
        cds_features = [f for f in features if f.featuretype == "CDS"]
        utr_features = [f for f in features if f.featuretype == "UTR"]
        exons[tx_id] = [f for f in features if f.featuretype == "exon"]

        # If no CDS features exist, set all lengths to 0， and compute transcript length by adding exons
        if not cds_features:
            cds_length = 0
            utr5_length = 0
            utr3_length = 0
            for exon in exons[tx_id]:
                tx_length += exon.end - exon.start + 1

        else:
            # Sort CDS features by start and end positions
            sorted_cds = sorted(cds_features, key=lambda x: (x.start, x.end))
            if utr_features[0].strand == "+":  # Positive strand
                start_codon_pos = sorted_cds[0].start  # First CDS start
                stop_codon_pos = sorted_cds[-1].end  # Last CDS end
            elif utr_features[0].strand == "-":  # Negative strand
                start_codon_pos = sorted_cds[-1].end  # Last CDS end
                stop_codon_pos = sorted_cds[0].start  # First CDS start

            # Calculate CDS length by summing the lengths of all CDS features
            for cds in sorted_cds:
                cds_length += cds.end - cds.start + 1

            # Handle UTRs and classify them as UTR5 or UTR3
            if utr_features:
                strand = utr_features[0].strand  # Assume all UTRs share the same strand
                for utr in utr_features:
                    utr_length = utr.end - utr.start + 1
                    if strand == "+":
                        if utr.end < start_codon_pos:  # 5' UTR is before start_codon
                            utr5_length += utr_length
                        elif utr.start > stop_codon_pos:  # 3' UTR is after stop_codon
                            utr3_length += utr_length
                    else:
                        if utr.start > start_codon_pos:  # 5' UTR is after start_codon on negative strand
                            utr5_length += utr_length
                        elif utr.end < stop_codon_pos:  # 3' UTR is before stop_codon on negative strand
                            utr3_length += utr_length
                
            tx_length = cds_length + utr3_length + utr5_length
    
        # Add the results for this transcript
        transcript_features.append({
            "transcript_id": tx_id,
            "tx_length": tx_length,
            "cds_length": cds_length,
            "utr5_length": utr5_length,
            "utr3_length": utr3_length
        })

    return pd.DataFrame(transcript_features).merge(transcript_data, on="transcript_id", how="left"), exons

# Function to convert transcript coordinates to genomic coordinates
def convert_to_genome_coordinates(tx_name, tx_pos, exonsdb, txfdb):
    """
    Convert transcript coordinates to genome coordinates and annotate regions.

    Args:
        tx_name (str): Transcript ID.
        tx_pos (int): Position in transcript coordinates.
        exonsdb (dict): Dictionary of transcript IDs mapped to exon features.
        txfdb (DataFrame): DataFrame containing transcript feature data.

    Returns:
        dict: Dictionary containing the genomic position, chromosome, exon junction distances, and region annotation.
    """
    
    if tx_name not in exonsdb:
        warnings.warn(f"Transcript {tx_name} not found in exons database.")
        return {
            "genome_pos": None,
            "chromosome": None,
            "dist_up_exon_junc": None,
            "dist_down_exon_junc": None,
            "region": None
        }
    
    tx_exons = sorted(exonsdb[tx_name], key=lambda x: x.start)  # Sort exons by start
    chrom = tx_exons[0].chrom
    strand = tx_exons[0].strand
    
    # Calculate cumulative exon lengths
    exon_lengths = [exon.end - exon.start + 1 for exon in tx_exons]
    cum_lengths = np.cumsum(exon_lengths)

    # Identify the exon containing the position
    idx = np.searchsorted(cum_lengths, tx_pos)
    
    if idx >= len(tx_exons):
        warnings.warn(f"Transcript position {tx_pos} is out of range for transcript {tx_name}")
        return {
            "genome_pos": None,
            "chromosome": None,
            "dist_up_exon_junc": None,
            "dist_down_exon_junc": None,
            "region": None
        }
    
    # Compute genomic position
    rel_pos = tx_pos - (cum_lengths[idx - 1] if idx > 0 else 0)
    exon = tx_exons[idx]
    if strand == "+":
        genome_pos = exon.start + rel_pos - 1
    else:
        genome_pos = exon.end - rel_pos + 1
    
    # Determine distance to exon junctions
    dist_upstream = rel_pos
    dist_downstream = exon.end - exon.start - rel_pos

    # Determine region (UTR, CDS, etc.)
    tx_feat = txfdb[txfdb["transcript_id"] == tx_name].iloc[0]
    if tx_feat["cds_length"] == 0:
        region = "ncRNA"
    elif tx_pos <= tx_feat["utr5_length"]:
        region = "UTR5"
    elif tx_pos <= tx_feat["utr5_length"] + tx_feat["cds_length"]:
        region = "CDS"
    elif tx_pos <= tx_feat["tx_length"]:
        region = "UTR3"
    else:
        region = None
    
    return {
        "genome_pos": genome_pos,
        "chromosome": chrom,
        "dist_up_exon_junc": dist_upstream,
        "dist_down_exon_junc": dist_downstream,
        "region": region
    }

def run_m6anet(eventalign_path, output_dir, n_processes=4, num_iterations=1000):

    try:
        import m6anet
    except ImportError:
        raise ImportError("m6anet is not installed. Please install it with `conda install m6anet` or `pip install m6anet`")
    
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