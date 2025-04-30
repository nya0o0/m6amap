#!/usr/bin/env python

import pandas as pd
import numpy as np
import gffutils
from collections import defaultdict
import warnings
import glob
import subprocess
import streamlit as st
import requests
from pyvis.network import Network
import math
import tempfile
from bs4 import BeautifulSoup



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

def load_m6add_data(datatype):
    """
    Load and preprocess m6ADD dataset files.
    """
    files = glob.glob(f"tests/data/GSE*_{datatype}.txt")
    dfs = []
    for f in files:
        df = pd.read_csv(f, sep="\t")
        dfs.append(df)
    df = pd.concat(dfs)
    df["chr"] = df["chr"].astype(str).str.lower().str.replace("chr", "")
    df["genome_pos"] = ((df["start"] + df["end"]) / 2).astype(int)
    return df

def fuzzy_merge(m6anet_df, m6add_df, window=3):
    """
    Perform fuzzy matching between m6anet and m6add dataframes based on genomic position.
    """
    matches = []
    for _, row in m6anet_df.iterrows():
        chr_match = m6add_df["chr"] == str(row["chromosome"])
        pos_match = (m6add_df["genome_pos"] >= row["genome_pos"] - window) & \
                    (m6add_df["genome_pos"] <= row["genome_pos"] + window)
        match_df = m6add_df[chr_match & pos_match]
        for _, match_row in match_df.iterrows():
            combined = row.to_dict()
            combined.update({f"m6add_{col}": match_row[col] for col in m6add_df.columns})
            matches.append(combined)
    return pd.DataFrame(matches)

def fetch_kegg_gene_list(organism="hsa"):
    url = f"https://rest.kegg.jp/list/{organism}"
    response = requests.get(url)
    symbol_to_kegg = {}
    kegg_to_symbol = {}
    if response.status_code == 200:
        for line in response.text.strip().split("\n"):
            parts = line.split("\t")
            if len(parts) >= 4:
                kegg_id = parts[0].split(":")[1]
                gene_info = parts[3]
                gene_symbol = gene_info.split(",")[0].strip()
                symbol_to_kegg[gene_symbol.upper()] = kegg_id
                kegg_to_symbol[kegg_id] = gene_symbol.upper()
    return symbol_to_kegg, kegg_to_symbol

def symbol_to_keggid(symbols, symbol_to_kegg):
    ids = []
    for sym in symbols:
        sym = sym.upper()
        if sym in symbol_to_kegg:
            ids.append(symbol_to_kegg[sym])
    return ids

def get_ko_ids(kegg_ids, organism="hsa"):
    ko_map = {}
    for kid in kegg_ids:
        url = f"http://rest.kegg.jp/link/ko/{organism}:{kid}"
        res = requests.get(url)
        if res.ok and res.text:
            lines = res.text.strip().split("\n")
            if lines:
                parts = lines[0].split("\t")
                if len(parts) == 2 and parts[1].startswith("ko:"):
                    ko_map[kid] = parts[1]
    return ko_map

def get_pathways_from_ko(ko_map):
    pathway_map = {}
    for kid, ko_id in ko_map.items():
        url = f"http://rest.kegg.jp/link/pathway/{ko_id}"
        res = requests.get(url)
        if res.ok and res.text:
            pathways = []
            for line in res.text.strip().split("\n"):
                parts = line.strip().split("\t")
                if len(parts) == 2 and parts[1].startswith("path:map"):
                    pathway_id = parts[1].split(":")[1]
                    pathways.append(pathway_id)
            pathway_map[kid] = pathways
    return pathway_map

def get_common_pathways(pathway_map):
    all_pathways = list(pathway_map.values())
    if not all_pathways:
        return []
    common = set(all_pathways[0])
    for p in all_pathways[1:]:
        common &= set(p)
    return list(common)

def get_kegg_pathway_image_url(pathway_id, ko_map):
    if not ko_map:
        return None, None

    ko_list = "+".join([ko.split(":")[1] for ko in ko_map.values()])
    url = f"https://www.kegg.jp/kegg-bin/show_pathway?{pathway_id}/{ko_list}%09red"
    res = requests.get(url)
    if res.ok:
        soup = BeautifulSoup(res.text, "html.parser")
        img_tag = soup.find("img")
        if img_tag and "src" in img_tag.attrs:
            img_src = img_tag["src"]
            full_img_url = "https://www.kegg.jp" + img_src
            return full_img_url, url
    return None, url

def get_go_annotations(gene_symbols):
    """
    Get GO term metadata for a list of gene symbols by using UniProt + QuickGO APIs.

    Args:
        gene_symbols (list): List of gene symbols (e.g., TP53, ACTB).

    Returns:
        dict: gene_symbol -> list of GO term dicts: {
            "id": GO ID,
            "name": term name,
            "aspect": BP/MF/CC,
            "definition": string
        }
    """
    import time

    gene_go = {}
    headers = {"Accept": "application/json"}

    def fetch_go_term_detail(go_id):
        """Fetch GO term name, aspect, and definition from QuickGO term API."""
        url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{go_id}"
        try:
            resp = requests.get(url, headers=headers)
            resp.raise_for_status()
            data = resp.json()["results"][0]
            return {
                "id": go_id,
                "name": data.get("name", ""),
                "aspect": data.get("aspect", ""),
                "definition": data.get("definition", {}).get("text", "")
            }
        except Exception as e:
            print(f"⚠️ Failed to fetch GO term detail for {go_id}: {e}")
            return {
                "id": go_id,
                "name": "",
                "aspect": "",
                "definition": ""
            }

    for gene in gene_symbols:
        try:
            # Step 1: Get UniProt ID
            uniprot_url = f"https://rest.uniprot.org/uniprotkb/search?query=gene_exact:{gene}+AND+organism_id:9606&fields=accession&format=json"
            response = requests.get(uniprot_url, headers=headers)
            response.raise_for_status()
            results = response.json().get("results", [])

            if not results:
                print(f"⚠️ No UniProt ID found for gene: {gene}")
                gene_go[gene] = []
                continue

            uniprot_id = results[0]["primaryAccession"]

            # Step 2: Get GO annotations for UniProt ID
            quickgo_url = f"https://www.ebi.ac.uk/QuickGO/services/annotation/search"
            params = {
                "geneProductId": f"UniProtKB:{uniprot_id}",
                "limit": 100
            }

            go_response = requests.get(quickgo_url, headers=headers, params=params)
            go_response.raise_for_status()
            go_results = go_response.json().get("results", [])

            go_ids = list(set([entry["goId"] for entry in go_results if "goId" in entry]))
            go_terms = []

            for go_id in go_ids:
                go_term = fetch_go_term_detail(go_id)
                go_terms.append(go_term)
                time.sleep(0.1)  # Avoid overloading API

            gene_go[gene] = go_terms
        except Exception as e:
            print(f":( Error fetching GO terms for {gene}: {e}")
            gene_go[gene] = []

    return gene_go


def summarize_go_terms(go_dict):
    """
    Count GO term frequencies from a gene-to-GO mapping.

    Args:
        go_dict (dict): gene -> list of GO term dicts.

    Returns:
        DataFrame: GO term frequencies with name and aspect.
    """
    from collections import Counter

    all_terms = []
    for terms in go_dict.values():
        for t in terms:
            term_label = f"{t['name']} ({t['aspect']})" if t['name'] else t['id']
            all_terms.append(term_label)

    counter = Counter(all_terms)
    df = pd.DataFrame(counter.items(), columns=["GO_term", "Count"])
    return df.sort_values(by="Count", ascending=False)

def fetch_expanded_network(gene_list, species=9606, required_score=400, add_nodes=10):
    url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": "%0d".join(gene_list),
        "species": species,
        "required_score": required_score,
        "add_nodes": add_nodes,
        "caller_identity": "m6amap_app"
    }
    response = requests.get(url, params=params)
    if response.status_code == 200:
        return response.json()
    else:
        return []

def get_ppi_enrichment(gene_list, species=9606):
    url = "https://string-db.org/api/json/enrichment"
    params = {
        "identifiers": "%0d".join(gene_list),
        "species": species
    }
    response = requests.get(url, params=params)
    if response.status_code == 200:
        enrichment = response.json()
        return next((e for e in enrichment if e['category'] == 'PPI enrichment'), None)
    return None

def generate_string_url(gene_list, species=9606):
    return f"https://string-db.org/cgi/network?species={species}&identifiers={'%0D'.join(gene_list)}"

def build_string_network(interactions, highlight_genes):
    net = Network(height="600px", width="100%", bgcolor="#111", font_color="white")
    net.repulsion(node_distance=150)

    seen_nodes = set()

    for entry in interactions:
        source = entry["preferredName_A"]
        target = entry["preferredName_B"]
        score = entry["score"]

        for gene in [source, target]:
            if gene not in seen_nodes:
                net.add_node(
                    gene,
                    label=gene,
                    color="#90ee90" if gene in highlight_genes else "#add8e6"
                )
                seen_nodes.add(gene)

        thickness = max(1, min(8, int(math.log(score + 1, 1.5)) - 5))

        net.add_edge(
            source,
            target,
            value=thickness,
            color='white',
            title=f"Confidence score: {score:.0f}"
        )

    tmp_file = tempfile.NamedTemporaryFile(delete=False, suffix=".html")
    net.write_html(tmp_file.name)
    return tmp_file.name
