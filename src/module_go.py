#!/usr/bin/env python

import pandas as pd
import requests



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