#!/usr/bin/env python

import requests
from pyvis.network import Network
import math
import tempfile


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
