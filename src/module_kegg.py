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