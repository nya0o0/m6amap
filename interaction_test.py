import streamlit as st
import requests
import math
from pyvis.network import Network
import tempfile
import os

# ------------------------------
# STRING API: Fetch expanded subnetwork
# ------------------------------
def fetch_expanded_network(gene_list, species=9606, required_score=400, add_nodes=10):
    url = "https://string-db.org/api/json/network"
    params = {
        "identifiers": "%0d".join(gene_list),
        "species": species,
        "required_score": required_score,
        "add_nodes": add_nodes,
        "caller_identity": "string_network_app"
    }
    response = requests.get(url, params=params)
    if response.status_code == 200:
        return response.json()
    else:
        return []

# ------------------------------
# STRING API: PPI Enrichment
# ------------------------------
def get_ppi_enrichment(gene_list, species=9606):
    url = "https://string-db.org/api/json/enrichment"
    params = {
        "identifiers": "%0d".join(gene_list),
        "species": species,
    }
    response = requests.get(url, params=params)
    if response.status_code == 200:
        enrichment = response.json()
        return next((e for e in enrichment if e['category'] == 'PPI enrichment'), None)
    return None

# ------------------------------
# Generate STRING.org link
# ------------------------------
def generate_string_url(gene_list, species=9606):
    return f"https://string-db.org/cgi/network?species={species}&identifiers={'%0D'.join(gene_list)}"

# ------------------------------
# Build and Save Network HTML
# ------------------------------
def build_network(interactions, highlight_genes):
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
                    color="#90ee90" if gene in highlight_genes else "#add8e6"  # light green / blue
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

# ------------------------------
# Streamlit UI
# ------------------------------
st.set_page_config(page_title="STRING Network Viewer", layout="wide")
st.title("🧬 STRING Network Analyzer")

with st.sidebar:
    st.header("Input Options")
    gene_input = st.text_area("Enter gene symbols (one per line or comma-separated):", "TP53, EGFR, BRCA1, AKT1, PTEN")
    required_score = st.slider("Minimum STRING confidence score", 400, 900, 700, step=100)
    extra_nodes = st.slider("Additional interactors from STRING", 0, 50, 10, step=5)
    submitted = st.button("Run Analysis")

if submitted:
    # Parse gene list
    gene_list = [g.strip() for g in gene_input.replace(",", "\n").split("\n") if g.strip()]
    if len(gene_list) < 2:
        st.error("Please enter at least two gene symbols.")
    else:
        with st.spinner("🔄 Fetching STRING network..."):
            interactions = fetch_expanded_network(gene_list, required_score=required_score, add_nodes=extra_nodes)

        if not interactions:
            st.warning("No interactions returned by STRING. Try lowering the confidence score or checking gene names.")
        else:
            st.success(f"✅ Network fetched: {len(interactions)} interactions.")
            html_path = build_network(interactions, highlight_genes=set(gene_list))

            st.subheader("🌐 Interactive STRING Network")
            with open(html_path, "r", encoding="utf-8") as f:
                html = f.read()
                st.components.v1.html(html, height=600, scrolling=True)

            with st.expander("📊 PPI Enrichment"):
                enrichment = get_ppi_enrichment(gene_list)
                if enrichment:
                    st.success(f"✅ PPI Enrichment p-value: {enrichment['p_value']:.2e}")
                    st.write(f"**Observed edges**: {enrichment['number_of_edges']}")
                    st.write(f"**Expected edges**: {enrichment['expected_number_of_edges']:.2f}")
                    st.write(f"**Enrichment score** (fold): {enrichment['enrichment']:.2f}")
                else:
                    st.warning("No PPI enrichment result found.")

            st.markdown("---")
            st.markdown("🔗 [View on STRING.org](%s)" % generate_string_url(gene_list), unsafe_allow_html=True)

            # Clean up temp file
            os.unlink(html_path)