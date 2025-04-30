import streamlit as st
import pandas as pd
from pyvis.network import Network
import tempfile
import math
import streamlit.components.v1 as components
import altair as alt
from src.module_annotation import process_files, run_m6anet
from src.module_match import load_m6add_data, fuzzy_merge
from src.module_interaction import fetch_expanded_network, get_ppi_enrichment, generate_string_url, build_string_network
from src.module_kegg import fetch_kegg_gene_list, symbol_to_keggid, get_ko_ids, get_pathways_from_ko, get_common_pathways, get_kegg_pathway_image_url
from src.module_go import get_go_annotations, summarize_go_terms

st.set_page_config(page_title="m6amap", layout="wide")
st.title("m6amap :D")

st.sidebar.header("1. Run m6anet")
eventalign_path = st.sidebar.text_input("Eventalign file", value="tests/data/eventalign.txt")
output_dir = st.sidebar.text_input("Output directory", value="m6anet_output")
n_processes = st.sidebar.number_input("Num processes", 1, 32, 4)
iterations = st.sidebar.number_input("Num iterations", 100, 10000, 1000)

if st.sidebar.button("Run m6anet"):
    try:
        run_m6anet(eventalign_path, output_dir, n_processes, iterations)
        st.success(":D m6anet run completed.")
    except Exception as e:
        st.error(f":( m6anet run failed: {e}")

st.sidebar.header("2. Annotate m6anet Results")
input_file = st.sidebar.text_input("m6anet output CSV", value="m6anet_output/data.site_proba.csv")
gtf_file = st.sidebar.text_input("GTF file", value="tests/data/gencode.v47.annotation.gtf")
output_prefix = st.sidebar.text_input("Output prefix", value="annotated_output")

if st.sidebar.button("Annotate with GTF"):
    try:
        process_files(input_file, gtf_file, output_prefix)
        st.success(":D Annotation completed.")
    except Exception as e:
        st.error(f":( Annotation failed: {e}")

st.sidebar.header("3. Match with m6ADD")
datatype = st.sidebar.selectbox("m6ADD datatype", ["MeTDiff", "RADARDiff"])
window = st.sidebar.slider("Matching window (nt)", 1, 50, 5)

if st.sidebar.button("Run Matching"):
    try:
        m6anet_df = pd.read_csv(f"{output_prefix}_{input_file.split('/')[-1]}")
        m6anet_df = m6anet_df.dropna(subset=["genome_pos", "chromosome"])
        m6anet_df["chromosome"] = m6anet_df["chromosome"].astype(str).str.lower().str.replace("chr", "")
        m6anet_df["genome_pos"] = m6anet_df["genome_pos"].astype(int)

        m6add_df = load_m6add_data(datatype)
        matched_df = fuzzy_merge(m6anet_df, m6add_df, window=window)

        # Drop the undesired column
        if "m6add_Symbol_id" in matched_df.columns:
            matched_df = matched_df.drop(columns=["m6add_Symbol_id"])

        st.success(f":D Matched {len(matched_df)} entries.")
        st.dataframe(matched_df.head(20))

        csv = matched_df.to_csv(index=False).encode("utf-8")
        st.download_button("📥 Download matched data", data=csv, file_name="matched_m6a.csv", mime="text/csv")
    except Exception as e:
        st.error(f":( Matching failed: {e}")

# Shared input: allow users to input custom gene list
st.sidebar.markdown("### Custom Gene List (Optional)")
use_custom_genes = st.sidebar.checkbox("Use custom input gene list")

custom_input = []
if use_custom_genes:
    gene_input_text = st.sidebar.text_area("Enter gene symbols (comma or newline separated)", height=150)
    if gene_input_text:
        custom_input = [g.strip() for g in gene_input_text.replace(",", "\n").splitlines() if g.strip()]

st.sidebar.header("4. STRING Interaction Network")
min_score = st.sidebar.slider("STRING required score", 400, 900, 700, step=100)
extra_nodes = st.sidebar.slider("Extra interactors", 0, 50, 10, step=5)
run_string = st.sidebar.button("Run STRING Network")

if run_string:
    try:
        if use_custom_genes and custom_input:
            gene_list = custom_input[:100]
        else:
            df = pd.read_csv(f"{output_prefix}_{input_file.split('/')[-1]}")
            gene_list = df["gene_name"].dropna().unique().tolist()[:100]

        if len(gene_list) < 2:
            st.error("Please provide at least two valid gene symbols.")
        else:
            with st.spinner("🔄 Fetching STRING interactions..."):
                interactions = fetch_expanded_network(gene_list, required_score=min_score, add_nodes=extra_nodes)

            if not interactions:
                st.warning("No STRING interactions found.")
            else:
                st.success(f":D STRING Network fetched with {len(interactions)} interactions.")
                html_path = build_string_network(interactions, highlight_genes=set(gene_list))

                st.subheader("🌐 STRING Interaction Network")
                with open(html_path, "r", encoding="utf-8") as f:
                    html = f.read()
                    components.html(html, height=600, scrolling=True)

                with st.expander("📊 PPI Enrichment"):
                    enrichment = get_ppi_enrichment(gene_list)
                    if enrichment:
                        st.success(f"PPI Enrichment p-value: {enrichment['p_value']:.2e}")
                        st.write(f"**Observed edges**: {enrichment['number_of_edges']}")
                        st.write(f"**Expected edges**: {enrichment['expected_number_of_edges']:.2f}")
                        st.write(f"**Enrichment score**: {enrichment['enrichment']:.2f}")
                    else:
                        st.warning("No enrichment result found.")

                st.markdown("---")
                st.markdown("🔗 [View full network on STRING.org](%s)" % generate_string_url(gene_list), unsafe_allow_html=True)
    except Exception as e:
        st.error(f":( STRING network failed: {e}")

st.sidebar.header("5. KEGG Pathway Viewer")
organism_code = st.sidebar.text_input("KEGG Organism Code", value="hsa")
run_kegg_pathway = st.sidebar.button("Analyze Pathways")

if "kegg_results" not in st.session_state:
    st.session_state.kegg_results = None

# Run the analysis and store results in session state
if run_kegg_pathway:
    try:
        if use_custom_genes and custom_input:
            gene_list = custom_input[:100]
        else:
            df = pd.read_csv(f"{output_prefix}_{input_file.split('/')[-1]}")
            gene_list = df["gene_name"].dropna().unique().tolist()[:100]

        st.info("🔍 Fetching KEGG gene list...")
        symbol_to_kegg, kegg_to_symbol = fetch_kegg_gene_list(organism_code)

        kegg_ids = symbol_to_keggid(gene_list, symbol_to_kegg)
        if not kegg_ids:
            st.error(":( No valid KEGG gene IDs found.")
        else:
            ko_map = get_ko_ids(kegg_ids, organism=organism_code)
            pathway_map = get_pathways_from_ko(ko_map)
            common_pathways = get_common_pathways(pathway_map)
            all_pathways = sorted({p for plist in pathway_map.values() for p in plist})
            pathways_to_display = common_pathways if common_pathways else all_pathways

            # Store results in session state
            st.session_state.kegg_results = {
                "gene_list": gene_list,
                "symbol_to_kegg": symbol_to_kegg,
                "kegg_ids": kegg_ids,
                "ko_map": ko_map,
                "pathway_map": pathway_map,
                "pathways_to_display": pathways_to_display
            }

            st.success(":D KEGG analysis complete. You can now select a pathway below.")

    except Exception as e:
        st.error(f"KEGG pathway analysis failed: {e}")

# If results already exist, show pathway selector and image
if st.session_state.kegg_results:
    results = st.session_state.kegg_results

    st.subheader("📊 Pathways associated with KO terms:")
    pathway_rows = []
    for gene_id, pathways in results["pathway_map"].items():
        for pid in pathways:
            pathway_rows.append({"KEGG Gene ID": gene_id, "Pathway ID": pid})
    st.dataframe(pd.DataFrame(pathway_rows))

    selected_pathway = st.selectbox("Select a pathway to view:", results["pathways_to_display"])

    if selected_pathway:
        st.subheader(f"🧪 Pathway: {selected_pathway}")
        img_url, viewer_url = get_kegg_pathway_image_url(selected_pathway, results["ko_map"])

        if img_url:
            try:
                st.image(img_url, caption=f"{selected_pathway} (highlighted)", use_container_width=True)
            except Exception as e:
                st.error(f"Image failed to load: {e}")
        else:
            st.warning("⚠️ Pathway image not available.")

        st.markdown(f"[🔗 Open in KEGG Viewer]({viewer_url})", unsafe_allow_html=True)


st.sidebar.header("6. GO Term Enrichment")

if "go_results" not in st.session_state:
    st.session_state.go_results = None

# Button to run GO analysis
if st.sidebar.button("Show GO Term Summary"):
    try:
        if use_custom_genes and custom_input:
            gene_symbols = custom_input[:50]
        else:
            df = pd.read_csv(f"{output_prefix}_{input_file.split('/')[-1]}")
            gene_symbols = df["gene_name"].dropna().unique().tolist()[:50]

        st.info("🔍 Fetching GO annotations ...")
        go_dict = get_go_annotations(gene_symbols)
        go_summary = summarize_go_terms(go_dict)

        if go_summary.empty:
            st.warning("XP No GO terms found for the selected genes.")
        else:
            st.success(":D GO terms fetched successfully.")

            # Store results
            st.session_state.go_results = {
                "go_dict": go_dict,
                "go_summary": go_summary
            }

    except Exception as e:
        st.error(f":() GO term annotation failed: {e}")

# If results exist, show chart and table
if st.session_state.go_results:
    go_dict = st.session_state.go_results["go_dict"]
    go_summary = st.session_state.go_results["go_summary"]

    st.subheader("📊 GO Term Enrichment Summary")

    # Slider only affects chart, not the whole computation
    top_n = st.slider("Top N GO Terms", 5, 30, 10)

    # Show chart
    chart = alt.Chart(go_summary.head(top_n)).mark_bar().encode(
        x=alt.X("Count:Q"),
        y=alt.Y("GO_term:N", sort='-x'),
        tooltip=["GO_term", "Count"]
    ).properties(
        height=400,
        title="Top GO Terms Across Annotated Genes"
    )
    st.altair_chart(chart, use_container_width=True)

    # Show full GO annotation table
    st.markdown("### 📋 Detailed GO Annotation Table")
    df_terms = []
    for gene, terms in go_dict.items():
        for t in terms:
            df_terms.append({
                "Gene": gene,
                "GO ID": t["id"],
                "Name": t["name"],
                "Aspect": t["aspect"],
                "Definition": t["definition"]
            })

    df_terms = pd.DataFrame(df_terms)
    st.dataframe(df_terms)

    # Download
    csv = df_terms.to_csv(index=False).encode("utf-8")
    st.download_button("📥 Download GO Annotations", data=csv, file_name="go_annotations.csv", mime="text/csv")

