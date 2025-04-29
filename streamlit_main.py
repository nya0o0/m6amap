import streamlit as st
import pandas as pd
from pyvis.network import Network
import tempfile
import math
import streamlit.components.v1 as components
import altair as alt
from src.module import process_files, load_m6add_data, fuzzy_merge, run_m6anet,  get_go_annotations, summarize_go_terms, fetch_expanded_network, get_ppi_enrichment, generate_string_url, build_string_network, fetch_kegg_gene_list, symbol_to_keggid, get_ko_ids, get_pathways_from_ko, get_common_pathways, get_kegg_pathway_image_url

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

st.sidebar.header("4. STRING Interaction Network")
min_score = st.sidebar.slider("STRING required score", 400, 900, 700, step=100)
extra_nodes = st.sidebar.slider("Extra interactors", 0, 50, 10, step=5)
run_string = st.sidebar.button("Run STRING Network")

if run_string:
    try:
        df = pd.read_csv(f"{output_prefix}_{input_file.split('/')[-1]}")
        gene_list = df["gene_name"].dropna().unique().tolist()[:100]  # limit for performance

        if len(gene_list) < 2:
            st.error("Please provide at least two genes from annotation.")
        else:
            with st.spinner("🔄 Fetching STRING interactions..."):
                interactions = fetch_expanded_network(gene_list, required_score=min_score, add_nodes=extra_nodes)

            if not interactions:
                st.warning("No STRING interactions found. Try lowering the threshold or checking gene names.")
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

if run_kegg_pathway:
    try:
        df = pd.read_csv(f"{output_prefix}_{input_file.split('/')[-1]}")
        gene_list = df["gene_name"].dropna().unique().tolist()[:100]

        st.info("🔍 Fetching KEGG gene list...")
        symbol_to_kegg, kegg_to_symbol = fetch_kegg_gene_list(organism_code)

        kegg_ids = symbol_to_keggid(gene_list, symbol_to_kegg)
        if not kegg_ids:
            st.error(":( No valid KEGG gene IDs found.")
        else:
            st.success(":D Mapped gene symbols to KEGG IDs.")
            st.json({sym: symbol_to_kegg.get(sym.upper()) for sym in gene_list if sym.upper() in symbol_to_kegg})

            ko_map = get_ko_ids(kegg_ids, organism=organism_code)
            if not ko_map:
                st.error(":( No KO (Orthology IDs) found.")
            else:
                pathway_map = get_pathways_from_ko(ko_map)
                st.subheader("📊 Pathways associated with KO terms:")
                st.json(pathway_map)

                common_pathways = get_common_pathways(pathway_map)
                all_pathways = sorted({p for plist in pathway_map.values() for p in plist})
                pathways_to_display = common_pathways if common_pathways else all_pathways

                if not pathways_to_display:
                    st.warning("XP No pathways found.")
                else:
                    st.success(f"✅ Found {len(pathways_to_display)} pathway(s).")
                    selected_pathway = st.selectbox("Select a pathway to view:", pathways_to_display)

                    if selected_pathway:
                        st.subheader(f"🧪 Pathway: {selected_pathway}")
                        img_url, viewer_url = get_kegg_pathway_image_url(selected_pathway, ko_map)

                        if img_url:
                            st.image(img_url, caption=f"{selected_pathway} (highlighted)", use_container_width=True)
                        else:
                            st.error(":( Could not retrieve image.")

                        st.markdown(f"[🔗 Open in KEGG Viewer]({viewer_url})", unsafe_allow_html=True)
    except Exception as e:
        st.error(f"KEGG pathway analysis failed: {e}")


st.sidebar.header("6. GO Term Enrichment")

if st.sidebar.button("Show GO Term Summary"):
    try:
        df = pd.read_csv(f"{output_prefix}_{input_file.split('/')[-1]}")
        gene_symbols = df["gene_name"].dropna().unique().tolist()[:50]  # Limit for performance

        st.info("🔍 Fetching GO annotations ...")
        go_dict = get_go_annotations(gene_symbols)
        go_summary = summarize_go_terms(go_dict)

        if go_summary.empty:
            st.warning("XP No GO terms found for the selected genes.")
        else:
            st.success(":D GO terms fetched successfully.")

            # Show bar chart of top GO terms
            top_n = st.slider("Top N GO Terms", 5, 30, 10)
            chart = alt.Chart(go_summary.head(top_n)).mark_bar().encode(
                x=alt.X("Count:Q"),
                y=alt.Y("GO_term:N", sort='-x'),
                tooltip=["GO_term", "Count"]
            ).properties(
                height=400,
                title="Top GO Terms Across Annotated Genes"
            )

            st.altair_chart(chart, use_container_width=True)

            # Show full GO term table
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

            # Optional: download button
            csv = df_terms.to_csv(index=False).encode("utf-8")
            st.download_button("📥 Download GO Annotations", data=csv, file_name="go_annotations.csv", mime="text/csv")

    except Exception as e:
        st.error(f":() GO term annotation failed: {e}")



