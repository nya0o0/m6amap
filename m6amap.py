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

# Set up Streamlit page
st.set_page_config(page_title="m6amap", layout="wide")
st.title("m6amap: m6A Modification Annotation and Pathway Mapping")

# Shared input: allow users to input custom gene list
st.markdown("### Custom Gene List (Optional)")
use_custom_genes = st.checkbox("Use custom input gene list")

custom_input = []
if use_custom_genes:
    gene_input_text = st.text_area("Enter gene symbols (comma or newline separated)", height=150)
    if gene_input_text:
        custom_input = [g.strip() for g in gene_input_text.replace(",", "\n").splitlines() if g.strip()]

# Step-by-step workflow using tabs
tabs = st.tabs(["README","1. Run m6anet", "2. Annotate Results", "3. Match with m6ADD", "4. STRING Network", "5. KEGG Pathways", "6. GO Enrichment"])

# Tab 0: README
with tabs[0]:
    st.header("README")
    try:
        with open("README.md", "r", encoding="utf-8") as readme_file:
            readme_content = readme_file.read()

        st.markdown(readme_content, unsafe_allow_html=True)
    except FileNotFoundError:
        st.error("README.md file not found. Please ensure it exists in the same directory.")

# Tab 1: Run m6anet
with tabs[1]:
    st.header("Step 1: Run m6anet")
    st.markdown("""
        Detect m6A modifications from eventalign files.
        Fill in the required parameters below to start the m6anet process.
    """)

    eventalign_path = st.text_input("Eventalign file path", value="tests/data/eventalign.txt", help=":)REQUIRED:) Path to the eventalign file (output from nanopolish).")
    output_dir = st.text_input("Output directory", value="m6anet_output", help=":)REQUIRED:) Directory to save the output files.")
    n_processes = st.number_input("Number of processes", 1, 32, 4, help= "Number of processes to run.")
    iterations = st.number_input("Number of iterations", 100, 10000, 1000, help="Number of times m6anet iterates through each potential m6a sites.")

    if st.button("Run m6anet"):
        try:
            run_m6anet(eventalign_path, output_dir, n_processes, iterations)
            st.success("✧｡٩(ˊᗜˋ )و✧*｡ m6anet run completed successfully.")
        except Exception as e:
            st.error(f"m6anet run failed: {e}")

# Tab 2: Annotate Results
with tabs[2]:
    st.header("Step 2: Annotate m6anet Results")
    st.markdown("""
        Annotate detected modification sites using a GTF file.
        Provide the paths to the m6anet output and GTF file below.
    """)

    input_file = st.text_input("m6anet output CSV", value="m6anet_output/data.site_proba.csv", help="Path to the m6anet output CSV file.")
    gtf_file = st.text_input("GTF file", value="tests/data/gencode.v47.annotation.gtf", help="Path to the GTF file for annotation.")
    output_prefix = st.text_input("Output prefix", value="annotated_output", help="Prefix for the annotated output file.")

    if st.button("Annotate Results"):
        try:
            process_files(input_file, gtf_file, output_prefix)
            st.success("✧｡٩(ˊᗜˋ )و✧*｡ Annotation completed successfully.")
        except Exception as e:
            st.error(f"Annotation failed: {e}")

# Tab 3: Match with m6ADD
with tabs[3]:
    st.header("Step 3: Match with m6ADD")
    st.markdown("""
        \Match annotated results with m6ADD data.
        Configure the matching parameters below.
    """)

    datatype = st.selectbox("m6ADD datatype", ["MeTDiff", "RADARDiff"], help="Choose the datatype for m6ADD.")
    window = st.slider("Matching window (nt)", 1, 50, 5, help="The size of the window for fuzzy matching (in nucleotides).")

    if st.button("Run Matching"):
        try:
            m6anet_df = pd.read_csv(f"{output_prefix}_{input_file.split('/')[-1]}")
            m6anet_df = m6anet_df.dropna(subset=["genome_pos", "chromosome"])
            m6anet_df["chromosome"] = m6anet_df["chromosome"].astype(str).str.lower().str.replace("chr", "")
            m6anet_df["genome_pos"] = m6anet_df["genome_pos"].astype(int)

            m6add_df = load_m6add_data(datatype)
            matched_df = fuzzy_merge(m6anet_df, m6add_df, window=window)

            if "m6add_Symbol_id" in matched_df.columns:
                matched_df = matched_df.drop(columns=["m6add_Symbol_id"])

            st.success(f"✧｡٩(ˊᗜˋ )و✧*｡ Matched {len(matched_df)} entries.")
            st.dataframe(matched_df.head(20))

            csv = matched_df.to_csv(index=False).encode("utf-8")
            st.download_button("( ._. )""↓↓↓↓ Download matched data", data=csv, file_name="matched_m6a.csv", mime="text/csv")
        except Exception as e:
            st.error(f"Matching failed: {e}")


# Tab 4: STRING Interaction Network
with tabs[4]:
    st.header("Step 4: STRING Interaction Network")
    st.markdown("""
        **Purpose**: Visualize interaction networks for genes.
        Use annotated gene symbols by default, or provide custom input.
    """)

    min_score = st.slider("STRING required score", 400, 900, 700, step=100, help="Minimum interaction confidence score.")
    extra_nodes = st.slider("Extra interactors", 0, 50, 10, step=5, help="Number of additional interactors to add to the network.")

    if st.button("Run STRING Network"):
        try:
            if use_custom_genes and custom_input:
                gene_list = custom_input[:100]
            else:
                df = pd.read_csv(f"{output_prefix}_{input_file.split('/')[-1]}")
                gene_list = df["gene_name"].dropna().unique().tolist()[:100]

            if len(gene_list) < 2:
                st.error("Please provide at least two valid gene symbols.")
            else:
                with st.spinner("( ◡̀_◡́)ᕤ Fetching STRING interactions..."):
                    interactions = fetch_expanded_network(gene_list, required_score=min_score, add_nodes=extra_nodes)

                if not interactions:
                    st.warning("No STRING interactions found.")
                else:
                    st.success(f"✧｡٩(ˊᗜˋ )و✧*｡ STRING Network fetched with {len(interactions)} interactions.")
                    html_path = build_string_network(interactions, highlight_genes=set(gene_list))

                    st.subheader(" STRING Interaction Network")
                    with open(html_path, "r", encoding="utf-8") as f:
                        html = f.read()
                        components.html(html, height=600, scrolling=True)

                    st.markdown("---")
                    st.markdown("🔗 [View full network on STRING.org](%s)" % generate_string_url(gene_list), unsafe_allow_html=True)
        except Exception as e:
            st.error(f":( STRING network failed: {e}")

# Tab 5: KEGG Pathways
with tabs[5]:
    st.header("Step 5: KEGG Pathway Viewer")
    st.markdown("""
        **Purpose**: Map genes to KEGG pathways.
        Use annotated gene symbols by default, or provide custom input.
    """)

    organism_code = st.text_input("KEGG Organism Code", value="hsa", help="e.g., 'hsa' for Homo sapiens.")
    if "kegg_results" not in st.session_state:
        st.session_state.kegg_results = None

    if st.button("Analyze Pathways"):
        try:
            if use_custom_genes and custom_input:
                gene_list = custom_input[:100]
            else:
                df = pd.read_csv(f"{output_prefix}_{input_file.split('/')[-1]}")
                gene_list = df["gene_name"].dropna().unique().tolist()[:100]

            st.info("( ◡̀_◡́)ᕤ Fetching KEGG gene list...")
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
            st.subheader(f"Pathway: {selected_pathway}")
            img_path, viewer_url = get_kegg_pathway_image_url(selected_pathway, results["ko_map"])

            if img_path:
                try:
                    st.image(img_path, caption=f"{selected_pathway} (highlighted)", use_container_width=True)
                except Exception as e:
                    st.error(f"Image failed to load: {e}")
            else:
                st.warning("Pathway image not available.")

            st.markdown(f"[🔗 Open in KEGG Viewer]({viewer_url})", unsafe_allow_html=True)
            

# Tab 6: GO Term Enrichment
with tabs[6]:
    st.header("Step 6: GO Enrichment Analysis")
    st.markdown("""
        **Purpose**: Perform Gene Ontology (GO) enrichment analysis.
        Use annotated gene symbols by default, or provide custom input.
    """)

    if "go_results" not in st.session_state:
        st.session_state.go_results = None

    if st.button("Show GO Term Summary"):
        try:
            if use_custom_genes and custom_input:
                gene_symbols = custom_input[:50]
            else:
                df = pd.read_csv(f"{output_prefix}_{input_file.split('/')[-1]}")
                gene_symbols = df["gene_name"].dropna().unique().tolist()[:50]

            st.info("( ◡̀_◡́)ᕤ Fetching GO annotations ...")
            go_dict = get_go_annotations(gene_symbols)
            go_summary = summarize_go_terms(go_dict)

            if go_summary.empty:
                st.warning("XP No GO terms found for the selected genes.")
            else:
                st.success(":D GO terms fetched successfully.")

                # Store data in session state to avoid recomputation
                st.session_state.go_results = {
                "go_dict": go_dict,
                "go_summary": go_summary
            }

        except Exception as e:
            st.error(f":( GO term annotation failed: {e}")

    # Display bar chart and table (only this part updates)
    if st.session_state.go_results:
        go_dict = st.session_state.go_results["go_dict"]
        go_summary = st.session_state.go_results["go_summary"]

        st.subheader("📊 GO Term Enrichment Summary")
        top_n = st.slider("Top N GO Terms", 5, 30, 10, key="go_slider")

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