import streamlit as st
import pandas as pd
from pyvis.network import Network
import tempfile
import streamlit.components.v1 as components
import altair as alt
from src.module import process_files, load_m6add_data, fuzzy_merge, run_m6anet, fetch_gene_annotations, build_interaction_graph, get_go_annotations, summarize_go_terms

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

st.sidebar.header("4. Gene Interaction Graph")

if st.sidebar.button("Build Enhanced Network Graph"):
    try:
        # Load annotated file
        df = pd.read_csv(f"{output_prefix}_{input_file.split('/')[-1]}")
        original_genes = df["gene_name"].dropna().unique().tolist()[:50]  # Limit for performance

        st.write("🔍 Fetching gene annotations from KEGG...")
        gene_info = fetch_gene_annotations(original_genes)

        st.write("🔗 Building gene-gene interaction graph...")
        G, interaction_sources = build_interaction_graph(original_genes)

        net = Network(height="600px", width="100%", bgcolor="#ffffff", font_color="black", notebook=False)
        net.barnes_hut()

        node_colors = {
            "query_gene": "skyblue",
            "interactor_gene": "orange"
        }

        # Add gene nodes with annotations
        for node in G.nodes():
            node_type = G.nodes[node].get("type", "unknown")
            color = node_colors.get(node_type, "gray")

            # Get annotation info for tooltip
            info = gene_info.get(node, {})
            pathways = info.get("pathways", [])
            diseases = info.get("diseases", [])

            tooltip = f"<b>{node}</b><br><br>"
            tooltip += f"<b>Pathways:</b><br>{'<br>'.join(pathways) if pathways else 'None'}<br><br>"
            tooltip += f"<b>Diseases:</b><br>{'<br>'.join(diseases) if diseases else 'None'}"

            net.add_node(node, label=node, color=color, title=tooltip)

        # Add edges with descriptions
        for source, target in G.edges():
            title = interaction_sources.get((source, target), "Interaction")
            net.add_edge(source, target, title=title)

        # Save and render the network graph
        with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as tmp_file:
            net.save_graph(tmp_file.name)
            components.html(open(tmp_file.name, "r", encoding="utf-8").read(), height=600)

        st.success(":D Interactive gene interaction network generated.")

    except Exception as e:
        st.error(f":() Graph generation failed: {e}")

st.sidebar.header("5. GO Term Enrichment")

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