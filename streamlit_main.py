import streamlit as st
import pandas as pd
from src.module import process_files, load_m6add_data, fuzzy_merge
from src.__main__ import run_m6anet

st.set_page_config(page_title="m6amap", layout="wide")
st.title("m6A Annotation & Matching Tool")

st.sidebar.header("1. Run m6anet")
eventalign_path = st.sidebar.text_input("Eventalign file", value="tests/data/eventalign.txt")
output_dir = st.sidebar.text_input("Output directory", value="m6anet_output")
n_processes = st.sidebar.number_input("Num processes", 1, 32, 4)
iterations = st.sidebar.number_input("Num iterations", 100, 10000, 1000)

if st.sidebar.button("Run m6anet"):
    try:
        run_m6anet(eventalign_path, output_dir, n_processes, iterations)
        st.success("m6anet run completed.")
    except Exception as e:
        st.error(f"Failed: {e}")

st.sidebar.header("2. Annotate m6anet Results")
input_file = st.sidebar.text_input("m6anet output CSV", value="m6anet_output/output.csv")
gtf_file = st.sidebar.text_input("GTF file", value="tests/data/gencode.v47.annotation.gtf")
output_prefix = st.sidebar.text_input("Output prefix", value="annotated_output")

if st.sidebar.button("Annotate with GTF"):
    try:
        process_files(input_file, gtf_file, output_prefix)
        st.success("Annotation completed.")
    except Exception as e:
        st.error(f"Annotation failed: {e}")

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

        st.success(f"✅ Matched {len(matched_df)} entries.")
        st.dataframe(matched_df.head(20))

        csv = matched_df.to_csv(index=False).encode("utf-8")
        st.download_button("📥 Download matched data", data=csv, file_name="matched_m6a.csv", mime="text/csv")
    except Exception as e:
        st.error(f"❌ Matching failed: {e}")