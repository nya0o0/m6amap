import streamlit as st
import requests
from bs4 import BeautifulSoup

# ------------------------------
# 1. KEGG data fetch functions
# ------------------------------
@st.cache_data(show_spinner=False)
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
        else:
            st.warning(f"Gene symbol '{sym}' not found in KEGG.")
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

# ------------------------------
# Streamlit UI
# ------------------------------
st.set_page_config(page_title="KEGG Pathway Viewer", layout="wide")
st.title("🧬 KEGG Pathway Viewer from Gene Symbols")

# Session state setup
if "ko_map" not in st.session_state:
    st.session_state.ko_map = {}
if "common_pathways" not in st.session_state:
    st.session_state.common_pathways = []
if "analysis_complete" not in st.session_state:
    st.session_state.analysis_complete = False

# Input
organism = st.text_input("Organism code (e.g., `hsa` for human):", value="hsa")
gene_input = st.text_area("Enter gene symbols (comma or newline separated):", height=150)
submitted = st.button("Analyze")

if submitted:
    # Reset session state before new analysis
    st.session_state.ko_map = {}
    st.session_state.common_pathways = []
    st.session_state.analysis_complete = False

    with st.spinner("🔄 Processing..."):
        symbol_to_kegg, kegg_to_symbol = fetch_kegg_gene_list(organism)
        gene_symbols = [x.strip() for x in gene_input.replace(",", "\n").split("\n") if x.strip()]
        kegg_ids = symbol_to_keggid(gene_symbols, symbol_to_kegg)

        if not kegg_ids:
            st.error("❌ No valid KEGG gene IDs found.")
        else:
            st.success("✅ Mapped gene symbols to KEGG IDs:")
            st.json({sym: symbol_to_kegg.get(sym.upper()) for sym in gene_symbols if sym.upper() in symbol_to_kegg})

            ko_map = get_ko_ids(kegg_ids, organism)
            if not ko_map:
                st.error("❌ No KO (Orthology) IDs found.")
            else:
                st.session_state.ko_map = ko_map
                pathway_map = get_pathways_from_ko(ko_map)
                st.subheader("📊 Pathways associated with each KO:")
                st.json(pathway_map)

                common_pathways = get_common_pathways(pathway_map)
                all_unique_pathways = sorted({p for plist in pathway_map.values() for p in plist})

                if not common_pathways:
                    st.warning("⚠️ No common pathways found. Displaying all available pathways instead.")
                    st.session_state.common_pathways = all_unique_pathways
                else:
                    st.success(f"✅ Found {len(common_pathways)} common pathway(s).")
                    st.session_state.common_pathways = common_pathways

                st.session_state.analysis_complete = True

# ------------------------------
# Pathway selection & display
# ------------------------------
if st.session_state.analysis_complete:
    st.success(f":D Found {len(st.session_state.common_pathways)} common pathways.")
    selected_pathway = st.selectbox("Select a pathway to view:", st.session_state.common_pathways)

    if selected_pathway:
        st.subheader(f"🧪 Pathway: {selected_pathway}")
        img_url, viewer_url = get_kegg_pathway_image_url(selected_pathway, st.session_state.ko_map)

        if img_url:
            st.image(img_url, caption=f"Pathway {selected_pathway} with KO terms highlighted", use_container_width=True)
        else:
            st.error(":( Could not retrieve pathway image. It may be missing or invalid.")

        st.markdown(f"[🔗 Open in KEGG Viewer]({viewer_url})", unsafe_allow_html=True)