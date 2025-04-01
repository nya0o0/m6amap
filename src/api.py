import requests

def get_kegg_gene_id(uniprot_id, organism_code="hsa"):
    """
    Convert UniProt ID to KEGG gene ID using KEGG 'conv' operation.
    Example: get_kegg_gene_id("P31749") => "hsa:5290"
    """
    url = f"https://rest.kegg.jp/conv/{organism_code}/uniprot:{uniprot_id}"
    response = requests.get(url)
    if response.status_code != 200 or not response.text:
        return None
    return response.text.strip().split('\t')[1]

def link_kegg_pathways(kegg_gene_id):
    """
    Get pathways linked to a KEGG gene ID using the 'link' operation.
    Example: link_kegg_pathways("hsa:5290")
    """
    url = f"https://rest.kegg.jp/link/pathway/{kegg_gene_id}"
    response = requests.get(url)
    if response.status_code != 200 or not response.text:
        return []
    lines = response.text.strip().split('\n')
    pathway_ids = [line.split('\t')[1].split(":")[1] for line in lines]
    return list(set(pathway_ids))

def get_pathway_names(pathway_ids):
    """
    Fetch readable names for a list of KEGG pathway IDs using the 'list' operation.
    Example: get_pathway_names(["hsa04151"]) => {"hsa04151": "PI3K-Akt signaling pathway"}
    """
    if not pathway_ids:
        return {}
    joined_ids = "+".join(pathway_ids)
    url = f"https://rest.kegg.jp/list/{joined_ids}"
    response = requests.get(url)
    if response.status_code != 200 or not response.text:
        return {}
    lines = response.text.strip().split('\n')
    return {line.split('\t')[0]: line.split('\t')[1] for line in lines}

#def query_go_terms(uniprot_id):

#def query_go_annotations(uniprot_id):


if __name__ == "__main__":
    # Step 1: Known UniProt ID (e.g., AKT1)
    uniprot_id = "P31749"  # AKT1 (well-studied, should be in KEGG)

    print(f"🔎 Querying KEGG gene ID for UniProt ID: {uniprot_id}")
    kegg_gene_id = get_kegg_gene_id(uniprot_id)

    if not kegg_gene_id:
        print("❌ Failed to retrieve KEGG gene ID.")
    else:
        print(f"✅ KEGG Gene ID: {kegg_gene_id}")

        # Step 2: Get linked pathways
        print("\n🔗 Fetching pathways linked to the gene...")
        pathways = link_kegg_pathways(kegg_gene_id)

        if not pathways:
            print("❌ No pathways found.")
        else:
            print(f"✅ Pathway IDs: {pathways}")

            # Step 3: Get pathway names
            print("\n🧬 Resolving pathway names...")
            pathway_names = get_pathway_names(pathways)

            for pid, name in pathway_names.items():
                print(f"🧠 {pid}: {name}")
