# **m6amap**  

related mini-project: https://github.com/nya0o0/m6alinker

## **1. Introduction**  

**What task will the project accomplish？**

This project aims to integrate [**m6Anet**](https://m6anet.readthedocs.io/en/latest/) and [**M6ADD**](http://m6add.edbc.org/) into a unified pipeline that takes **nanopore sequencing data** as input and outputs **predicted RNA modification sites (m6A)**, their **associated diseases, genes, and pathways**, with interactive visualizations for better interpretability.  

**Why is this useful?**  

- **Biomedical Insights**: This project links m6A sites to diseases, genes, and pathways, helping researchers better understand a more comperhansive map of RNA modifications in disease progression.  
- **Interactive Data Exploration**: Users will explore results dynamically using interactive graphs, making it easier to identify patterns and correlations.  

---

## **2. User Input and Data Requirements**  

**User Input**  
The user will provide:  
- **Nanopore direct RNA sequencing data** (BAM format)  
- **Reference transcriptome** (FA format)  
- **ONT reads file** (FASTA format)
- **Refernece anotated genome data** (GTF format)

**Data Sources**  
- **Nanopore sequencing data** (provided by the user)  
- **GENCODE refernece anotated genome data** from gencode release.([Example](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.annotation.gtf.gz))
- **M6ADD** database for linking modifications to diseases, genes, and pathways ([M6ADD Database](https://m6add.org))  
- **Gene Ontology (GO) and KEGG Pathway Databases** 
  - **Gene Ontology API (via QuickGO)**  [https://www.ebi.ac.uk/QuickGO/](https://www.ebi.ac.uk/QuickGO/)  
  - **KEGG API**  [https://www.kegg.jp/kegg/rest/keggapi.html](https://www.kegg.jp/kegg/rest/keggapi.html)  

---

## **3. How Users Will Interact with the Program**   

1. **Command-Line Interface (CLI)** for running the programme
2. **Web-based Interactive Dashboard** for exploring results  

**Example Usage:**  

```bash
m6amap --reads reads.fa --bam alignments.bam --genome genome.fa -G gencode.v47.annotation.gtf -O output_perfix
```

**Expected Outputs:**  
- **Annotated m6Anet Output (CSV)**:  
  - RNA modification sites with confidence scores  and annotation
- **Interactive Visualization**:  
  - **Network graphs** linking m6A sites to diseases and pathways. 

---

## **4. Expected Output Format & Interactive Visualization**  

### **Text Output Example (CSV)**  

**Annotated m6Anet results (CSV):**

| transcript_id | gene_id | transcript_position | n_reads | probability_modified | kmer | mod_ratio | sample_id | group_id | genome_pos | chromosome | dist_up_exon_junc | dist_down_exon_junc | region | gene_name | transcript_type |Associated Disease (from M6ADD) | Pathways |
|--------------|--------|---------------------|---------|----------------------|------|-----------|-----------|---------|------------|------------|-------------------|--------------------|--------|-----------|-----------------|---------| -------|
| ENST00000263741 | ENSG00000078808 | 1466 | 42 | 0.9288 | GGACT | 0.7619 | s2 | caud | 1228941 | chr1 | 6 | 472 | UTR3 | SDF4 | protein_coding | Neurodevelopmental Disorders | mTOR Signaling
| ENST00000354700 | ENSG00000131584 | 3222 | 78 | 0.9559 | AGACA | 0.7949 | s2 | caud | 1300627 | chr1 | 66 | 117 | UTR3 | ACAP3 | protein_coding | Cancer |	RNA Degradation |
| ENST00000379198 | ENSG00000176022 | 2612 | 31 | 0.9527 | GGACT | 0.6774 | s2 | caud | 1234848 | chr1 | 2612 | 192 | UTR3 | B3GALT6 | protein_coding | Connective Tissue Disorders | Glycosaminoglycan Biosynthesis
| ENST00000416718 | ENSG00000198744 | 341 | 801 | 0.9272 | GGACA | 0.6667 | s1 | caud | 634716 | chr1 | 341 | 205 | ncRNA | MTCO3P12 | unprocessed_pseudogene | N/A | N/A |

### **Interactive Visualization**  

**1. Network Graph (m6A Sites → Genes → Diseases → Pathways)**  
- **Nodes**: genes
- **Edges**: Relationships between them  
- **Click**: View detailed annotations about the related pathways and diseases

**2. Possible method:**
  - **Interactive Visualization Using Python (Dash + Plotly)**
 A python-based dashboard.
  
---
## **5. What other tools currently exist to do this task, or something similar?**

There are several tools that can predict m6a modification sites and also databases that annotate modification sites, but this project try to link this two kinds of tools.

---

## **6. In Development** 

### Installation
```bash
git clone https://github.com/nya0o0/project
cd  src
pip install -e .
```