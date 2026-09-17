# 🧬 VAXINTAIC — Vaccine Intelligence for Antigenic Insight and Construction
![Python](https://img.shields.io/badge/python-3.12-blue)
![Status](https://img.shields.io/badge/status-active-success)
![Platform](https://img.shields.io/badge/platform-bioinformatics-green)
VAXINTAIC is an AI-driven computational bioinformatics platform designed to analyze viral genomic data, identify conserved antigenic regions, and evaluate candidate vaccine constructs using automated pipelines and machine learning.

---

# 🚀 Overview

VAXINTAIC provides an end-to-end system that transforms raw viral sequence data into ranked vaccine candidate insights through:

* Sequence acquisition and preprocessing
* Multiple sequence alignment and diversity analysis
* Epitope prediction and ranking
* Synthetic construct simulation (mRNA-level)
* Machine learning-based efficacy prediction
* Interactive dashboard visualization

---

# 🧠 System Architecture

```text
External Sequence Data
        ↓
[01] Sequence Fetching
        ↓
[02] Alignment & Processing
        ↓
[03] Epitope Prediction
        ↓
[04] Epitope Ranking
        ↓
[05] Construct Design
        ↓
[06] Structural Modeling
        ↓
[07] AI Efficacy Prediction
        ↓
[08] Reports & Storage
        ↓
[Dashboard Visualization]
```

---

# 📂 Project Structure

```text
vaxintaic/
│
├── pipeline/                 # Core computational pipeline
│   ├── 01_fetch_sequences.py
│   ├── 02_align_sequences.py
│   ├── 03_epitope_prediction.py
│   ├── 04_epitope_ranking.py
│   ├── 05_mrna_design.py
│   ├── 06_nanoparticle_model.py
│   ├── 07_efficacy_predictor.py
│   ├── 09_diversity_analysis.py
│   ├── 10_immunogenicity_score.py
│   ├── 11_codon_adaptation.py
│   └── 13_generate_advanced_report.py
│
├── dashboard/               # Streamlit dashboard
│   ├── app.py
│   ├── part1_overview.py
│   ├── part2_epitope.py
│   ├── part3_diversity.py
│   └── part4_reports.py
│
├── data/                    # Generated outputs
│   ├── sequences/
│   ├── alignments/
│   ├── epitopes/
│   ├── candidates/
│   ├── mrna/
│   ├── models/
│   ├── analysis/
│   └── reports/
│
├── config/                  # Configuration files
├── docs/                    # Images/screenshots
│
├── run_vaxintaic_full.py    # Master pipeline controller
├── requirements.txt
├── README.md
└── .gitignore
```

---

# ⚙️ Installation Guide

## 1. Clone Repository

```bash
git clone https://github.com/YOUR_USERNAME/vaxintaic.git
cd vaxintaic
```

## 2. Create Virtual Environment

```bash
python3 -m venv venv
source venv/bin/activate
```

## 3. Install Dependencies

```bash
pip install -r requirements.txt
```

---

# ▶️ Running the Pipeline

```bash
python run_vaxintaic_full.py
```

This executes all pipeline stages sequentially:

* Sequence fetching
* Alignment
* Epitope prediction
* Ranking
* Construct design
* ML prediction
* Report generation

---

# 🌐 Running the Dashboard

```bash
streamlit run dashboard/app_master.py
```

Open in browser:

```text
http://localhost:8501
```

---

# 📊 Dashboard Features

* 🧭 Pipeline overview and metrics
* 📈 Construct performance analytics
* 🧫 Epitope visualization
* 🧬 Diversity and alignment insights
* 📄 Report and log viewer
* 🧠 JSON construct intelligence explorer

---

# 🧬 Pipeline Inspection Guide

## Inspect Epitope Prediction

```bash
nano pipeline/03_epitope_prediction.py
```

Key logic:

* Sliding window peptide scanning
* Immunogenic scoring

---

## Inspect mRNA Design

```bash
nano pipeline/05_mrna_design.py
```

Key logic:

* Construct assembly
* Codon optimization (CAI, GC%)

---

## Inspect AI Model

```bash
nano pipeline/07_efficacy_predictor.py
```

Features used:

* GC content
* Codon adaptation index
* Stability index
* Immunogenicity

---

# 🗄️ Database (MySQL)

VAXINTAIC uses MySQL for structured storage.

### Main Tables

#### pipeline_metrics_v3

Tracks pipeline runs:

```sql
timestamp_utc
num_sequences
construct_count
duration_sec
```

#### construct_intelligence_v3

Stores construct-level intelligence:

```sql
construct_id
prrsv_type
timestamp_utc
data (JSON)
```

---

# 📊 Data Flow

```text
Pipeline Scripts
      ↓
CSV / FASTA Outputs
      ↓
MySQL Database
      ↓
Streamlit Dashboard
```

---

# 🌍 Deployment

## Local Deployment

```bash
streamlit run dashboard/app_master.py --server.port 8090
```

---

## Server Deployment (Example)

```bash
http://YOUR_SERVER_IP:8090
```

---

## Live System (Example)

```text
https://prrsv-vaccine.aiconceptlimited.com.ng/
```

---

# 🧠 Technologies Used

* Python
* Pandas / NumPy
* BioPython
* Scikit-learn
* Streamlit
* Plotly
* MySQL

---

# ⚠️ Disclaimer

This platform is intended for computational research, simulation, and educational purposes only.

---

# 👨‍🔬 Author

**Abubakar**
VAXINTAIC — AI for Vaccine Intelligence

