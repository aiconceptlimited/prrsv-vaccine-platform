import os
#!/usr/bin/env python3
# =========================================================
# 🧬 VAXINTAIC — PRRSV Vaccine Intelligence Dashboard
# Author: Abubakar | 2026
# Disease Focus: Porcine Reproductive and Respiratory Syndrome Virus (PRRSV)
# Genotypes: Type I (European) & Type II (North American)
# =========================================================

import streamlit as st
import pandas as pd
import mysql.connector
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
# VAXINTAIC database configuration.
# Credentials are supplied through the environment and must never be committed.
DB_HOST = os.environ.get("VAXINTAIC_DB_HOST", "localhost")
DB_PORT = int(os.environ.get("VAXINTAIC_DB_PORT", "3306"))
DB_USER = os.environ.get("VAXINTAIC_DB_USER", "vaxuser")
DB_PASSWORD = os.environ.get("VAXINTAIC_DB_PASSWORD")
DB_NAME = os.environ.get("VAXINTAIC_DB_NAME", "vaxintaic")

if not DB_PASSWORD:
    raise RuntimeError(
        "VAXINTAIC_DB_PASSWORD is not set; refusing database connection."
    )

DB_URL = (
    f"mysql+mysqlconnector://{DB_USER}:{DB_PASSWORD}"
    f"@{DB_HOST}:{DB_PORT}/{DB_NAME}"
)

# =========================================================
# PAGE CONFIG
# =========================================================
st.set_page_config(
    page_title="VAXINTAIC v3 — PRRSV Vaccine Intelligence Platform",
    layout="wide",
    page_icon="🧬",
    initial_sidebar_state="expanded"
)

# =========================================================
# NEON DARK THEME - PROFESSIONAL DASHBOARD STYLE
# =========================================================
st.markdown("""
<style>
/* Main theme colors */
:root {
    --bg-primary: #0a0e17;
    --bg-secondary: #0f1525;
    --bg-container: #1a2238;
    --bg-header: #0d1117;
    --bg-footer: #0a0e17;
    --neon-blue: #00d4ff;
    --neon-purple: #8a2be2;
    --neon-green: #00ff88;
    --neon-cyan: #00f7ff;
    --text-primary: #ffffff;
    --text-secondary: #b0b7c3;
    --border-color: #2a3552;
    --card-shadow: 0 8px 32px rgba(0, 212, 255, 0.1);
    --neon-glow: 0 0 20px rgba(0, 212, 255, 0.3);
}

/* Apply dark theme to all elements */
.stApp {
    background: linear-gradient(135deg, var(--bg-primary) 0%, var(--bg-secondary) 100%);
    color: var(--text-primary);
}

/* Headers */
h1, h2, h3, h4, h5, h6 {
    color: var(--text-primary) !important;
    font-weight: 600 !important;
    background: linear-gradient(90deg, var(--neon-blue), var(--neon-purple));
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
    border-bottom: 1px solid var(--border-color);
    padding-bottom: 10px;
    margin-bottom: 20px !important;
}

/* Professional containers */
.main .block-container {
    padding-top: 2rem;
    padding-bottom: 2rem;
}

/* Custom container styling */
.st-emotion-cache-1r6slb0, .st-emotion-cache-1y4p8pa, .element-container {
    background-color: var(--bg-container) !important;
    border: 1px solid var(--border-color) !important;
    border-radius: 12px !important;
    padding: 20px !important;
    margin-bottom: 20px !important;
    box-shadow: var(--card-shadow) !important;
}

/* Metrics styling */
[data-testid="metric-container"] {
    background-color: var(--bg-container) !important;
    border: 1px solid var(--border-color) !important;
    border-radius: 10px !important;
    padding: 15px !important;
    box-shadow: var(--card-shadow) !important;
}

[data-testid="metric-container"]:hover {
    border-color: var(--neon-cyan) !important;
    box-shadow: var(--neon-glow) !important;
}

[data-testid="stMetricLabel"], [data-testid="stMetricValue"] {
    color: var(--text-primary) !important;
}

[data-testid="stMetricDelta"] {
    color: var(--neon-green) !important;
}

/* Sidebar styling */
section[data-testid="stSidebar"] {
    background-color: var(--bg-header) !important;
    border-right: 1px solid var(--border-color) !important;
}

section[data-testid="stSidebar"] h1, 
section[data-testid="stSidebar"] h2, 
section[data-testid="stSidebar"] h3,
section[data-testid="stSidebar"] .stRadio > label {
    color: var(--text-primary) !important;
}

/* Dataframes */
.stDataFrame {
    background-color: var(--bg-container) !important;
    border: 1px solid var(--border-color) !important;
    border-radius: 8px !important;
}

/* Radio buttons and selection */
.stRadio > div {
    background-color: var(--bg-container) !important;
    border-radius: 8px !important;
    padding: 10px !important;
}

.stRadio > div > label {
    color: var(--text-primary) !important;
}

/* Text elements */
p, li, td, th, div, span {
    color: var(--text-primary) !important;
}

/* Warning/Info boxes */
.stAlert, .stWarning, .stInfo, .stSuccess, .stError {
    background-color: var(--bg-container) !important;
    border: 1px solid var(--border-color) !important;
    color: var(--text-primary) !important;
    border-radius: 8px !important;
}

/* Expander styling */
.streamlit-expanderHeader {
    background-color: var(--bg-container) !important;
    color: var(--text-primary) !important;
    border: 1px solid var(--border-color) !important;
    border-radius: 8px !important;
}

.streamlit-expanderContent {
    background-color: var(--bg-secondary) !important;
    color: var(--text-primary) !important;
}

/* Badge styling */
.badge {
    padding: 6px 12px;
    border-radius: 20px;
    font-size: 0.8rem;
    font-weight: 500;
    margin: 2px 5px 2px 0;
    display: inline-block;
    background: linear-gradient(135deg, var(--bg-secondary), var(--bg-container));
    border: 1px solid var(--border-color);
    color: var(--text-primary) !important;
}

.badge-success {
    background: linear-gradient(135deg, rgba(0, 255, 136, 0.1), rgba(0, 255, 136, 0.2));
    border-color: var(--neon-green);
    color: var(--neon-green) !important;
}

.badge-warning {
    background: linear-gradient(135deg, rgba(255, 193, 7, 0.1), rgba(255, 193, 7, 0.2));
    border-color: #ffc107;
    color: #ffc107 !important;
}

.badge-info {
    background: linear-gradient(135deg, rgba(0, 212, 255, 0.1), rgba(0, 212, 255, 0.2));
    border-color: var(--neon-blue);
    color: var(--neon-blue) !important;
}

/* Footer styling */
.footer-container {
    background-color: var(--bg-footer) !important;
    border-top: 2px solid var(--border-color) !important;
    padding: 20px !important;
    margin-top: 40px !important;
    border-radius: 12px !important;
}

.footer {
    text-align: center;
    font-size: 0.9rem;
    color: var(--text-secondary) !important;
}

.footer strong {
    color: var(--neon-cyan) !important;
    background: linear-gradient(90deg, var(--neon-blue), var(--neon-cyan));
    -webkit-background-clip: text;
    -webkit-text-fill-color: transparent;
}

/* Plotly chart background fix */
.js-plotly-plot, .plotly, .modebar {
    background-color: transparent !important;
}

/* Custom scrollbar */
::-webkit-scrollbar {
    width: 8px;
    height: 8px;
}

::-webkit-scrollbar-track {
    background: var(--bg-secondary);
}

::-webkit-scrollbar-thumb {
    background: var(--border-color);
    border-radius: 4px;
}

::-webkit-scrollbar-thumb:hover {
    background: var(--neon-blue);
}
</style>
""", unsafe_allow_html=True)

# =========================================================
# DATABASE LOADER
# =========================================================
def load_mysql():
    conn = mysql.connector.connect(
        host="localhost",
        user="vaxuser",
        password=DB_PASSWORD,
        database="vaxintaic"
    )

    profiles = pd.read_sql("SELECT * FROM v_construct_profiles_valid", conn)
    logical = pd.read_sql("SELECT * FROM v_logical_constructs", conn)
    ranking = pd.read_sql("SELECT * FROM v_construct_ranking", conn)
    explanations = pd.read_sql("SELECT * FROM v_construct_explanations", conn)
    exclusions = pd.read_sql("SELECT * FROM v_construct_exclusion_reasons", conn)
    pipeline = pd.read_sql(
        "SELECT * FROM pipeline_metrics_v3 ORDER BY timestamp_utc DESC",
        conn
    )

    conn.close()
    return profiles, logical, ranking, explanations, exclusions, pipeline

profiles_df, logical_df, ranking_df, explanations_df, exclusions_df, pipeline_df = load_mysql()

# =========================================================
# SIDEBAR
# =========================================================
st.sidebar.title("🧬 VAXINTAIC v3")
st.sidebar.markdown("""
<div style='color: #00d4ff; font-size: 0.9rem; margin-bottom: 20px;'>
<strong>PRRSV Vaccine Intelligence Platform</strong><br/>
<em style='color: #b0b7c3;'>Porcine Reproductive and Respiratory Syndrome Virus</em>
</div>
""", unsafe_allow_html=True)

section = st.sidebar.radio(
    "Navigate",
    [
        "📊 Overview",
        "🧠 Construct Intelligence",
        "🧫 PRRSV Type Comparison",
        "🚚 Delivery & Computational Scoring",
        "🏆 Construct Ranking",
        "📈 Pipeline Trends",
        "ℹ️ Scientific Scope",
    ]
)

# =========================================================
# OVERVIEW
# =========================================================
if section == "📊 Overview":
    st.title("🧬 VAXINTAIC v3 — PRRSV Vaccine Intelligence Platform")

    st.markdown("""
<div style='background-color: #1a2238; padding: 20px; border-radius: 12px; border: 1px solid #2a3552; margin-bottom: 20px;'>
<strong style='color: #00d4ff;'>Disease Focus:</strong> Porcine Reproductive and Respiratory Syndrome Virus (PRRSV)<br/>
<strong style='color: #00d4ff;'>Genotypes Covered:</strong> Type I (European) and Type II (North American)
</div>
""", unsafe_allow_html=True)

    latest = pipeline_df.iloc[0]

    c1, c2, c3, c4 = st.columns(4)
    with c1:
        st.metric("Logical Vaccine Constructs", len(logical_df))
    with c2:
        st.metric("Valid Construct Profiles", len(profiles_df))
    with c3:
        st.metric("Pipeline Runtime (s)", round(latest["duration_sec"], 2))
    with c4:
        st.metric("Last Pipeline Run", str(latest["timestamp_utc"])[:19])

    # Apply dark theme to plot
    fig = px.bar(
        logical_df,
        x="construct_id",
        y="prrsv_coverage",
        title="PRRSV Genotype Coverage per Vaccine Construct",
        labels={
            "construct_id": "Vaccine Construct",
            "prrsv_coverage": "Number of PRRSV Types Covered"
        }
    )
    
    # Update plot layout for dark theme
    fig.update_layout(
        plot_bgcolor='rgba(26, 34, 56, 0.8)',
        paper_bgcolor='rgba(10, 14, 23, 0.9)',
        font_color='#ffffff',
        title_font_color='#00d4ff',
        xaxis_title_font_color='#b0b7c3',
        yaxis_title_font_color='#b0b7c3',
        xaxis=dict(
            gridcolor='rgba(42, 53, 82, 0.5)',
            tickfont=dict(color='#ffffff')
        ),
        yaxis=dict(
            gridcolor='rgba(42, 53, 82, 0.5)',
            tickfont=dict(color='#ffffff')
        )
    )
    
    fig.update_traces(marker_color='#00d4ff', marker_line_color='#0088cc')
    
    st.plotly_chart(fig, use_container_width=True)

# =========================================================
# CONSTRUCT INTELLIGENCE
# =========================================================
elif section == "🧠 Construct Intelligence":
    st.header("🧠 Computational Vaccine Construct Intelligence")

    if profiles_df.empty:
        st.warning("No constructs currently meet strict scientific validation criteria.")
    else:
        st.dataframe(profiles_df, use_container_width=True)

        st.subheader("🧠 Why These Constructs Were Selected")
        for _, r in explanations_df.iterrows():
            st.markdown(
                f"<div style='background-color: #1a2238; padding: 15px; border-radius: 8px; border: 1px solid #2a3552; margin-bottom: 10px;'>"
                f"<strong style='color: #00d4ff;'>{r.construct_id} ({r.prrsv_type})</strong> — {r.explanation}"
                f"</div>",
                unsafe_allow_html=True
            )

        def badge_row(row):
            badges = []
            if row.delivery_index >= 0.65:
                badges.append("<span class='badge badge-success'>🟢 Delivery-stable</span>")
            if 40 <= row.gc_percent <= 60:
                badges.append("<span class='badge badge-success'>🟢 Balanced GC</span>")
            return " ".join(badges)

        profiles_df["Robustness"] = profiles_df.apply(badge_row, axis=1)

        st.subheader("🛡 Robustness Indicators")
        st.markdown(
            profiles_df[["construct_id", "prrsv_type", "Robustness"]].to_html(
                escape=False, 
                index=False,
                classes='dataframe',
                border=0
            ),
            unsafe_allow_html=True
        )

        if not exclusions_df.empty:
            with st.expander("⚠️ Excluded Constructs & Reasons"):
                st.dataframe(exclusions_df, use_container_width=True)

# =========================================================
# PRRSV TYPE COMPARISON
# =========================================================
elif section == "🧫 PRRSV Type Comparison":
    st.header("🧫 PRRSV Type I vs Type II Performance")

    fig = px.box(
        profiles_df,
        x="prrsv_type",
        y="predicted_model_score",
        points="all",
        title="In-Silico Model Score by PRRSV Type"
    )
    
    # Update plot layout for dark theme
    fig.update_layout(
        plot_bgcolor='rgba(26, 34, 56, 0.8)',
        paper_bgcolor='rgba(10, 14, 23, 0.9)',
        font_color='#ffffff',
        title_font_color='#00d4ff',
        xaxis_title_font_color='#b0b7c3',
        yaxis_title_font_color='#b0b7c3',
        xaxis=dict(
            gridcolor='rgba(42, 53, 82, 0.5)',
            tickfont=dict(color='#ffffff')
        ),
        yaxis=dict(
            gridcolor='rgba(42, 53, 82, 0.5)',
            tickfont=dict(color='#ffffff')
        )
    )
    
    st.plotly_chart(fig, use_container_width=True)

# =========================================================
# DELIVERY & EFFICACY
# =========================================================
elif section == "🚚 Delivery & Computational Scoring":
    st.header("🚚 Nanoparticle Delivery vs Computational Model Score")

    fig = px.scatter(
        profiles_df,
        x="delivery_index",
        y="predicted_model_score",
        size="mean_immunogenicity",
        color="prrsv_type",
        hover_data=["construct_id", "num_epitopes", "gc_percent"],
        title="Delivery–Computational-Score Relationship (PRRSV)"
    )
    
    # Update plot layout for dark theme
    fig.update_layout(
        plot_bgcolor='rgba(26, 34, 56, 0.8)',
        paper_bgcolor='rgba(10, 14, 23, 0.9)',
        font_color='#ffffff',
        title_font_color='#00d4ff',
        xaxis_title_font_color='#b0b7c3',
        yaxis_title_font_color='#b0b7c3',
        xaxis=dict(
            gridcolor='rgba(42, 53, 82, 0.5)',
            tickfont=dict(color='#ffffff')
        ),
        yaxis=dict(
            gridcolor='rgba(42, 53, 82, 0.5)',
            tickfont=dict(color='#ffffff')
        ),
        legend=dict(
            bgcolor='rgba(26, 34, 56, 0.8)',
            bordercolor='#2a3552',
            borderwidth=1
        )
    )
    
    st.plotly_chart(fig, use_container_width=True)

# =========================================================
# CONSTRUCT RANKING
# =========================================================
elif section == "🏆 Construct Ranking":
    st.header("🏆 PRRSV Vaccine Construct Ranking")

    st.dataframe(ranking_df, use_container_width=True)

    fig = px.bar(
        ranking_df,
        x="construct_id",
        y="composite_score",
        title="Overall Construct Ranking (Composite Score)"
    )
    
    # Update plot layout for dark theme
    fig.update_layout(
        plot_bgcolor='rgba(26, 34, 56, 0.8)',
        paper_bgcolor='rgba(10, 14, 23, 0.9)',
        font_color='#ffffff',
        title_font_color='#00d4ff',
        xaxis_title_font_color='#b0b7c3',
        yaxis_title_font_color='#b0b7c3',
        xaxis=dict(
            gridcolor='rgba(42, 53, 82, 0.5)',
            tickfont=dict(color='#ffffff')
        ),
        yaxis=dict(
            gridcolor='rgba(42, 53, 82, 0.5)',
            tickfont=dict(color='#ffffff')
        )
    )
    
    fig.update_traces(marker_color='#00ff88', marker_line_color='#00cc66')
    
    st.plotly_chart(fig, use_container_width=True)

# =========================================================
# PIPELINE TRENDS
# =========================================================
elif section == "📈 Pipeline Trends":
    st.header("📈 Pipeline Runtime History")

    fig = px.line(
        pipeline_df,
        x="timestamp_utc",
        y="duration_sec",
        markers=True,
        title="Pipeline Execution Time Over Runs"
    )
    
    # Update plot layout for dark theme
    fig.update_layout(
        plot_bgcolor='rgba(26, 34, 56, 0.8)',
        paper_bgcolor='rgba(10, 14, 23, 0.9)',
        font_color='#ffffff',
        title_font_color='#00d4ff',
        xaxis_title_font_color='#b0b7c3',
        yaxis_title_font_color='#b0b7c3',
        xaxis=dict(
            gridcolor='rgba(42, 53, 82, 0.5)',
            tickfont=dict(color='#ffffff')
        ),
        yaxis=dict(
            gridcolor='rgba(42, 53, 82, 0.5)',
            tickfont=dict(color='#ffffff')
        )
    )
    
    fig.update_traces(line_color='#00d4ff', line_width=3)
    
    st.plotly_chart(fig, use_container_width=True)

# =========================================================
# SCIENTIFIC SCOPE
# =========================================================
elif section == "ℹ️ Scientific Scope":
    st.header("ℹ️ Scientific Scope & Methodology")

    st.markdown("""
<div style='background-color: #1a2238; padding: 25px; border-radius: 12px; border: 1px solid #2a3552;'>
<h3 style='color: #00d4ff; margin-top: 0;'>Disease Target</h3>
<strong style='color: #00ff88;'>Porcine Reproductive and Respiratory Syndrome Virus (PRRSV)</strong>

<h3 style='color: #00d4ff; margin-top: 20px;'>Genotypes Covered</h3>
<ul style='color: #ffffff;'>
<li>Type I (European)</li>
<li>Type II (North American)</li>
</ul>

<h3 style='color: #00d4ff; margin-top: 20px;'>Platform Capabilities</h3>
<ul style='color: #ffffff;'>
<li>Epitope discovery & ranking</li>
<li>mRNA construct design</li>
<li>Codon optimization</li>
<li>Nanoparticle delivery modeling</li>
<li>Computational prioritization & immunogenicity scoring</li>
<li>Cross-genotype robustness analysis</li>
</ul>

<h3 style='color: #00d4ff; margin-top: 20px;'>Scientific Integrity</h3>
<ul style='color: #ffffff;'>
<li>No placeholder data</li>
<li>Strict validation filters</li>
<li>Transparent exclusion criteria</li>
</ul>
</div>
""", unsafe_allow_html=True)

# =========================================================
# PROFESSIONAL FOOTER
# =========================================================
st.markdown("""
<div class='footer-container'>
<div class='footer'>
© 2026 Abubakar Research Labs — 
<strong>VAXINTAIC v3 | PRRSV Vaccine Intelligence Platform</strong><br/>
Porcine Reproductive and Respiratory Syndrome Virus (Type I & Type II)<br/>
<small style='color: #8a2be2;'>Advanced Vaccine Intelligence Platform</small>
</div>
</div>
""", unsafe_allow_html=True)
