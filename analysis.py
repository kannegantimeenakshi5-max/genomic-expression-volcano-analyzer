import numpy as np
import pandas as pd
import plotly.express as px
from scipy import stats
import streamlit as st

st.set_page_config(
    page_title="Genomic Expression Profiler", layout="wide"
)
st.title(
    "🧬 Differential Gene Expression (DGE) & Volcano Plotter"
)


@st.cache_data
def generate_synthetic_microarray_data(
    n_genes=600, n_samples_ctrl=6, n_samples_treat=6
):
  np.random.seed(42)
  genes = [f"GENE_{i:04d}" for i in range(1, n_genes + 1)]

  # Baseline normal distributions
  ctrl = np.random.normal(
      loc=8.0, scale=1.2, size=(n_genes, n_samples_ctrl)
  )
  treat = np.random.normal(
      loc=8.0, scale=1.2, size=(n_genes, n_samples_treat)
  )

  # Introduce targeted differential expression for 60 genes
  treat[:30] += np.random.uniform(
      1.5, 3.5, size=(30, n_samples_treat)
  )  # Upregulated
  treat[30:60] -= np.random.uniform(
      1.5, 3.5, size=(30, n_samples_treat)
  )  # Downregulated

  ctrl_cols = [f"Control_{i+1}" for i in range(n_samples_ctrl)]
  treat_cols = [f"Treated_{i+1}" for i in range(n_samples_treat)]

  df = pd.DataFrame(
      np.hstack([ctrl, treat]),
      columns=ctrl_cols + treat_cols,
      index=genes,
  )
  return df, ctrl_cols, treat_cols


def run_dge_analysis(
    df, ctrl_cols, treat_cols, log2fc_thresh=1.0, pval_thresh=0.05
):
  mean_ctrl = df[ctrl_cols].mean(axis=1)
  mean_treat = df[treat_cols].mean(axis=1)

  log2fc = mean_treat - mean_ctrl
  p_vals = [
      stats.ttest_ind(df.loc[g, treat_cols], df.loc[g, ctrl_cols])[
          1
      ]
      for g in df.index
  ]

  results = pd.DataFrame({
      "Gene": df.index,
      "Log2FC": log2fc,
      "p_value": p_vals,
      "neg_log10_pval": -np.log10(p_vals),
  })

  conditions = [
      (results["Log2FC"] >= log2fc_thresh)
      & (results["p_value"] < pval_thresh),
      (results["Log2FC"] <= -log2fc_thresh)
      & (results["p_value"] < pval_thresh),
  ]
  choices = ["Upregulated", "Downregulated"]
  results["Regulation"] = np.select(
      conditions, choices, default="Not Significant"
  )

  return results


st.sidebar.header("Filter Thresholds")
fc_cut = st.sidebar.slider(
    "Log2 Fold-Change Cutoff",
    min_value=0.5,
    max_value=2.5,
    value=1.0,
    step=0.1,
)
p_cut = st.sidebar.slider(
    "P-value Significance Threshold",
    min_value=0.001,
    max_value=0.05,
    value=0.01,
    step=0.005,
)

data, c_cols, t_cols = generate_synthetic_microarray_data()
dge_results = run_dge_analysis(
    data, c_cols, t_cols, log2fc_thresh=fc_cut, pval_thresh=p_cut
)

col1, col2 = st.columns([2, 1])
with col1:
  st.subheader("🌋 Interactive Volcano Plot")
  fig = px.scatter(
      dge_results,
      x="Log2FC",
      y="neg_log10_pval",
      color="Regulation",
      hover_data=["Gene", "p_value"],
      color_discrete_map={
          "Upregulated": "#ef4444",
          "Downregulated": "#3b82f6",
          "Not Significant": "#94a3b8",
      },
      labels={
          "neg_log10_pval": "-log10(P-value)",
          "Log2FC": "Log2 Fold Change",
      },
  )
  fig.add_hline(
      y=-np.log10(p_cut), line_dash="dash", line_color="gray"
  )
  fig.add_vline(x=fc_cut, line_dash="dash", line_color="gray")
  fig.add_vline(x=-fc_cut, line_dash="dash", line_color="gray")
  st.plotly_chart(fig, use_container_width=True)

with col2:
  st.subheader("📋 Top Significant Biomarkers")
  sig_genes = dge_results[
      dge_results["Regulation"] != "Not Significant"
  ].sort_values(by="p_value")
  st.dataframe(
      sig_genes[["Gene", "Log2FC", "p_value", "Regulation"]].head(
          15
      ),
      use_container_width=True,
  )
