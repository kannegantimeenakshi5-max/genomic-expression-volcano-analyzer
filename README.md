Genomic Expression Volcano Analyzer
![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)
![Python](https://img.shields.io/badge/python-3.x-blue)
![SciPy](https://img.shields.io/badge/SciPy-stats-orange)
![Plotly](https://img.shields.io/badge/Plotly-interactive-lightgrey)
An interactive bioinformatics dashboard for differential gene expression analysis. It runs Student's t-tests and calculates log2 fold changes on a transcriptomic matrix, then shows the results as a dynamic, publication-ready volcano plot so candidate biomarker genes stand out.
Live demo: https://genomic-expression-volcano-analyzer-syylmldmnztjmkywn36whs.streamlit.app/
<!-- TODO: add a screenshot of a volcano plot. Save it as docs/volcano.png, then uncomment the next line.
![Volcano plot](docs/volcano.png)
---
Background
Differential expression analysis asks which genes change between two conditions, for example tumour vs normal tissue or treated vs control. Two numbers answer that:
Log2 fold change (log2FC): how large the change is. A value of +1 means expression doubled, and -1 means it halved.
p-value: how likely a change that large would be to appear by chance.
A volcano plot puts both on one chart: log2FC on the x-axis and -log10(p-value) on the y-axis. Genes in the upper left and upper right corners change a lot and are statistically significant, so they are the candidate biomarkers. Typical cut-offs are |log2FC| ≥ 1 and p < 0.05, though the right choice depends on the study.
Features
Differential expression testing with Student's t-test (SciPy)
Log2 fold-change calculation between two groups
Interactive volcano plot built with Plotly (hover, zoom, pan)
Highlights significantly up- and down-regulated genes as candidate biomarkers
Runs in the browser through Streamlit
<!-- TODO: confirm against analysis.py and add or remove bullets:
     - input format (CSV of genes x samples? how are the two groups specified?)
     - whether log2FC / p-value thresholds are adjustable in the UI
     - whether a ranked table of significant genes is shown or downloadable
How it works
Load an expression matrix (genes × samples).
Split samples into the two groups being compared.
For each gene, compute the mean expression in each group and the log2 fold change.
For each gene, run a Student's t-test between the groups.
Plot log2FC against -log10(p-value) and colour genes by significance.
<!-- TODO: edit these steps so they match analysis.py exactly. -->
Run it locally
```bash
git clone https://github.com/kannegantimeenakshi5-max/genomic-expression-volcano-analyzer.git
cd genomic-expression-volcano-analyzer
pip install -r requirements.txt
streamlit run analysis.py
```
Limitations
A t-test per gene is a screening step. Testing thousands of genes at once produces false positives unless p-values are corrected for multiple testing (for example Benjamini-Hochberg FDR).
Student's t-test assumes roughly normal data with similar variance in both groups. For raw count data, dedicated tools such as DESeq2 or edgeR are the standard choice.
Results identify candidates for follow-up, not confirmed biomarkers.
<!-- TODO: if analysis.py already applies FDR correction, delete the first bullet and say so in Features.
Tech stack
Python · SciPy · Plotly · Streamlit
Author
Kanneganti Meenakshi Naga Prasanna
B.Tech Bioinformatics, Vignan University
GitHub
