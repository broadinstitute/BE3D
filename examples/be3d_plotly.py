"""
BEClust3D Interactive (Plotly) Result Viewers
-----------------------------------------------
Reads the same tabular outputs the pipeline scripts already write to
`output_dir` (no changes to the core matplotlib plotting in the beclust3d
package) and renders interactive Plotly equivalents of the BE3D result
figures: hover tooltips, zoom/pan, and legend toggling.

Every function returns a `plotly.graph_objects.Figure`, or `None` (with a
printed warning) if the expected data file isn't present -- e.g. because a
step was skipped or output paths differ from the standard layout.

Usage in a notebook cell:
    from be3d_plotly import plot_hypothesis_qa
    fig = plot_hypothesis_qa(output_dir, pthr=0.05)
    if fig is not None:
        fig.show()

To lay two or more already-built figures out side by side instead of
stacked, use show_side_by_side():
    from be3d_plotly import plot_score_scatter, show_side_by_side
    fig1 = plot_score_scatter(output_dir, input_gene, screen_name, score_type='LFC')
    fig2 = plot_score_scatter(output_dir, input_gene, screen_name, score_type='LFC3D')
    show_side_by_side(fig1, fig2)
"""

import os
import glob
import pickle

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import ipywidgets as widgets
from IPython.display import display

# Default figure size (px), kept 1:1 square. Plotly figures otherwise
# auto-fit to the notebook cell/container width, which shrinks/stretches
# unpredictably between Colab and Jupyter. Every plotting function below
# accepts width=/height= overrides if you want a different size for a
# specific plot.
DEFAULT_WIDTH = 700
DEFAULT_HEIGHT = 700


def _find_one(pattern):
    matches = sorted(glob.glob(pattern))
    return matches[0] if matches else None


def _load_tsv(path):
    if path is None or not os.path.exists(path):
        return None
    return pd.read_csv(path, sep="\t")


def _warn(msg):
    print(f"[be3d_plotly] {msg}")


def show_side_by_side(*figs, width=None, height=None):
    """
    Display multiple already-built Plotly figures in one row instead of
    stacked, each staying independently interactive (hover/zoom/legend
    toggle) via ipywidgets.HBox + go.FigureWidget.

    Any None entries (e.g. a plot_* call that returned None because its
    input data was missing) are skipped. If only one figure remains after
    that, it's just shown normally at DEFAULT_WIDTH x DEFAULT_HEIGHT.

    width / height : optional per-panel override (px). Defaults to
    splitting DEFAULT_WIDTH evenly across however many figures are shown
    (min 350px per panel) and DEFAULT_HEIGHT for the height.
    """
    figs = [f for f in figs if f is not None]
    if not figs:
        return
    if len(figs) == 1:
        figs[0].update_layout(width=width or DEFAULT_WIDTH, height=height or DEFAULT_HEIGHT,
                               autosize=False)
        figs[0].show()
        return

    panel_width = width or max(DEFAULT_WIDTH // len(figs), 350)
    panel_height = height or DEFAULT_HEIGHT
    panels = []
    for fig in figs:
        fig.update_layout(width=panel_width, height=panel_height, autosize=False)
        panels.append(go.FigureWidget(fig))
    display(widgets.HBox(panels))


# ── 1. BE-QA hypothesis test scatter ──────────────────────────────────────────

def plot_hypothesis_qa(output_dir, pthr=0.05, test="MannWhitney",
                        width=DEFAULT_WIDTH, height=DEFAULT_HEIGHT):
    """
    Interactive scatter of hypothesis-test p-values by screen id.
    Reads `hypothesis_qa/{test}_hypothesis2.tsv` (test = "MannWhitney" or
    "KolmogorovSmirnov").

    The p-value column is named dynamically by the pipeline as
    `p_{cases}_vs_{controls}` (e.g. `p_Nonsense_vs_No Mutation`), not a
    fixed `p_CaseVsControl`, so it's detected by prefix instead of a
    hardcoded name.
    """
    path = os.path.join(output_dir, "hypothesis_qa", f"{test}_hypothesis2.tsv")
    df = _load_tsv(path)
    if df is None:
        _warn(f"Could not find {path}; skipping hypothesis QA plot.")
        return None

    p_cols = [c for c in df.columns if c.startswith("p_")]
    if not p_cols:
        _warn(f"No 'p_*' column found in {path}; skipping plot.")
        return None
    pcol = p_cols[0]
    comp_name = pcol[len("p_"):]

    df = df.copy()
    df[pcol] = pd.to_numeric(df[pcol], errors="coerce")
    df["-log10(p)"] = -np.log10(df[pcol].clip(lower=1e-300))
    df["significant"] = np.where(df[pcol] < pthr, f"p < {pthr}", f"p >= {pthr}")

    hover_cols = [c for c in ("gene_name", "num_of_cases", "num_of_controls") if c in df.columns]
    fig = px.scatter(
        df, x="screenid", y="-log10(p)", color="significant",
        hover_data=hover_cols + [pcol],
        title=f"BE-QA: {test} test, {comp_name} by screen",
        color_discrete_map={f"p < {pthr}": "#e74c3c", f"p >= {pthr}": "#7f8c8d"},
    )
    fig.add_hline(y=-np.log10(pthr), line_dash="dash", line_color="gray",
                   annotation_text=f"p = {pthr}")
    fig.update_layout(xaxis_title="Screen", yaxis_title="-log10(p-value)",
                       width=width, height=height, autosize=False)
    return fig


# ── 2. Violin plot of raw scores by mutation type ─────────────────────────────

def plot_violin_by_muttype(screen_dir, screen_file, mut_col="Mut_type", val_col="sgRNA_score",
                            gene_col=None, gene_name=None,
                            width=DEFAULT_WIDTH, height=DEFAULT_HEIGHT):
    """
    Interactive violin plot of the raw per-guide score distribution,
    split by mutation category. Reads the raw screen TSV directly
    (this data is never written back out by the pipeline).
    """
    path = os.path.join(screen_dir, screen_file)
    df = _load_tsv(path)
    if df is None:
        _warn(f"Could not find raw screen file {path}; skipping violin plot.")
        return None
    if mut_col not in df.columns or val_col not in df.columns:
        _warn(f"Expected columns '{mut_col}'/'{val_col}' not found in {path}; skipping plot.")
        return None

    if gene_col and gene_name and gene_col in df.columns:
        df = df[df[gene_col] == gene_name]

    fig = px.violin(
        df, x=mut_col, y=val_col, color=mut_col, box=True, points="outliers",
        title=f"Score distribution by mutation category — {os.path.basename(screen_file)}",
    )
    fig.update_layout(xaxis_title="Mutation category", yaxis_title=val_col, showlegend=False,
                       width=width, height=height, autosize=False)
    return fig


# ── 3. LFC / LFC3D cutoff scatter along sequence position ────────────────────

def plot_score_scatter(output_dir, input_gene, screen_name, score_type="LFC", pthr_str="05",
                        width=DEFAULT_WIDTH, height=DEFAULT_HEIGHT):
    """
    Interactive scatter of per-residue LFC or LFC3D scores along sequence
    position (unipos), colored by significance at the given p-value cutoff.
    Reads `{score_type}/*_{input_gene}_NonAggr_{score_type}.tsv`.
    """
    pattern = os.path.join(output_dir, score_type, f"*_{input_gene}_NonAggr_{score_type}.tsv")
    path = _find_one(pattern)
    df = _load_tsv(path)
    if df is None:
        _warn(f"Could not find NonAggr {score_type} table ({pattern}); skipping plot.")
        return None

    pos_val_col = f"{screen_name}_{score_type}_pos"
    neg_val_col = f"{screen_name}_{score_type}_neg"
    pos_sig_col = f"{screen_name}_{score_type}_pos_{pthr_str}_psig"
    neg_sig_col = f"{screen_name}_{score_type}_neg_{pthr_str}_psig"
    needed = [pos_val_col, neg_val_col, pos_sig_col, neg_sig_col]
    if "unipos" not in df.columns or not all(c in df.columns for c in needed):
        _warn(f"Expected columns {['unipos'] + needed} not all found in {path}; skipping plot.")
        return None

    # The pipeline uses the literal string '-' (not NaN) to mark residues
    # with no value in the pos/neg split; coerce to numeric so those become
    # proper NaN instead of tripping up abs()/comparisons below.
    df[pos_val_col] = pd.to_numeric(df[pos_val_col], errors="coerce")
    df[neg_val_col] = pd.to_numeric(df[neg_val_col], errors="coerce")

    rows = []
    for _, r in df.iterrows():
        if pd.notna(r[pos_val_col]) and r[pos_val_col] != 0:
            rows.append({"unipos": r["unipos"], "value": r[pos_val_col], "direction": "positive",
                          "significant": r[pos_sig_col]})
        if pd.notna(r[neg_val_col]) and r[neg_val_col] != 0:
            rows.append({"unipos": r["unipos"], "value": -abs(r[neg_val_col]), "direction": "negative",
                         "significant": r[neg_sig_col]})
    if not rows:
        _warn(f"No non-zero {score_type} values found for {screen_name}; skipping plot.")
        return None
    long_df = pd.DataFrame(rows)
    long_df["group"] = long_df["direction"] + " / " + long_df["significant"].astype(str)

    fig = px.scatter(
        long_df, x="unipos", y="value", color="group",
        title=f"{input_gene} {score_type} by sequence position — {screen_name}",
        color_discrete_map={
            f"positive / p<{pthr_str}": "#c0392b", f"positive / p>={pthr_str}": "#f5b7b1",
            f"negative / p<{pthr_str}": "#1f618d", f"negative / p>={pthr_str}": "#aed6f1",
        },
    )
    fig.add_hline(y=0, line_color="black", line_width=1)
    fig.update_layout(xaxis_title="Residue position", yaxis_title=score_type,
                       width=width, height=height, autosize=False)
    return fig


# ── 4. 3-D structural cluster viewer (replaces the static dendrogram) ────────

def plot_cluster_3d(output_dir, input_gene, screen_name, score_type="LFC3D",
                     direction="Positive", pthr_str="001", dist="6A",
                     width=DEFAULT_WIDTH, height=DEFAULT_HEIGHT):
    """
    Interactive 3-D scatter of residues colored by spatial cluster
    assignment, using the residue coordinates + cluster ids already saved
    to `cluster_{score_type}/{input_gene}_{screen_name}_Aggr_Hits.tsv`.
    This is a more explorable substitute for the static dendrogram image
    (rotate/zoom/hover instead of a flat tree).
    """
    path = os.path.join(output_dir, f"cluster_{score_type}", f"{input_gene}_{screen_name}_Aggr_Hits.tsv")
    df = _load_tsv(path)
    if df is None:
        _warn(f"Could not find {path}; skipping 3D cluster plot.")
        return None

    sign = "pos" if direction.lower().startswith("pos") else "neg"
    cluster_col = f"{score_type}_{sign}_{pthr_str}_psig_Clust_{dist}"
    coord_cols = ["x_coord", "y_coord", "z_coord"]
    if cluster_col not in df.columns or not all(c in df.columns for c in coord_cols):
        candidates = [c for c in df.columns if "_Clust_" in c]
        _warn(f"Expected cluster column '{cluster_col}' not found in {path}. "
              f"Available cluster columns: {candidates}. Skipping plot.")
        return None

    df = df.dropna(subset=coord_cols)
    df[cluster_col] = df[cluster_col].fillna("none").astype(str)

    fig = px.scatter_3d(
        df, x="x_coord", y="y_coord", z="z_coord", color=cluster_col,
        hover_data=[c for c in ("unipos", "unires", "chain") if c in df.columns],
        title=f"{input_gene} {direction} {score_type} clusters (p<{pthr_str}, {dist}) — {screen_name}",
    )
    fig.update_traces(marker=dict(size=4))
    fig.update_layout(width=width, height=height, autosize=False)
    return fig


# ── 5a. LFC vs LFC3D scatter ──────────────────────────────────────────────────

def plot_lfc_lfc3d_scatter(output_dir, input_gene, screen_name, pthr_str="05",
                            width=DEFAULT_WIDTH, height=DEFAULT_HEIGHT):
    lfc_path = _find_one(os.path.join(output_dir, "LFC", f"*_{input_gene}_LFC_dis_wght.tsv"))
    lfc3d_path = _find_one(os.path.join(output_dir, "LFC3D", f"*_{input_gene}_LFC3D_dis_wght.tsv"))
    psig_path = _find_one(os.path.join(output_dir, "LFC3D", f"*_{input_gene}_NonAggr_LFC3D.tsv"))
    lfc_df, lfc3d_df, psig_df = _load_tsv(lfc_path), _load_tsv(lfc3d_path), _load_tsv(psig_path)
    if lfc_df is None or lfc3d_df is None or psig_df is None:
        _warn("Could not find LFC/LFC3D dis_wght or NonAggr tables; skipping LFC vs LFC3D scatter.")
        return None

    lfc_col, lfc3d_col = f"{screen_name}_LFC", f"{screen_name}_LFC3D"
    pos_sig_col, neg_sig_col = f"{screen_name}_LFC3D_pos_{pthr_str}_psig", f"{screen_name}_LFC3D_neg_{pthr_str}_psig"
    if lfc_col not in lfc_df.columns or lfc3d_col not in lfc3d_df.columns:
        _warn(f"Expected columns '{lfc_col}'/'{lfc3d_col}' not found; skipping plot.")
        return None

    # Non-conserved / skipped residues are written as the literal string
    # '-' rather than NaN in these tables; coerce to numeric.
    lfc_df = lfc_df.copy()
    lfc3d_df = lfc3d_df.copy()
    lfc_df[lfc_col] = pd.to_numeric(lfc_df[lfc_col], errors="coerce")
    lfc3d_df[lfc3d_col] = pd.to_numeric(lfc3d_df[lfc3d_col], errors="coerce")

    merged = lfc_df[["unipos", lfc_col]].merge(lfc3d_df[["unipos", lfc3d_col]], on="unipos")
    merged = merged.dropna(subset=[lfc_col, lfc3d_col])
    if pos_sig_col in psig_df.columns and neg_sig_col in psig_df.columns:
        merged = merged.merge(psig_df[["unipos", pos_sig_col, neg_sig_col]], on="unipos", how="left")

        def _label(r):
            pos_hit = str(r.get(pos_sig_col, "")).startswith(f"p<{pthr_str}")
            neg_hit = str(r.get(neg_sig_col, "")).startswith(f"p<{pthr_str}")
            if pos_hit and neg_hit:
                return "Pos + Neg Hit"
            if pos_hit:
                return "Pos Hit"
            if neg_hit:
                return "Neg Hit"
            return "Not a Hit"

        merged["hit_type"] = merged.apply(_label, axis=1)
    else:
        merged["hit_type"] = "Not a Hit"

    fig = px.scatter(
        merged, x=lfc_col, y=lfc3d_col, color="hit_type", hover_data=["unipos"],
        title=f"{input_gene} LFC vs. LFC3D — {screen_name}",
        color_discrete_map={"Not a Hit": "#bdc3c7", "Pos Hit": "#e74c3c",
                             "Neg Hit": "#2980b9", "Pos + Neg Hit": "#8e44ad"},
    )
    fig.update_layout(xaxis_title="LFC", yaxis_title="LFC3D",
                       width=width, height=height, autosize=False)
    return fig


# ── 5b. pLDDT vs RSA scatter ──────────────────────────────────────────────────

def plot_plddt_rsa_scatter(output_dir, input_gene, screen_name,
                            width=DEFAULT_WIDTH, height=DEFAULT_HEIGHT):
    struc_path = _find_one(os.path.join(output_dir, "sequence_structure", "*_coord_struc_features.tsv"))
    lfc3d_path = _find_one(os.path.join(output_dir, "LFC3D", f"*_{input_gene}_LFC3D_dis_wght.tsv"))
    struc_df, lfc3d_df = _load_tsv(struc_path), _load_tsv(lfc3d_path)
    if struc_df is None or lfc3d_df is None:
        _warn("Could not find structural features or LFC3D dis_wght tables; skipping pLDDT vs RSA scatter.")
        return None

    lfc3d_col, wght_col = f"{screen_name}_LFC3D", f"{screen_name}_LFC3D_wght"
    needed = ["bfactor_pLDDT", "RSA"]
    if not all(c in struc_df.columns for c in needed) or lfc3d_col not in lfc3d_df.columns:
        _warn(f"Expected columns {needed + [lfc3d_col]} not all found; skipping plot.")
        return None

    # Non-conserved / skipped residues are written as the literal string
    # '-' rather than NaN in these tables; coerce to numeric.
    lfc3d_df = lfc3d_df.copy()
    lfc3d_df[lfc3d_col] = pd.to_numeric(lfc3d_df[lfc3d_col], errors="coerce")
    if wght_col in lfc3d_df.columns:
        lfc3d_df[wght_col] = pd.to_numeric(lfc3d_df[wght_col], errors="coerce")

    cols = ["unipos"] + [c for c in (lfc3d_col, wght_col) if c in lfc3d_df.columns]
    merged = struc_df[["unipos"] + needed].merge(lfc3d_df[cols], on="unipos")
    merged = merged.dropna(subset=[lfc3d_col])
    merged["direction"] = np.where(merged[lfc3d_col] >= 0, "positive", "negative")
    size_col = None
    if wght_col in merged.columns:
        merged["marker_size"] = merged[wght_col].abs() * 100
        size_col = "marker_size"

    fig = px.scatter(
        merged, x="bfactor_pLDDT", y="RSA", color="direction", size=size_col,
        hover_data=["unipos", lfc3d_col],
        title=f"{input_gene} pLDDT vs. RSA — {screen_name}",
        color_discrete_map={"positive": "#c0392b", "negative": "#1f618d"},
    )
    fig.update_layout(xaxis_title="pLDDT", yaxis_title="RSA",
                       width=width, height=height, autosize=False)
    return fig


# ── 5c. Hit-count bar plots by domain / pLDDT disorder category ──────────────

def _hit_counts_by_category(output_dir, input_gene, screen_name, category_df, category_col, pthr_str):
    nonaggr_path = _find_one(os.path.join(output_dir, "LFC3D", f"*_{input_gene}_NonAggr_LFC3D.tsv"))
    nonaggr_df = _load_tsv(nonaggr_path)
    if nonaggr_df is None or category_col not in category_df.columns:
        return None

    pos_sig_col, neg_sig_col = f"{screen_name}_LFC3D_pos_{pthr_str}_psig", f"{screen_name}_LFC3D_neg_{pthr_str}_psig"
    if pos_sig_col not in nonaggr_df.columns or neg_sig_col not in nonaggr_df.columns:
        return None

    merged = category_df[["unipos", category_col]].merge(
        nonaggr_df[["unipos", pos_sig_col, neg_sig_col]], on="unipos", how="left")
    pos_hits = merged[merged[pos_sig_col].astype(str).str.startswith(f"p<{pthr_str}")]
    neg_hits = merged[merged[neg_sig_col].astype(str).str.startswith(f"p<{pthr_str}")]

    counts = pd.concat([
        pos_hits.groupby(category_col).size().rename("count").reset_index().assign(hit_type="POS"),
        neg_hits.groupby(category_col).size().rename("count").reset_index().assign(hit_type="NEG"),
    ], ignore_index=True)
    return counts


def plot_domain_barplot(output_dir, input_gene, input_uniprot, screen_name, pthr_str="05",
                         width=DEFAULT_WIDTH, height=DEFAULT_HEIGHT):
    domains_path = os.path.join(output_dir, "sequence_structure", f"{input_gene}_{input_uniprot}_domains.tsv")
    domains_df = _load_tsv(domains_path)
    if domains_df is None:
        _warn(f"Could not find {domains_path}; skipping domain barplot.")
        return None

    # domains.tsv uses 'Position' for the residue index (not 'unipos' like
    # the other structure/score tables); normalize so it can be merged.
    if "unipos" not in domains_df.columns and "Position" in domains_df.columns:
        domains_df = domains_df.rename(columns={"Position": "unipos"})

    counts = _hit_counts_by_category(output_dir, input_gene, screen_name, domains_df, "Domain", pthr_str)
    if counts is None or counts.empty:
        _warn("Could not assemble domain hit counts; skipping domain barplot.")
        return None

    fig = px.bar(
        counts, x="Domain", y="count", color="hit_type", barmode="group",
        title=f"{input_gene} LFC3D hit count by domain (p<{pthr_str}) — {screen_name}",
        color_discrete_map={"POS": "#c0392b", "NEG": "#1f618d"},
    )
    fig.update_layout(xaxis_title="Domain", yaxis_title="Hit count",
                       width=width, height=height, autosize=False)
    return fig


def plot_plddt_dis_barplot(output_dir, input_gene, screen_name, pthr_str="05",
                            width=DEFAULT_WIDTH, height=DEFAULT_HEIGHT):
    struc_path = _find_one(os.path.join(output_dir, "sequence_structure", "*_coord_struc_features.tsv"))
    struc_df = _load_tsv(struc_path)
    if struc_df is None:
        _warn(f"Could not find structural features table; skipping pLDDT-disorder barplot.")
        return None

    counts = _hit_counts_by_category(output_dir, input_gene, screen_name, struc_df, "pLDDT_dis", pthr_str)
    if counts is None or counts.empty:
        _warn("Could not assemble pLDDT-disorder hit counts; skipping barplot.")
        return None

    fig = px.bar(
        counts, x="pLDDT_dis", y="count", color="hit_type", barmode="group",
        title=f"{input_gene} LFC3D hit count by pLDDT-disorder category (p<{pthr_str}) — {screen_name}",
        color_discrete_map={"POS": "#c0392b", "NEG": "#1f618d"},
    )
    fig.update_layout(xaxis_title="pLDDT disorder category", yaxis_title="Hit count",
                       width=width, height=height, autosize=False)
    return fig


# ── 5d. Enrichment test forest plot ──────────────────────────────────────────

def plot_enrichment_test(output_dir, input_gene, screen_name, score_type="LFC3D", pthr_str="05",
                          width=DEFAULT_WIDTH, height=DEFAULT_HEIGHT):
    """
    Reads the enrichment_test() pickle (list of dicts with score_type,
    odds_ratio [log2], ci [log2 low/high], p_value) and renders a forest
    plot with hoverable confidence intervals.
    """
    path = os.path.join(output_dir, "characterization",
                         f"{input_gene}_enrichment_test_{score_type}_{pthr_str}_{screen_name}.pickle")
    if not os.path.exists(path):
        _warn(f"Could not find {path}; skipping enrichment test plot.")
        return None
    with open(path, "rb") as fh:
        records = pickle.load(fh)
    df = pd.DataFrame(records)
    if df.empty or "odds_ratio" not in df.columns:
        _warn(f"No usable records in {path}; skipping enrichment test plot.")
        return None

    df["ci_low"] = df["ci"].apply(lambda c: c[0])
    df["ci_high"] = df["ci"].apply(lambda c: c[1])
    y_labels = df["score_type"] if "score_type" in df.columns else df.index.astype(str)

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=df["odds_ratio"], y=y_labels,
        mode="markers",
        error_x=dict(type="data", symmetric=False,
                     array=df["ci_high"] - df["odds_ratio"], arrayminus=df["odds_ratio"] - df["ci_low"]),
        marker=dict(size=10, color="#2c3e50"),
        text=[f"p = {p:.3g}" for p in df.get("p_value", [None] * len(df))],
        hovertemplate="%{y}: log2(OR) = %{x:.2f}<br>%{text}<extra></extra>",
    ))
    fig.add_vline(x=0, line_dash="dash", line_color="gray")
    fig.update_layout(
        title=f"{input_gene} enrichment test ({score_type}, p<{pthr_str}) — {screen_name}",
        xaxis_title="log2(odds ratio)", yaxis_title="",
        width=width, height=height, autosize=False,
    )
    return fig
