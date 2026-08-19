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
        display(fig)

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
from plotly.subplots import make_subplots
from scipy.cluster.hierarchy import dendrogram as _scipy_dendrogram
from sklearn.cluster import AgglomerativeClustering
from IPython.display import display

# Default figure size (px), kept 1:1 square. Plotly figures otherwise
# auto-fit to the notebook cell/container width, which shrinks/stretches
# unpredictably between Colab and Jupyter. Every plotting function below
# accepts width=/height= overrides if you want a different size for a
# specific plot.
DEFAULT_WIDTH = 1400
DEFAULT_HEIGHT = 600

# Shared color scheme, used consistently for every p-value / direction
# based coloring in this module: positive selection is blue, negative
# selection is red. Wherever a plot splits the same direction further by
# significance (p<thr vs p>=thr), the non-significant half uses the pale
# version of that direction's color rather than a separate hue.
COLOR_POS = "#1f618d"        # blue: positive direction, significant
COLOR_POS_PALE = "#aed6f1"   # pale blue: positive direction, not significant
COLOR_NEG = "#c0392b"        # red: negative direction, significant
COLOR_NEG_PALE = "#f5b7b1"   # pale red: negative direction, not significant
COLOR_NEUTRAL = "#bdc3c7"    # gray: no clear direction / not a hit
COLOR_COMBINED = "#8e44ad"   # purple: significant in both directions

_HIT_TYPE_COLORS = {"Not a Hit": COLOR_NEUTRAL, "Pos Hit": COLOR_POS,
                     "Neg Hit": COLOR_NEG, "Pos + Neg Hit": COLOR_COMBINED}


def _find_one(pattern):
    matches = sorted(glob.glob(pattern))
    return matches[0] if matches else None


def _load_tsv(path):
    if path is None or not os.path.exists(path):
        return None
    return pd.read_csv(path, sep="\t")


def _warn(msg):
    print(f"[be3d_plotly] {msg}")


def _pthr_label(pthr_str):
    """
    The pipeline's psig columns store significance labels as the literal
    string 'p<{threshold}' / 'p>={threshold}' using the full decimal
    representation (e.g. 'p<0.05', written by calculate_stats()), while
    column *names* use only the fractional digits (e.g.
    '..._pos_05_psig', written by nonaggregate()). pthr_str here is that
    column-name fragment ('05', '001', ...); this reconstructs the
    decimal string ('05' -> '0.05') so label comparisons actually match
    what's stored in the data instead of silently never matching.
    """
    return f"0.{pthr_str}"


def _axis_key(prefix, n):
    return prefix if n == 1 else f"{prefix}{n}"


def _fig_ncols(fig):
    """Counts a figure's own subplot columns (1 for an ordinary single-panel figure,
    2+ for one already built with make_subplots, e.g. _lfc_lfc3d_two_panel_figure's
    'No LFC' strip + main scatter) by finding how many numbered xaxes it has."""
    n = 1
    while _axis_key("xaxis", n + 1) in fig.layout:
        n += 1
    return n


def show_side_by_side(*figs, width=None, height=None, spacing=0.04):
    """
    Display multiple already-built Plotly figures in one row, merged into a single
    make_subplots figure and shown with display() -- deliberately NOT
    go.FigureWidget + ipywidgets.HBox (the previous approach), since FigureWidget
    requires the 'jupyterlab-plotly' comm widget to be registered with the notebook
    frontend; when it isn't (observed failing in a VS Code Jupyter session: "Failed to
    load model class 'FigureModel' ... No version of module jupyterlab-plotly is
    registered"), every side-by-side plot in the notebook breaks. A merged static
    figure only needs the plain Plotly renderer, which was already working.

    Each input figure's traces, axis titles/ranges/tick visibility, and 0-reference
    lines (added via add_hline/add_vline) are carried over into their own subplot
    column(s). An input figure may itself already be a multi-column subplot figure
    (e.g. plot_lfc_lfc3d_scatter's 'No LFC' strip + main scatter) -- its own columns
    are kept as a contiguous block with their original relative widths preserved,
    rather than assuming every input is a single panel. Traces sharing the same name
    across panels (e.g. 'positive / p<0.05' in two panels being compared) share one
    legendgroup, so there's one legend entry, not a duplicate one per panel, and
    toggling it hides/shows that name in every panel at once.

    Any None entries (e.g. a plot_* call that returned None because its input data
    was missing) are skipped. If only one figure remains after that, it's just shown
    normally at DEFAULT_WIDTH x DEFAULT_HEIGHT.

    width / height : optional override (px) applied uniformly to every input figure's
    overall share of the row, same as before this handled multi-column inputs. When
    left None (the common case), each figure's OWN pre-set width (whatever its plot_*
    call actually produced, e.g. plot_lfc_lfc3d_scatter's intentional half-width
    default) is honored instead of being silently replaced by an even DEFAULT_WIDTH
    split -- the combined row's total width is just the sum of those.

    spacing : fraction (0-1) of the total row width reserved as blank gap BETWEEN
    separate input figures (not within one figure's own internal columns, e.g. the
    'No LFC' strip stays tight against its main scatter regardless of this value) --
    default 0.04. Spread evenly across however many figure-to-figure gaps there are.
    """
    figs = [f for f in figs if f is not None]
    if not figs:
        return
    if len(figs) == 1:
        figs[0].update_layout(width=width or figs[0].layout.width or DEFAULT_WIDTH,
                               height=height or figs[0].layout.height or DEFAULT_HEIGHT,
                               autosize=False)
        display(figs[0])
        return

    ncols_list = [_fig_ncols(f) for f in figs]
    n_gaps = len(figs) - 1
    fig_widths = [width or fig.layout.width or max(DEFAULT_WIDTH // len(figs), 350) for fig in figs]
    total_width = sum(fig_widths)
    # All panels share one row, so they must share one height -- take the largest of each
    # figure's own set height rather than a blanket DEFAULT_HEIGHT, so a plot_* call with
    # an explicit height=... still gets it once composed into a row via show_side_by_side.
    panel_height = height or max((fig.layout.height or DEFAULT_HEIGHT) for fig in figs)

    # Each input figure keeps its own overall width share (fig_widths, above) out of the
    # (1 - spacing) content budget; within that share, its own columns keep their
    # original relative proportions (from their own domain widths). A figure built
    # directly via go.Figure() (e.g. plot_enrichment_test) rather than make_subplots
    # never had its xaxis.domain set at all (None, not (0, 1)), so fall back to a
    # full-width (0, 1) domain for those single-column figures. A dedicated empty spacer
    # column (hidden axes, no traces) sits between each figure's block -- reusing
    # make_subplots' own horizontal_spacing for this instead would also widen the gap
    # inside a multi-column figure's own panels, which should stay tight.
    if not 0.0 <= spacing < 1.0:
        raise ValueError(
            f"spacing must be a fraction of the total row width in [0, 1), got {spacing!r} -- "
            f"it's not pixels (unlike width/height); try something like spacing=0.08."
        )
    content_budget = 1.0 - spacing
    gap_share = (spacing / n_gaps) if n_gaps else 0.0

    column_widths = []
    block_start_col = []  # 1-indexed column each figure's own block starts at
    col = 1
    for i, (fig, n, fig_width) in enumerate(zip(figs, ncols_list, fig_widths)):
        block_start_col.append(col)
        raw_widths = []
        for c in range(1, n + 1):
            domain = fig.layout[_axis_key("xaxis", c)].domain or (0, 1)
            raw_widths.append(domain[1] - domain[0])
        total_raw = sum(raw_widths) or 1.0
        share = fig_width / total_width * content_budget
        column_widths.extend([w / total_raw * share for w in raw_widths])
        col += n
        if i < n_gaps:
            column_widths.append(gap_share)
            col += 1

    total_cols = len(column_widths)
    spacer_cols = [block_start_col[i] + ncols_list[i] for i in range(n_gaps)]

    subplot_titles = [""] * total_cols
    for fig, n, start in zip(figs, ncols_list, block_start_col):
        title = (fig.layout.title.text or "") if fig.layout.title else ""
        subplot_titles[start - 1] = title

    combined = make_subplots(rows=1, cols=total_cols, column_widths=column_widths,
                              subplot_titles=subplot_titles, horizontal_spacing=0.02)
    for spacer_col in spacer_cols:
        combined.update_xaxes(visible=False, row=1, col=spacer_col)
        combined.update_yaxes(visible=False, row=1, col=spacer_col)

    seen_legend_names = set()
    for fig, n, col_start in zip(figs, ncols_list, block_start_col):
        col_offset = col_start - 1
        for trace in fig.data:
            xaxis_id = trace.xaxis or "x"
            local_col = 1 if xaxis_id == "x" else int(xaxis_id[1:])
            target_col = col_offset + local_col

            t = trace.__class__(trace.to_plotly_json())  # copy -- don't mutate the caller's figure
            if t.name:
                t.legendgroup = t.name
                t.showlegend = t.name not in seen_legend_names
                seen_legend_names.add(t.name)
            combined.add_trace(t, row=1, col=target_col)

        for local_col in range(1, n + 1):
            target_col = col_offset + local_col
            xaxis = fig.layout[_axis_key("xaxis", local_col)]
            yaxis = fig.layout[_axis_key("yaxis", local_col)]
            if xaxis.title and xaxis.title.text:
                combined.update_xaxes(title_text=xaxis.title.text, row=1, col=target_col)
            if xaxis.showticklabels is False:
                combined.update_xaxes(showticklabels=False, row=1, col=target_col)
            if xaxis.range:
                combined.update_xaxes(range=list(xaxis.range), row=1, col=target_col)
            if yaxis.title and yaxis.title.text:
                combined.update_yaxes(title_text=yaxis.title.text, row=1, col=target_col)

        for shape in (fig.layout.shapes or []):
            xref = (shape.xref or "x").replace(" domain", "")
            local_col = 1 if xref in ("x", "") else int(xref[1:])
            target_col = col_offset + local_col
            if shape.y0 == shape.y1:
                combined.add_hline(y=shape.y0, row=1, col=target_col, line=shape.line)
            elif shape.x0 == shape.x1:
                combined.add_vline(x=shape.x0, row=1, col=target_col, line=shape.line)

    combined.update_layout(width=total_width, height=panel_height, autosize=False)
    display(combined)


def show_stacked(*figs):
    """
    Display multiple already-built Plotly figures one per row (plain sequential
    display() calls) instead of side by side. Used for dendrograms, which are wide and
    have many small leaf tick labels -- squeezing two into a shared row makes both
    unreadable, so each gets the page's full width and its own row.

    Any None entries (e.g. a plot_* call that returned None because its input data was
    missing) are skipped.
    """
    for fig in figs:
        if fig is not None:
            display(fig)


def show_picker(options, description="Select:"):
    """
    Display ONE already-built Plotly figure at a time out of `options`
    ({label: Figure or None}), with an ipywidgets Dropdown to switch between them --
    e.g. for dendrograms, where showing every score-type/direction/mode combination at
    once (stacked or side by side) is a lot to scan just to compare two of them, but a
    single "pick one to look at, then pick another" viewer stays legible. A label whose
    figure is None (e.g. a direction with no significant residues) shows a plain
    "not available" message instead of a blank plot.
    """
    import ipywidgets as widgets

    def _show(choice):
        fig = options.get(choice)
        if fig is not None:
            display(fig)
        else:
            print(f"[not available] {choice}")

    widgets.interact(_show, choice=widgets.Dropdown(options=list(options), description=description))


def show_figure_dropdown(options, description="Select:", width=None, height=None):
    """
    Display ONE already-built Plotly figure at a time out of `options`
    ({label: Figure or None}), with a dropdown to switch between them -- same purpose as
    show_picker(), but the dropdown is Plotly's OWN layout.updatemenus control rather than
    an ipywidgets Dropdown, and switching is handled client-side by toggling trace
    visibility. Nothing is re-executed, so no cell re-run is needed to change the
    selection.

    This exists because show_picker() (and ipywidgets generally) cannot render Plotly
    figures in Colab: interact() captures its callback's display() calls into an
    ipywidgets Output widget, and Colab's Plotly renderer only injects its JS into the
    cell's normal output stream, never into a widget's separate output area -- so those
    figures came up blank or flashed once and vanished. An updatemenus dropdown lives
    inside the figure itself, so it needs no widget frontend at all and behaves the same
    in Colab, Jupyter, and VS Code.

    All variants' traces are added to one figure up front (only the first visible), so
    every variant is computed once when this is called. Per-variant layout that actually
    differs between dendrograms -- title, and the x tick positions/labels, since each
    score-type/direction has its own set of significant residues as leaves -- is carried
    in each button's own layout update.

    Labels whose figure is None (e.g. a direction with no significant residues) are left
    out of the dropdown and reported once in a printed note, rather than being selectable
    and showing an empty plot.
    """
    available = {label: fig for label, fig in options.items() if fig is not None}
    missing = [label for label, fig in options.items() if fig is None]
    if missing:
        _warn(f"not available (no residues pass the significance threshold): {', '.join(missing)}")
    if not available:
        _warn("nothing to show.")
        return
    if len(available) == 1:
        display(next(iter(available.values())))
        return

    labels = list(available)
    figs = [available[label] for label in labels]
    trace_counts = [len(fig.data) for fig in figs]
    total_traces = sum(trace_counts)

    combined = go.Figure()
    for fig_idx, fig in enumerate(figs):
        for trace in fig.data:
            t = trace.__class__(trace.to_plotly_json())  # copy -- don't mutate the caller's figure
            t.visible = (fig_idx == 0)
            combined.add_trace(t)

    # Base layout comes from the first variant (shared: axis titles, threshold line/annotation,
    # size), then each button overrides only what genuinely differs per variant.
    first = figs[0]
    combined.update_layout(first.layout)
    combined.update_layout(
        width=width or first.layout.width or DEFAULT_WIDTH,
        height=height or first.layout.height or DEFAULT_HEIGHT,
        autosize=False, showlegend=False,
    )

    buttons = []
    trace_start = 0
    for fig_idx, (label, fig) in enumerate(zip(labels, figs)):
        visible = [False] * total_traces
        for i in range(trace_start, trace_start + trace_counts[fig_idx]):
            visible[i] = True
        trace_start += trace_counts[fig_idx]

        layout_update = {"title.text": (fig.layout.title.text or "") if fig.layout.title else ""}
        xaxis = fig.layout.xaxis
        if xaxis is not None:
            if xaxis.tickvals is not None:
                layout_update["xaxis.tickvals"] = list(xaxis.tickvals)
            if xaxis.ticktext is not None:
                layout_update["xaxis.ticktext"] = list(xaxis.ticktext)
            if xaxis.tickangle is not None:
                layout_update["xaxis.tickangle"] = xaxis.tickangle

        buttons.append(dict(label=label, method="update", args=[{"visible": visible}, layout_update]))

    combined.update_layout(updatemenus=[dict(
        active=0, buttons=buttons, direction="down", showactive=True,
        x=0.0, xanchor="left", y=1.18, yanchor="top",
    )], margin=dict(t=140))
    combined.add_annotation(
        text=description, showarrow=False, xref="paper", yref="paper",
        x=0.0, xanchor="left", y=1.26, yanchor="bottom", font=dict(size=13),
    )
    display(combined)


# ── 1. BE-QA hypothesis test scatter ──────────────────────────────────────────

# stat_prefix: the pipeline names the test-statistic column '{prefix}_{cases}_vs_{controls}',
# alongside 'p_{cases}_vs_{controls}' for the p-value -- D for Kolmogorov-Smirnov, U for
# Mann-Whitney (calculate_stats()'s own naming, scipy's ks_2samp/mannwhitneyu statistic names).
_TEST_INFO = {
    "MannWhitney": {"stat_prefix": "U", "label": "Mann-Whitney U test"},
    "KolmogorovSmirnov": {"stat_prefix": "D", "label": "Kolmogorov-Smirnov test"},
}


def plot_hypothesis_qa(output_dir, pthr=0.05, test="MannWhitney", hypothesis=2,
                        width=DEFAULT_WIDTH, height=DEFAULT_HEIGHT):
    """
    Interactive scatter of the hypothesis test's own statistic (KS's D, or MW's U) on
    the x-axis against -log10(p) on the y-axis, one point per screen/gene comparison.
    Reads `hypothesis_qa/{test}_hypothesis{hypothesis}.tsv` (test = "MannWhitney" or
    "KolmogorovSmirnov"; hypothesis = 1 or 2).

    The p-value and statistic columns are named dynamically by the pipeline as
    `p_{cases}_vs_{controls}` / `{D,U}_{cases}_vs_{controls}` (e.g.
    `p_Nonsense_vs_No Mutation`), not a fixed name, so the p column is detected by
    prefix and the statistic column derived from it rather than hardcoded.
    """
    path = os.path.join(output_dir, "hypothesis_qa", f"{test}_hypothesis{hypothesis}.tsv")
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

    test_info = _TEST_INFO.get(test, {"stat_prefix": None, "label": test})
    stat_col = f"{test_info['stat_prefix']}_{comp_name}" if test_info["stat_prefix"] else None
    if stat_col is None or stat_col not in df.columns:
        _warn(f"Expected statistic column '{stat_col}' not found in {path}; skipping plot.")
        return None

    df = df.copy()
    df[pcol] = pd.to_numeric(df[pcol], errors="coerce")
    df[stat_col] = pd.to_numeric(df[stat_col], errors="coerce")
    df["-log10(p)"] = -np.log10(df[pcol].clip(lower=1e-300))
    df["significant"] = np.where(df[pcol] < pthr, f"p < {pthr}", f"p >= {pthr}")

    hover_cols = [c for c in ("screenid", "gene_name", "num_of_cases", "num_of_controls") if c in df.columns]
    fig = px.scatter(
        df, x=stat_col, y="-log10(p)", color="significant",
        hover_data=hover_cols + [pcol],
        title=f"{test_info['label']} — Hypothesis {hypothesis}",
        color_discrete_map={f"p < {pthr}": COLOR_NEG, f"p >= {pthr}": COLOR_NEG_PALE},
    )
    fig.add_hline(y=-np.log10(pthr), line_dash="dash", line_color="gray",
                   annotation_text=f"p = {pthr}")
    fig.update_layout(xaxis_title=f"{test_info['stat_prefix']} statistic", yaxis_title="-log10(p-value)",
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


def plot_violin_by_processed_muttype(output_dir, input_gene, screen_name,
                                      width=DEFAULT_WIDTH, height=DEFAULT_HEIGHT):
    """
    Interactive violin plot of the PROCESSED per-guide LFC distribution, split by
    mutation category -- unlike plot_violin_by_muttype (which reads the raw screen file's
    mut_col as-is), this reads parse_be_data's own output,
    `screendata/{input_gene}_{screen_name}_{category}.tsv` (one file per category), so the
    categories/values shown already reflect mutation_priority collapsing (a guide with a
    semicolon-joined multi-category edit list resolved to one category) and each
    category's own filtering (e.g. Missense requires refAA != altAA, Nonsense requires
    altAA == '*') -- exactly what the rest of the pipeline (parse_be_data onward) sees.
    """
    pattern = os.path.join(output_dir, "screendata", f"{input_gene}_{screen_name}_*.tsv")
    paths = sorted(glob.glob(pattern))
    prefix = f"{input_gene}_{screen_name}_"
    rows = []
    for path in paths:
        category = os.path.basename(path)[len(prefix):-len(".tsv")]
        df = _load_tsv(path)
        if df is None or "LFC" not in df.columns:
            continue
        rows.append(pd.DataFrame({"category": category, "LFC": df["LFC"]}))
    if not rows:
        _warn(f"Could not find any {pattern}; skipping processed violin plot.")
        return None

    long_df = pd.concat(rows, ignore_index=True)
    fig = px.violin(
        long_df, x="category", y="LFC", color="category", box=True, points="outliers",
        title=f"{input_gene} processed LFC distribution by mutation category",
    )
    fig.update_layout(xaxis_title="Mutation category", yaxis_title="LFC", showlegend=False,
                       width=width, height=height, autosize=False)
    return fig


# ── 3. LFC / LFC3D cutoff scatter along sequence position ────────────────────

def plot_score_scatter(output_dir, input_gene, screen_name, score_type="LFC", direction=None, pthr_str="05",
                        width=DEFAULT_WIDTH, height=DEFAULT_HEIGHT):
    """
    Interactive scatter of per-residue LFC or LFC3D scores along sequence
    position (unipos), colored by significance at the given p-value cutoff.
    Reads `{score_type}/*_{input_gene}_NonAggr_{score_type}.tsv`.

    direction : "positive", "negative", or None (both, the previous default) --
        BE-Clust3D shows these as two separate plots rather than one combined
        scatter, so a hit in one direction doesn't visually compete for attention
        with (or get mistaken for) a hit in the other.
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
        if direction in (None, "positive") and pd.notna(r[pos_val_col]) and r[pos_val_col] != 0:
            rows.append({"unipos": r["unipos"], "value": r[pos_val_col], "direction": "positive",
                          "significant": r[pos_sig_col]})
        if direction in (None, "negative") and pd.notna(r[neg_val_col]) and r[neg_val_col] != 0:
            rows.append({"unipos": r["unipos"], "value": -abs(r[neg_val_col]), "direction": "negative",
                         "significant": r[neg_sig_col]})
    if not rows:
        _warn(f"No non-zero {direction or ''} {score_type} values found for {screen_name}; skipping plot.")
        return None
    long_df = pd.DataFrame(rows)
    long_df["group"] = long_df["direction"] + " / " + long_df["significant"].astype(str)

    # color_discrete_map keys must match the actual 'p<0.05'/'p>=0.05'
    # label text stored in the psig columns, not the '05' file-name suffix.
    pthr_label = _pthr_label(pthr_str)
    title_suffix = f" ({direction})" if direction else ""
    fig = px.scatter(
        long_df, x="unipos", y="value", color="group",
        title=f"{input_gene} {score_type}{title_suffix} by sequence position",
        color_discrete_map={
            f"positive / p<{pthr_label}": COLOR_POS, f"positive / p>={pthr_label}": COLOR_POS_PALE,
            f"negative / p<{pthr_label}": COLOR_NEG, f"negative / p>={pthr_label}": COLOR_NEG_PALE,
        },
    )
    fig.add_hline(y=0, line_color="black", line_width=1)
    fig.update_layout(xaxis_title="Residue position", yaxis_title=score_type,
                       width=width, height=height, autosize=False)
    return fig


# ── 4. 3-D structural cluster viewer (replaces the static dendrogram) ────────

def plot_cluster_3d(output_dir, input_gene, screen_name, score_type="LFC3D",
                     direction="Positive", pthr_str="05", dist="6A", column_prefix=None,
                     width=DEFAULT_WIDTH, height=DEFAULT_HEIGHT):
    """
    Interactive 3-D scatter of residues colored by spatial cluster
    assignment, using the residue coordinates + cluster ids already saved
    to `cluster_{score_type}/{input_gene}_{screen_name}_Aggr_Hits.tsv`.
    This is a more explorable substitute for the static dendrogram image
    (rotate/zoom/hover instead of a flat tree).

    column_prefix : str, optional
        Prefix used when clustering() built the psig column it clustered on
        (clustering.py names the output column '{psig_column}_Clust_{dist}A').
        For a per-screen call this is the screen name itself (the default,
        when left None), but for the meta-aggregate call it's the
        aggregation function name instead (e.g. 'SUM', 'mean') -- pass
        function_for_meta explicitly for screen_name='Meta'.
    """
    path = os.path.join(output_dir, f"cluster_{score_type}", f"{input_gene}_{screen_name}_Aggr_Hits.tsv")
    df = _load_tsv(path)
    if df is None:
        _warn(f"Could not find {path}; skipping 3D cluster plot.")
        return None

    if column_prefix is None:
        column_prefix = screen_name
    sign = "pos" if direction.lower().startswith("pos") else "neg"
    cluster_col = f"{column_prefix}_{score_type}_{sign}_{pthr_str}_psig_Clust_{dist}"
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
        title=f"{input_gene} {direction} {score_type} clusters (p<{pthr_str}, {dist})",
    )
    fig.update_traces(marker=dict(size=4))
    fig.update_layout(width=width, height=height, autosize=False)
    return fig


# ── 4b. Interactive dendrogram (the actual hierarchical-clustering view) ─────
#
# The pipeline never persists the full merge tree, only the rendered
# matplotlib image (beclust3d.lfc3d.clustering_plot.plot_dendrogram) -- the
# saved Aggr_Hits.tsv only has cluster *labels* at a handful of discrete
# distance cutoffs. So this recomputes the same single-linkage Euclidean
# AgglomerativeClustering(distance_threshold=dist) fit on the same
# significant-residue coordinates the pipeline used, which is deterministic
# and reproduces the identical tree/labels, then renders it as an actual
# Plotly dendrogram (not a 3-D scatter substitute).

def _agglomerative_dendrogram_figure(coords, leaf_labels, dist, title,
                                      width=DEFAULT_WIDTH, height=600):
    n = coords.shape[0]
    if n < 2:
        _warn("Not enough significant residues to build a dendrogram.")
        return None

    model = AgglomerativeClustering(n_clusters=None, metric="euclidean", linkage="single",
                                     distance_threshold=dist).fit(coords)

    # Reconstruct the linkage matrix (children_, distances_, leaf counts per
    # merge) exactly the way clustering_plot.plot_dendrogram does, since
    # scipy's dendrogram() wants a standard linkage matrix, not the raw
    # sklearn AgglomerativeClustering attributes.
    counts = np.zeros(model.children_.shape[0])
    for i, merge in enumerate(model.children_):
        c = 0
        for child_idx in merge:
            c += 1 if child_idx < n else counts[child_idx - n]
        counts[i] = c
    linkage_matrix = np.column_stack([model.children_, model.distances_, counts]).astype(float)

    def node_distance(node_id):
        if node_id < n:
            return 0.0
        idx = int(node_id - n)
        return linkage_matrix[idx, 2] if idx < len(linkage_matrix) else 0.0

    def cluster_for_node(node_id):
        if node_id < n:
            return model.labels_[node_id]
        idx = int(node_id - n)
        if idx < len(linkage_matrix):
            return cluster_for_node(int(linkage_matrix[idx, 0]))
        return -1

    palette = px.colors.qualitative.Dark24
    cluster_color_map = {}

    def color_for(cluster_id):
        cluster_id = int(cluster_id)
        if cluster_id not in cluster_color_map:
            cluster_color_map[cluster_id] = palette[len(cluster_color_map) % len(palette)]
        return cluster_color_map[cluster_id]

    def link_color_func(node_id):
        # Above the clustering threshold: gray, matching the matplotlib
        # version's above_threshold_color. Below it: colored by which
        # cluster this merge belongs to (traced down to its left leaf).
        if node_distance(node_id) > dist:
            return "#CCCCCC"
        return color_for(cluster_for_node(node_id))

    # link_color_func is called by scipy internally during the recursive
    # plotting traversal, which visits links in a different order than
    # linkage_matrix's merge order -- its returned color_list is guaranteed
    # to already be paired index-for-index with icoord/dcoord, so don't try
    # to re-derive per-link color from a naive zip against linkage_matrix.
    ddata = _scipy_dendrogram(
        linkage_matrix, no_plot=True, labels=leaf_labels,
        color_threshold=dist, link_color_func=link_color_func,
    )

    # log2 y-axis: dendrogram distances span a wide range where the small, meaningful
    # merges near the bottom get visually crushed against the near-threshold ones on a
    # linear axis. log2(0) is undefined (every leaf's own baseline is y=0, and the leaf
    # marker trace below sits at y=0 too), so floor at 1.0 Å -- any merge/leaf at or under
    # 1 Å reads as 0, never negative (a tiny data-derived epsilon floor let sub-Ångström
    # merges swing the axis far below 0, which is more precision than a structural distance
    # in Å needs to show).
    floor = 1.0

    def log2_floor(vals):
        return np.log2(np.maximum(np.asarray(vals, dtype=float), floor))

    fig = go.Figure()
    for xs, ys, color in zip(ddata["icoord"], ddata["dcoord"], ddata["color_list"]):
        fig.add_trace(go.Scatter(x=xs, y=log2_floor(ys), mode="lines", line=dict(color=color, width=1.5),
                                  hoverinfo="skip", showlegend=False))

    # scipy places leaves at x = 5, 15, 25, ... in the order given by ivl
    # (the reordered leaf labels) -- this is the standard convention for
    # converting a scipy dendrogram result into a custom-drawn plot.
    ivl = ddata["ivl"]
    leaf_x = [5 + 10 * i for i in range(len(ivl))]

    # Each leaf's OWN branch color must come from the same above/below-threshold test as
    # the lines (link_color_func), not straight from model.labels_ -- a leaf whose very
    # first (lowest) merge is already above the clustering threshold has a gray branch
    # (it never joins a below-threshold cluster with anything), but model.labels_ still
    # assigns it its own singleton cluster id, which would otherwise paint its marker a
    # "real" cluster color that doesn't match the gray line under it. Building a
    # leaf -> its-own-first-merge-distance map lets each marker use exactly the color its
    # branch would use.
    leaf_first_merge_dist = {}
    for i, merge in enumerate(model.children_):
        d = linkage_matrix[i, 2]
        for child in merge:
            if child < n and child not in leaf_first_merge_dist:
                leaf_first_merge_dist[int(child)] = d

    def leaf_marker_color(orig_idx):
        if leaf_first_merge_dist.get(orig_idx, float("inf")) > dist:
            return "#CCCCCC"
        return color_for(cluster_for_node(orig_idx))

    # ddata["leaves"] gives each plotted leaf's ORIGINAL index (same order as ivl/leaf_x),
    # which is what model.labels_ is indexed by -- add an invisible-line, visible-marker
    # trace at each leaf position so hovering shows which residue it is and which cluster
    # it landed in (the line traces above are hoverinfo="skip" since a merge line spans many
    # residues at once and has no single residue to report). Cluster id shown is still
    # model.labels_ (the real flat-clustering id) even for gray/ungrouped leaves, since
    # that's still an accurate singleton "cluster of one" label.
    leaf_clusters = [int(model.labels_[orig_idx]) for orig_idx in ddata["leaves"]]
    leaf_colors = [leaf_marker_color(orig_idx) for orig_idx in ddata["leaves"]]
    fig.add_trace(go.Scatter(
        x=leaf_x, y=[log2_floor([0])[0]] * len(leaf_x), mode="markers",
        marker=dict(color=leaf_colors, size=6),
        customdata=list(zip(ivl, leaf_clusters)),
        hovertemplate="Residue: %{customdata[0]}<br>Cluster: %{customdata[1]}<extra></extra>",
        showlegend=False,
    ))

    fig.add_hline(y=log2_floor([dist])[0], line_dash="dash", line_color="red",
                  annotation_text=f"Threshold: {dist}Å", annotation_position="top left")
    # Tick text carries residue + cluster visibly (not just on hover) -- "31-B (C15)".
    ticktext = [f"{label} (C{cluster})" for label, cluster in zip(ivl, leaf_clusters)]
    fig.update_xaxes(tickmode="array", tickvals=leaf_x, ticktext=ticktext, tickangle=90,
                      tickfont=dict(size=14))
    fig.update_yaxes(rangemode="tozero")  # floor is 1.0 Å (log2 = 0) -- never dip below 0
    fig.update_layout(title=title, yaxis_title="log2(Distance (Å))",
                       width=width, height=height, autosize=False, showlegend=False)
    return fig


def _significant_residue_coords(struc_df, sig_df, sig_col, pthr_label):
    coord_cols = ["x_coord", "y_coord", "z_coord"]
    merged = struc_df[["unipos", "unires", "chain"] + coord_cols].merge(
        sig_df[["unipos", sig_col]], on="unipos")
    merged = merged[merged[sig_col] == f"p<{pthr_label}"]
    merged = merged[~merged[coord_cols].isin(["-"]).any(axis=1)]
    if merged.empty:
        return None, None
    coords = merged[coord_cols].astype(float).to_numpy()
    leaf_labels = [f"{row.unipos}-{row.chain}" for row in merged.itertuples()]
    return coords, leaf_labels


def plot_dendrogram(output_dir, input_gene, screen_name, score_type="LFC3D",
                     direction="Positive", pthr_str="05", dist=6.0,
                     width=DEFAULT_WIDTH, height=600):
    """
    Interactive Plotly dendrogram of the spatial hierarchical clustering for
    a single screen -- a genuine merge-tree view (rotate/pan/hover), not the
    3-D cluster scatter. dist is the clustering radius in Angstroms (the
    yaml's clustering_radius, typically 6.0).
    """
    struc_path = _find_one(os.path.join(output_dir, "sequence_structure", "*_coord_struc_features.tsv"))
    nonaggr_path = _find_one(os.path.join(output_dir, score_type, f"*_{input_gene}_NonAggr_{score_type}.tsv"))
    struc_df, sig_df = _load_tsv(struc_path), _load_tsv(nonaggr_path)
    if struc_df is None or sig_df is None:
        _warn(f"Could not find structural features or NonAggr {score_type} table; skipping dendrogram.")
        return None

    sign = "pos" if direction.lower().startswith("pos") else "neg"
    sig_col = f"{screen_name}_{score_type}_{sign}_{pthr_str}_psig"
    if sig_col not in sig_df.columns:
        _warn(f"Expected column '{sig_col}' not found in {nonaggr_path}; skipping dendrogram.")
        return None

    pthr_label = _pthr_label(pthr_str)
    coords, leaf_labels = _significant_residue_coords(struc_df, sig_df, sig_col, pthr_label)
    if coords is None:
        _warn(f"No residues pass {sig_col} == p<{pthr_label}; skipping dendrogram.")
        return None

    title = f"{input_gene} {score_type} {direction} Clusters (p<{pthr_label}, {dist}Å)"
    return _agglomerative_dendrogram_figure(coords, leaf_labels, dist, title, width, height)


def plot_meta_dendrogram(output_dir, input_gene, function_for_meta="SUM", score_type="LFC3D",
                          direction="Positive", pthr_str="05", dist=6.0,
                          conservation_run=False, width=DEFAULT_WIDTH, height=600):
    """Meta-aggregate equivalent of plot_dendrogram."""
    struc_path = _find_one(os.path.join(output_dir, "sequence_structure", "*_coord_struc_features.tsv"))
    psig_path = _meta_agg_path(output_dir, input_gene, f"MetaAggr_{score_type}.tsv", conservation_run)
    struc_df, sig_df = _load_tsv(struc_path), _load_tsv(psig_path)
    if struc_df is None or sig_df is None:
        _warn(f"Could not find structural features or meta-aggregate {score_type} table; skipping dendrogram.")
        return None

    sign = "pos" if direction.lower().startswith("pos") else "neg"
    sig_col = f"{function_for_meta}_{score_type}_{sign}_{pthr_str}_psig"
    if sig_col not in sig_df.columns:
        _warn(f"Expected column '{sig_col}' not found in {psig_path}; skipping dendrogram.")
        return None

    pthr_label = _pthr_label(pthr_str)
    coords, leaf_labels = _significant_residue_coords(struc_df, sig_df, sig_col, pthr_label)
    if coords is None:
        _warn(f"No residues pass {sig_col} == p<{pthr_label}; skipping dendrogram.")
        return None

    title = f"{input_gene} {score_type} {direction} Clusters (p<{pthr_label}, {dist}Å) — Meta ({function_for_meta})"
    return _agglomerative_dendrogram_figure(coords, leaf_labels, dist, title, width, height)


# ── 5a. LFC vs LFC3D scatter ──────────────────────────────────────────────────

def _lfc_lfc3d_two_panel_figure(merged, lfc_col, lfc3d_col, title, width, height):
    """
    Shared by plot_lfc_lfc3d_scatter and plot_meta_lfc_lfc3d_scatter: a residue can have
    an LFC3D value (smoothed over structural neighbors) with no LFC value of its own (no
    direct screen data at that position) -- plotting it at some fake x position on the
    main LFC axis (the old approach) is misleading, since it has no real x coordinate at
    all. Instead it gets its own narrow left-hand panel (a vertical strip, x is just
    jitter for readability) sharing the main panel's y-axis (LFC3D), so it stays visually
    comparable without pretending it has an LFC value.
    """
    has_lfc = merged[lfc_col].notna()
    main_df = merged[has_lfc]
    no_lfc_df = merged[~has_lfc]

    fig = make_subplots(rows=1, cols=2, column_widths=[0.15, 0.85], shared_yaxes=True,
                         horizontal_spacing=0.03, subplot_titles=["No LFC", ""])

    rng = np.random.default_rng(0)  # fixed seed -- same points land at the same jittered x every call
    seen_legend_names = set()

    def _add(df_sub, col, x_col):
        for hit_type, color in _HIT_TYPE_COLORS.items():
            sub = df_sub[df_sub["hit_type"] == hit_type]
            if sub.empty:
                continue
            if x_col is None:
                x_vals = rng.uniform(-0.3, 0.3, size=len(sub))
                hovertemplate = "residue %{customdata[0]}<br>LFC3D: %{y:.3f}<extra></extra>"
            else:
                x_vals = sub[x_col]
                hovertemplate = "residue %{customdata[0]}<br>LFC: %{x:.3f}<br>LFC3D: %{y:.3f}<extra></extra>"
            fig.add_trace(go.Scatter(
                x=x_vals, y=sub[lfc3d_col], mode="markers", marker=dict(color=color, size=6),
                name=hit_type, legendgroup=hit_type, showlegend=hit_type not in seen_legend_names,
                customdata=sub[["unipos"]], hovertemplate=hovertemplate,
            ), row=1, col=col)
            seen_legend_names.add(hit_type)

    _add(no_lfc_df, 1, None)
    _add(main_df, 2, lfc_col)

    fig.update_xaxes(showticklabels=False, range=[-1, 1], row=1, col=1)
    fig.update_xaxes(title_text="LFC", row=1, col=2)
    fig.update_yaxes(title_text="LFC3D", row=1, col=1)
    fig.update_layout(title=title, width=width, height=height, autosize=False)
    return fig


def plot_lfc_lfc3d_scatter(output_dir, input_gene, screen_name, pthr_str="05",
                            width=DEFAULT_WIDTH // 2, height=DEFAULT_HEIGHT):
    # Half DEFAULT_WIDTH: unlike the other plots here, this one is normally shown alone
    # (not through show_side_by_side, since it's already a 2-column figure internally),
    # so it shouldn't stretch to the same width meant for a side-by-side pair.
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
    # Only LFC3D is required to plot at all -- LFC may legitimately be missing (see
    # _lfc_lfc3d_two_panel_figure), so don't drop those rows here.
    merged = merged.dropna(subset=[lfc3d_col])
    if pos_sig_col in psig_df.columns and neg_sig_col in psig_df.columns:
        merged = merged.merge(psig_df[["unipos", pos_sig_col, neg_sig_col]], on="unipos", how="left")

        # Compare against the actual 'p<0.05' label text, not the '05'
        # file-name suffix -- otherwise this never matches and everything
        # falls through to "Not a Hit".
        pthr_label = _pthr_label(pthr_str)

        def _label(r):
            pos_hit = str(r.get(pos_sig_col, "")).startswith(f"p<{pthr_label}")
            neg_hit = str(r.get(neg_sig_col, "")).startswith(f"p<{pthr_label}")
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

    return _lfc_lfc3d_two_panel_figure(merged, lfc_col, lfc3d_col,
                                        title=f"{input_gene} LFC vs. LFC3D",
                                        width=width, height=height)


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
    # '-' rather than NaN in these tables; coerce to numeric. bfactor_pLDDT/RSA in
    # particular must be numeric BEFORE px.scatter sees them, or else the '-' entries
    # make the whole column dtype=object, which Plotly then treats as a categorical
    # axis (sorted lexicographically, e.g. "100" before "20") instead of a numeric one.
    struc_df = struc_df.copy()
    struc_df["bfactor_pLDDT"] = pd.to_numeric(struc_df["bfactor_pLDDT"], errors="coerce")
    struc_df["RSA"] = pd.to_numeric(struc_df["RSA"], errors="coerce")

    lfc3d_df = lfc3d_df.copy()
    lfc3d_df[lfc3d_col] = pd.to_numeric(lfc3d_df[lfc3d_col], errors="coerce")
    if wght_col in lfc3d_df.columns:
        lfc3d_df[wght_col] = pd.to_numeric(lfc3d_df[wght_col], errors="coerce")

    cols = ["unipos"] + [c for c in (lfc3d_col, wght_col) if c in lfc3d_df.columns]
    merged = struc_df[["unipos"] + needed].merge(lfc3d_df[cols], on="unipos")
    merged = merged.dropna(subset=["bfactor_pLDDT", "RSA", lfc3d_col])
    merged["direction"] = np.where(merged[lfc3d_col] >= 0, "positive", "negative")
    size_col = None
    if wght_col in merged.columns:
        merged["marker_size"] = merged[wght_col].abs() * 100
        size_col = "marker_size"

    fig = px.scatter(
        merged, x="bfactor_pLDDT", y="RSA", color="direction", size=size_col,
        hover_data=["unipos", lfc3d_col],
        title=f"{input_gene} pLDDT vs. RSA",
        color_discrete_map={"positive": COLOR_POS, "negative": COLOR_NEG},
    )
    fig.update_layout(xaxis_title="pLDDT", yaxis_title="RSA",
                       width=width, height=height, autosize=False)
    return fig


# ── 5c. Hit-count bar plots by domain / pLDDT disorder category ──────────────

def _hit_counts_by_category(output_dir, input_gene, screen_name, category_df, category_col, pthr_str, label=None):
    """
    Returns None (after printing the SPECIFIC reason via _warn) rather than a generic
    "could not assemble" message, since the actual cause matters to the reader: e.g. a gene
    with no queried domain annotations at all (category column entirely NaN -- a data
    availability fact about that gene/UniProt entry, not a bug) looks nothing like a
    missing NonAggr table or zero residues clearing the significance threshold.
    """
    label = label or category_col
    nonaggr_path = _find_one(os.path.join(output_dir, "LFC3D", f"*_{input_gene}_NonAggr_LFC3D.tsv"))
    nonaggr_df = _load_tsv(nonaggr_path)
    if nonaggr_df is None:
        _warn(f"Could not find NonAggr LFC3D table (pattern: *_{input_gene}_NonAggr_LFC3D.tsv in "
              f"{output_dir}/LFC3D); skipping {label} barplot.")
        return None
    if category_col not in category_df.columns:
        _warn(f"Column '{category_col}' not found in the structural features table; skipping {label} barplot.")
        return None
    if category_df[category_col].nunique(dropna=True) == 0:
        _warn(f"No {label} annotations available for {input_gene} ('{category_col}' is empty/NaN for every "
              f"residue -- e.g. UniProt has no queried domain features for this entry) — skipping {label} barplot.")
        return None

    pos_sig_col, neg_sig_col = f"{screen_name}_LFC3D_pos_{pthr_str}_psig", f"{screen_name}_LFC3D_neg_{pthr_str}_psig"
    if pos_sig_col not in nonaggr_df.columns or neg_sig_col not in nonaggr_df.columns:
        _warn(f"Expected columns '{pos_sig_col}'/'{neg_sig_col}' not found in the NonAggr LFC3D table; "
              f"skipping {label} barplot.")
        return None

    merged = category_df[["unipos", category_col]].merge(
        nonaggr_df[["unipos", pos_sig_col, neg_sig_col]], on="unipos", how="left")
    # Compare against the actual 'p<0.05' label text, not the '05'
    # file-name suffix -- otherwise this never matches and every category
    # comes back with zero hits.
    pthr_label = _pthr_label(pthr_str)
    pos_hits = merged[merged[pos_sig_col].astype(str).str.startswith(f"p<{pthr_label}")]
    neg_hits = merged[merged[neg_sig_col].astype(str).str.startswith(f"p<{pthr_label}")]
    if pos_hits.empty and neg_hits.empty:
        _warn(f"No residues pass p<{pthr_label} (LFC3D, {screen_name}) in either direction; "
              f"skipping {label} barplot.")
        return None

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

    counts = _hit_counts_by_category(output_dir, input_gene, screen_name, domains_df, "Domain", pthr_str, label="domain")
    if counts is None:
        return None  # _hit_counts_by_category already printed the specific reason

    fig = px.bar(
        counts, x="Domain", y="count", color="hit_type", barmode="group",
        title=f"{input_gene} LFC3D hit count by domain (p<{_pthr_label(pthr_str)})",
        color_discrete_map={"POS": COLOR_POS, "NEG": COLOR_NEG},
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

    counts = _hit_counts_by_category(output_dir, input_gene, screen_name, struc_df, "pLDDT_dis", pthr_str, label="pLDDT-disorder")
    if counts is None:
        return None  # _hit_counts_by_category already printed the specific reason

    fig = px.bar(
        counts, x="pLDDT_dis", y="count", color="hit_type", barmode="group",
        title=f"{input_gene} LFC3D hit count by pLDDT-disorder category (p<{_pthr_label(pthr_str)})",
        color_discrete_map={"POS": COLOR_POS, "NEG": COLOR_NEG},
    )
    fig.update_layout(xaxis_title="pLDDT disorder category", yaxis_title="Hit count",
                       width=width, height=height, autosize=False)
    return fig


# ── 5d. Enrichment test forest plot ──────────────────────────────────────────

def plot_enrichment_test(output_dir, input_gene, screen_name=None, score_type="LFC3D", pthr_str="05",
                          feature="pLDDT_dis", width=DEFAULT_WIDTH, height=DEFAULT_HEIGHT):
    """
    Reads the enrichment_test() pickle (list of dicts with score_type,
    odds_ratio [log2], ci [log2 low/high], p_value) and renders a forest
    plot with hoverable confidence intervals.

    be3d_local.py's main() runs enrichment_test separately per feature
    ('Domain' and 'pLDDT_dis'), naming the pickle
    '{input_gene}_enrichment_test_{screen_name}_{feature}_{score_type}_{pthr_str}.pickle'
    (per-screen) or '{input_gene}_enrichment_test_{feature}_{score_type}_{pthr_str}.pickle'
    (meta-aggregate, screen_name=None) -- e.g. 'Domain' is skipped whenever a gene has no
    queried domain annotations (that feature column ends up single-valued), so feature
    defaults to 'pLDDT_dis', which is always available.
    """
    screen_part = f"{screen_name}_" if screen_name else ""
    fname = f"{input_gene}_enrichment_test_{screen_part}{feature}_{score_type}_{pthr_str}.pickle"
    path = os.path.join(output_dir, "characterization", fname)
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

    # Gray out points that don't clear the significance threshold;
    # significant points are colored by the direction of the odds ratio
    # (enriched vs. depleted), consistent with the rest of the module.
    pthr_value = float(_pthr_label(pthr_str))
    p_values = df["p_value"] if "p_value" in df.columns else pd.Series([None] * len(df))

    def _dot_color(p_value, odds_ratio):
        if p_value is None or pd.isna(p_value) or p_value >= pthr_value:
            return COLOR_NEUTRAL
        return COLOR_POS if odds_ratio >= 0 else COLOR_NEG

    marker_colors = [_dot_color(p, o) for p, o in zip(p_values, df["odds_ratio"])]

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=df["odds_ratio"], y=y_labels,
        mode="markers",
        error_x=dict(type="data", symmetric=False,
                     array=df["ci_high"] - df["odds_ratio"], arrayminus=df["odds_ratio"] - df["ci_low"]),
        marker=dict(size=10, color=marker_colors),
        text=[f"p = {p:.3g}" for p in df.get("p_value", [None] * len(df))],
        hovertemplate="%{y}: log2(OR) = %{x:.2f}<br>%{text}<extra></extra>",
    ))
    fig.add_vline(x=0, line_dash="dash", line_color="gray")
    fig.update_layout(
        title=f"{input_gene} enrichment test ({score_type}, p<{_pthr_label(pthr_str)})",
        xaxis_title="log2(odds ratio)", yaxis_title="",
        width=width, height=height, autosize=False,
    )
    return fig


# ── 6. Meta-aggregate (BE-MetaClust3D) equivalents ───────────────────────────
#
# bemetaclust3d_metaaggregate.py / bemetaclust3d_characterization.py write a
# genuinely different file/column layout than the per-screen pipeline: there's
# one shared `meta-aggregate/{input_gene}_*.tsv` per score type (no per-screen
# suffix), and the value columns are prefixed by the aggregation function name
# (function_for_meta, e.g. 'SUM', 'mean') rather than by a screen name. These
# can't reuse the per-screen functions above without silently reading the
# wrong column, so they get their own small set of equivalents instead.

def _meta_agg_path(output_dir, input_gene, suffix, conservation_run=False):
    prefix = "Merged" if conservation_run else input_gene
    return os.path.join(output_dir, "meta-aggregate", f"{prefix}_{suffix}")


def plot_meta_score_scatter(output_dir, input_gene, function_for_meta="SUM", score_type="LFC3D", direction=None,
                             pthr_str="05", conservation_run=False, width=DEFAULT_WIDTH, height=DEFAULT_HEIGHT):
    """
    Meta-aggregate equivalent of plot_score_scatter (per-residue score along
    sequence position, colored by significance). Reads
    `meta-aggregate/{input_gene}_MetaAggr_{score_type}.tsv`; value/psig
    columns are prefixed by the aggregation function name (average_split_bin_plots
    builds them as '{func}_{score_type}_pos'/'_neg' when screen_name==''),
    not by a screen name like the per-screen NonAggr table.

    direction : "positive", "negative", or None (both) -- see plot_score_scatter.
    """
    path = _meta_agg_path(output_dir, input_gene, f"MetaAggr_{score_type}.tsv", conservation_run)
    df = _load_tsv(path)
    if df is None:
        _warn(f"Could not find meta-aggregate {score_type} table ({path}); skipping plot.")
        return None

    pos_val_col = f"{function_for_meta}_{score_type}_pos"
    neg_val_col = f"{function_for_meta}_{score_type}_neg"
    pos_sig_col = f"{function_for_meta}_{score_type}_pos_{pthr_str}_psig"
    neg_sig_col = f"{function_for_meta}_{score_type}_neg_{pthr_str}_psig"
    needed = [pos_val_col, neg_val_col, pos_sig_col, neg_sig_col]
    if "unipos" not in df.columns or not all(c in df.columns for c in needed):
        _warn(f"Expected columns {['unipos'] + needed} not all found in {path}; skipping plot.")
        return None

    df = df.copy()
    df[pos_val_col] = pd.to_numeric(df[pos_val_col], errors="coerce")
    df[neg_val_col] = pd.to_numeric(df[neg_val_col], errors="coerce")

    rows = []
    for _, r in df.iterrows():
        if direction in (None, "positive") and pd.notna(r[pos_val_col]) and r[pos_val_col] != 0:
            rows.append({"unipos": r["unipos"], "value": r[pos_val_col], "direction": "positive",
                          "significant": r[pos_sig_col]})
        if direction in (None, "negative") and pd.notna(r[neg_val_col]) and r[neg_val_col] != 0:
            rows.append({"unipos": r["unipos"], "value": -abs(r[neg_val_col]), "direction": "negative",
                         "significant": r[neg_sig_col]})
    if not rows:
        _warn(f"No non-zero meta {direction or ''} {score_type} values found; skipping plot.")
        return None
    long_df = pd.DataFrame(rows)
    long_df["group"] = long_df["direction"] + " / " + long_df["significant"].astype(str)

    pthr_label = _pthr_label(pthr_str)
    title_suffix = f" ({direction})" if direction else ""
    fig = px.scatter(
        long_df, x="unipos", y="value", color="group",
        title=f"{input_gene} {score_type}{title_suffix} by sequence position — Meta ({function_for_meta})",
        color_discrete_map={
            f"positive / p<{pthr_label}": COLOR_POS, f"positive / p>={pthr_label}": COLOR_POS_PALE,
            f"negative / p<{pthr_label}": COLOR_NEG, f"negative / p>={pthr_label}": COLOR_NEG_PALE,
        },
    )
    fig.add_hline(y=0, line_color="black", line_width=1)
    fig.update_layout(xaxis_title="Residue position", yaxis_title=score_type,
                       width=width, height=height, autosize=False)
    return fig


def plot_meta_lfc_lfc3d_scatter(output_dir, input_gene, function_for_meta="SUM", pthr_str="05",
                                 conservation_run=False, width=DEFAULT_WIDTH // 2, height=DEFAULT_HEIGHT):
    """Meta-aggregate equivalent of plot_lfc_lfc3d_scatter."""
    lfc_path = _meta_agg_path(output_dir, input_gene, "LFC_dis_wght.tsv", conservation_run)
    lfc3d_path = _meta_agg_path(output_dir, input_gene, "LFC3D_dis_wght.tsv", conservation_run)
    lfc3d_psig_path = _meta_agg_path(output_dir, input_gene, "MetaAggr_LFC3D.tsv", conservation_run)
    lfc_df, lfc3d_df, psig_df = _load_tsv(lfc_path), _load_tsv(lfc3d_path), _load_tsv(lfc3d_psig_path)
    if lfc_df is None or lfc3d_df is None or psig_df is None:
        _warn("Could not find meta-aggregate LFC/LFC3D dis_wght or MetaAggr tables; skipping LFC vs LFC3D scatter.")
        return None

    lfc_col, lfc3d_col = f"{function_for_meta}_LFC", f"{function_for_meta}_LFC3D"
    pos_sig_col, neg_sig_col = f"{function_for_meta}_LFC3D_pos_{pthr_str}_psig", f"{function_for_meta}_LFC3D_neg_{pthr_str}_psig"
    if lfc_col not in lfc_df.columns or lfc3d_col not in lfc3d_df.columns:
        _warn(f"Expected columns '{lfc_col}'/'{lfc3d_col}' not found; skipping plot.")
        return None

    lfc_df, lfc3d_df = lfc_df.copy(), lfc3d_df.copy()
    lfc_df[lfc_col] = pd.to_numeric(lfc_df[lfc_col], errors="coerce")
    lfc3d_df[lfc3d_col] = pd.to_numeric(lfc3d_df[lfc3d_col], errors="coerce")

    merged = lfc_df[["unipos", lfc_col]].merge(lfc3d_df[["unipos", lfc3d_col]], on="unipos")
    # Only LFC3D is required to plot at all -- see _lfc_lfc3d_two_panel_figure.
    merged = merged.dropna(subset=[lfc3d_col])
    if pos_sig_col in psig_df.columns and neg_sig_col in psig_df.columns:
        merged = merged.merge(psig_df[["unipos", pos_sig_col, neg_sig_col]], on="unipos", how="left")
        pthr_label = _pthr_label(pthr_str)

        def _label(r):
            pos_hit = str(r.get(pos_sig_col, "")).startswith(f"p<{pthr_label}")
            neg_hit = str(r.get(neg_sig_col, "")).startswith(f"p<{pthr_label}")
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

    return _lfc_lfc3d_two_panel_figure(merged, lfc_col, lfc3d_col,
                                        title=f"{input_gene} LFC vs. LFC3D — Meta ({function_for_meta})",
                                        width=width, height=height)


def plot_meta_plddt_rsa_scatter(output_dir, input_gene, function_for_meta="SUM", score_type="LFC3D",
                                 conservation_run=False, width=DEFAULT_WIDTH, height=DEFAULT_HEIGHT):
    """Meta-aggregate equivalent of plot_plddt_rsa_scatter."""
    struc_path = _find_one(os.path.join(output_dir, "sequence_structure", "*_coord_struc_features.tsv"))
    wght_path = _meta_agg_path(output_dir, input_gene, f"{score_type}_dis_wght.tsv", conservation_run)
    struc_df, wght_df = _load_tsv(struc_path), _load_tsv(wght_path)
    if struc_df is None or wght_df is None:
        _warn("Could not find structural features or meta-aggregate dis_wght table; skipping pLDDT vs RSA scatter.")
        return None

    wght_col = f"{function_for_meta}_{score_type}_wght"
    needed = ["bfactor_pLDDT", "RSA"]
    if not all(c in struc_df.columns for c in needed) or wght_col not in wght_df.columns:
        _warn(f"Expected columns {needed + [wght_col]} not all found; skipping plot.")
        return None

    # See plot_plddt_rsa_scatter's comment: '-' placeholders make these columns
    # dtype=object unless coerced, which breaks numeric x-axis ordering in Plotly.
    struc_df = struc_df.copy()
    struc_df["bfactor_pLDDT"] = pd.to_numeric(struc_df["bfactor_pLDDT"], errors="coerce")
    struc_df["RSA"] = pd.to_numeric(struc_df["RSA"], errors="coerce")

    wght_df = wght_df.copy()
    wght_df[wght_col] = pd.to_numeric(wght_df[wght_col], errors="coerce")

    merged = struc_df[["unipos"] + needed].merge(wght_df[["unipos", wght_col]], on="unipos")
    merged = merged.dropna(subset=["bfactor_pLDDT", "RSA", wght_col])
    merged = merged[merged[wght_col] != 0]
    merged["direction"] = np.where(merged[wght_col] > 0, "positive", "negative")
    merged["marker_size"] = merged[wght_col].abs() * 100

    fig = px.scatter(
        merged, x="bfactor_pLDDT", y="RSA", color="direction", size="marker_size",
        hover_data=["unipos", wght_col],
        title=f"{input_gene} pLDDT vs. RSA — Meta ({function_for_meta})",
        color_discrete_map={"positive": COLOR_POS, "negative": COLOR_NEG},
    )
    fig.update_layout(xaxis_title="pLDDT", yaxis_title="RSA",
                       width=width, height=height, autosize=False)
    return fig


def _meta_hit_counts_by_category(output_dir, input_gene, function_for_meta, category_df, category_col,
                                  pthr_str, score_type="LFC3D", conservation_run=False, label=None):
    """Meta-aggregate equivalent of _hit_counts_by_category -- see its docstring for why
    each failure mode gets its own specific _warn() rather than one generic message."""
    label = label or category_col
    psig_path = _meta_agg_path(output_dir, input_gene, f"MetaAggr_{score_type}.tsv", conservation_run)
    psig_df = _load_tsv(psig_path)
    if psig_df is None:
        _warn(f"Could not find meta-aggregate {score_type} table ({psig_path}); skipping {label} barplot.")
        return None
    if category_col not in category_df.columns:
        _warn(f"Column '{category_col}' not found in the structural features table; skipping {label} barplot.")
        return None
    if category_df[category_col].nunique(dropna=True) == 0:
        _warn(f"No {label} annotations available for {input_gene} ('{category_col}' is empty/NaN for every "
              f"residue -- e.g. UniProt has no queried domain features for this entry) — skipping {label} barplot.")
        return None

    pos_sig_col, neg_sig_col = f"{function_for_meta}_{score_type}_pos_{pthr_str}_psig", f"{function_for_meta}_{score_type}_neg_{pthr_str}_psig"
    if pos_sig_col not in psig_df.columns or neg_sig_col not in psig_df.columns:
        _warn(f"Expected columns '{pos_sig_col}'/'{neg_sig_col}' not found in the meta-aggregate table; "
              f"skipping {label} barplot.")
        return None

    merged = category_df[["unipos", category_col]].merge(
        psig_df[["unipos", pos_sig_col, neg_sig_col]], on="unipos", how="left")
    pthr_label = _pthr_label(pthr_str)
    pos_hits = merged[merged[pos_sig_col].astype(str).str.startswith(f"p<{pthr_label}")]
    neg_hits = merged[merged[neg_sig_col].astype(str).str.startswith(f"p<{pthr_label}")]
    if pos_hits.empty and neg_hits.empty:
        _warn(f"No residues pass p<{pthr_label} (meta {score_type}) in either direction; skipping {label} barplot.")
        return None

    counts = pd.concat([
        pos_hits.groupby(category_col).size().rename("count").reset_index().assign(hit_type="POS"),
        neg_hits.groupby(category_col).size().rename("count").reset_index().assign(hit_type="NEG"),
    ], ignore_index=True)
    return counts


def plot_meta_domain_barplot(output_dir, input_gene, input_uniprot, function_for_meta="SUM", pthr_str="05",
                              conservation_run=False, width=DEFAULT_WIDTH, height=DEFAULT_HEIGHT):
    """Meta-aggregate equivalent of plot_domain_barplot."""
    domains_path = os.path.join(output_dir, "sequence_structure", f"{input_gene}_{input_uniprot}_domains.tsv")
    domains_df = _load_tsv(domains_path)
    if domains_df is None:
        _warn(f"Could not find {domains_path}; skipping domain barplot.")
        return None
    if "unipos" not in domains_df.columns and "Position" in domains_df.columns:
        domains_df = domains_df.rename(columns={"Position": "unipos"})

    counts = _meta_hit_counts_by_category(output_dir, input_gene, function_for_meta, domains_df, "Domain",
                                           pthr_str, conservation_run=conservation_run, label="domain")
    if counts is None:
        return None  # _meta_hit_counts_by_category already printed the specific reason

    fig = px.bar(
        counts, x="Domain", y="count", color="hit_type", barmode="group",
        title=f"{input_gene} LFC3D hit count by domain (p<{_pthr_label(pthr_str)}) — Meta ({function_for_meta})",
        color_discrete_map={"POS": COLOR_POS, "NEG": COLOR_NEG},
    )
    fig.update_layout(xaxis_title="Domain", yaxis_title="Hit count",
                       width=width, height=height, autosize=False)
    return fig


def plot_meta_plddt_dis_barplot(output_dir, input_gene, function_for_meta="SUM", pthr_str="05",
                                 conservation_run=False, width=DEFAULT_WIDTH, height=DEFAULT_HEIGHT):
    """Meta-aggregate equivalent of plot_plddt_dis_barplot."""
    struc_path = _find_one(os.path.join(output_dir, "sequence_structure", "*_coord_struc_features.tsv"))
    struc_df = _load_tsv(struc_path)
    if struc_df is None:
        _warn("Could not find structural features table; skipping pLDDT-disorder barplot.")
        return None

    counts = _meta_hit_counts_by_category(output_dir, input_gene, function_for_meta, struc_df, "pLDDT_dis",
                                           pthr_str, conservation_run=conservation_run, label="pLDDT-disorder")
    if counts is None:
        return None  # _meta_hit_counts_by_category already printed the specific reason

    fig = px.bar(
        counts, x="pLDDT_dis", y="count", color="hit_type", barmode="group",
        title=f"{input_gene} LFC3D hit count by pLDDT-disorder category (p<{_pthr_label(pthr_str)}) — Meta ({function_for_meta})",
        color_discrete_map={"POS": COLOR_POS, "NEG": COLOR_NEG},
    )
    fig.update_layout(xaxis_title="pLDDT disorder category", yaxis_title="Hit count",
                       width=width, height=height, autosize=False)
    return fig
