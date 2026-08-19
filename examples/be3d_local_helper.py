"""
File: be3d_local_helper.py
Description:
    Small display/plotting helpers for BE3D_local.ipynb. The pipeline (be3d_local.py)
    already writes most plots (residue dot-plots, LFC-vs-LFC3D scatter, dendrograms) as
    SVG/PNG files under each run's output_dir -- these helpers just display them, plus a
    few comparison plots (PPI vs no-PPI, blind-target) the pipeline itself doesn't produce.
"""

import os
import pandas as pd
import plotly.graph_objects as go
from IPython.display import display, Image, SVG


def show_svgs(paths):
    for path in paths:
        if os.path.exists(path):
            display(SVG(filename=path))
        else:
            print(f'[not found] {path}')


def show_images(paths):
    for path in paths:
        if os.path.exists(path):
            display(Image(filename=path))
        else:
            print(f'[not found] {path}')


def plot_residue_dot(df, value_col, title):
    """Residue position (x) vs. value (y), colored by sign, hover shows unires+unipos --
    for score tables that don't already have a pipeline-generated scatter_cutoff.svg (e.g.
    blind_target's output)."""
    plot_df = df[['unipos', 'unires', value_col]].copy()
    plot_df[value_col] = pd.to_numeric(plot_df[value_col].replace('-', pd.NA), errors='coerce').fillna(0.0)

    colors = ['#d62728' if v < 0 else '#1f77b4' if v > 0 else '#cccccc' for v in plot_df[value_col]]
    fig = go.Figure(go.Scatter(
        x=plot_df['unipos'], y=plot_df[value_col],
        mode='markers', marker=dict(color=colors, size=6, opacity=0.75),
        customdata=plot_df[['unipos', 'unires']],
        hovertemplate='residue: %{customdata[1]}%{customdata[0]}<br>value: %{y:.3f}<extra></extra>',
    ))
    fig.add_hline(y=0, line=dict(color='black', width=0.5))
    fig.update_layout(
        title=title, xaxis_title='Residue position', yaxis_title=value_col,
        plot_bgcolor='white', paper_bgcolor='white', width=900, height=320,
    )
    fig.show()


def plot_ppi_vs_noppi_scatter(df, score_label, highlight_list=None):
    """x = no-PPI score, y = PPI score, one point per residue, hover shows unires+unipos+gene;
    optionally highlights a list of unipos values (e.g. the top-N by |delta|)."""
    plot_df = df.dropna(subset=['noppi_score', 'ppi_score'])
    highlight_list = highlight_list or []
    is_highlighted = plot_df['unipos'].isin(highlight_list)

    lo = min(plot_df['noppi_score'].min(), plot_df['ppi_score'].min())
    hi = max(plot_df['noppi_score'].max(), plot_df['ppi_score'].max())
    pad = (hi - lo) * 0.05
    lo, hi = lo - pad, hi + pad

    fig = go.Figure()
    fig.add_trace(go.Scatter(
        x=[lo, hi], y=[lo, hi], mode='lines',
        line=dict(color='#9a9a94', width=1.5, dash='dash'),
        name='y = x', hoverinfo='skip', showlegend=False,
    ))
    for (chain, gene), sub in plot_df[~is_highlighted].groupby(['chain', 'gene']):
        fig.add_trace(go.Scatter(
            x=sub['noppi_score'], y=sub['ppi_score'], mode='markers',
            marker=dict(size=6, opacity=0.6),
            name=f'{gene} (chain {chain})',
            customdata=sub[['unipos', 'unires']],
            hovertemplate='residue: %{customdata[1]}%{customdata[0]}<br>'
                          f'no-PPI: %{{x:.3f}}<br>PPI: %{{y:.3f}}<extra>{gene}</extra>',
        ))
    if highlight_list:
        sub = plot_df[is_highlighted]
        fig.add_trace(go.Scatter(
            x=sub['noppi_score'], y=sub['ppi_score'], mode='markers',
            marker=dict(color='#e34948', size=9, opacity=0.9, line=dict(width=1, color='white')),
            name='highlighted', customdata=sub[['unipos', 'unires', 'gene']],
            hovertemplate='%{customdata[2]} residue: %{customdata[1]}%{customdata[0]}<br>'
                          'no-PPI: %{x:.3f}<br>PPI: %{y:.3f}<extra></extra>',
        ))

    fig.update_layout(
        title=f'no-PPI vs. PPI {score_label}',
        xaxis=dict(title=f'no-PPI {score_label}', range=[lo, hi], gridcolor='#e8e8e5'),
        yaxis=dict(title=f'PPI {score_label}', range=[lo, hi], gridcolor='#e8e8e5', scaleanchor='x', scaleratio=1),
        plot_bgcolor='white', paper_bgcolor='white', width=560, height=560,
    )
    fig.show()


def value_to_rgb(value, vmax=2.0):
    """Diverging colormap: red (negative min) -> white (0) -> blue (positive max)."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return {'r': 255, 'g': 255, 'b': 255}
    t = max(-1.0, min(1.0, float(value) / vmax))
    if t < 0:
        f = 1 + t
        return {'r': 255, 'g': round(255 * f), 'b': round(255 * f)}
    else:
        f = 1 - t
        return {'r': round(255 * f), 'g': round(255 * f), 'b': 255}


def load_molstar_pdb(widget, pdb_path):
    """Loads a local PDB file into an ipymolstar.PDBeMolstar widget (structure only, no
    coloring) -- split out from color_molstar so a single widget backing a dropdown of
    coloring options (e.g. no-PPI/PPI/delta, all the same underlying structure) only has
    to load the structure once, then just recolor on each dropdown change."""
    with open(pdb_path) as f:
        pdb_text = f.read()
    widget.custom_data = {'data': pdb_text, 'format': 'pdb', 'binary': False}
    widget.visual_style = 'cartoon'


def color_molstar(widget, chain_values, vmax=2.0, highlight_top_n=None):
    """
    Colors an already-loaded PDBeMolstar widget's residues by score via value_to_rgb --
    chain_values: {chain: {unipos: value}}. If highlight_top_n is set, the N residues with
    the largest |value| (across all chains) are additionally drawn as spheres (spacefill)
    so they stand out against the cartoon.

    NOTE: deliberately does NOT use PDBeMolstar.color() -- that helper sets color_data then
    immediately resets it to None. custom_data's structure load is asynchronous in the
    browser (see pdbemolstar.js's "change:custom_data" -> visual.update(), and the
    "change:color_data" listener that re-fires on the widget's loadComplete event); by the
    time loading finishes, color_data would already be back to None and nothing gets
    colored. Setting color_data directly and leaving it set lets the loadComplete-triggered
    re-fire pick up the real query, regardless of whether the structure was still loading.
    """
    flat = [(chain, unipos, value) for chain, pos_values in chain_values.items() for unipos, value in pos_values.items()]
    top_positions = set()
    if highlight_top_n:
        ranked = sorted(flat, key=lambda x: abs(x[2]), reverse=True)[:highlight_top_n]
        top_positions = {(chain, unipos) for chain, unipos, _ in ranked}

    query = []
    for chain, unipos, value in flat:
        entry = {
            'struct_asym_id': chain,
            'residue_number': int(unipos),
            'color': value_to_rgb(value, vmax=vmax),
        }
        if (chain, unipos) in top_positions:
            entry['representation'] = 'spacefill'
            entry['representationColor'] = value_to_rgb(value, vmax=vmax)
        query.append(entry)

    widget.color_data = {
        'data': query,
        'nonSelectedColor': {'r': 220, 'g': 220, 'b': 220},
        'keepColors': False,
        'keepRepresentations': False,
    }


def render_molstar(widget, pdb_path, chain_values, vmax=2.0, highlight_top_n=None):
    """One-shot: load a PDB into a widget and color it -- see load_molstar_pdb/color_molstar
    for the split version used by dropdown-driven single-viewer displays."""
    load_molstar_pdb(widget, pdb_path)
    color_molstar(widget, chain_values, vmax=vmax, highlight_top_n=highlight_top_n)
    return widget


def chain_values_from_df(df, value_col, chain_col='chain', pos_col='unipos'):
    values = pd.to_numeric(df[value_col].replace('-', pd.NA), errors='coerce')
    out = {}
    for chain, sub_idx in df.groupby(chain_col).groups.items():
        sub = df.loc[sub_idx]
        out[chain] = {int(p): float(v) for p, v in zip(sub[pos_col], values.loc[sub_idx]) if pd.notna(v)}
    return out
