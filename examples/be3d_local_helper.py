"""
File: be3d_local_helper.py
Description:
    Small display/plotting helpers for BE3D_local.ipynb. The pipeline (be3d_local.py)
    already writes most plots (residue dot-plots, LFC-vs-LFC3D scatter, dendrograms) as
    SVG/PNG files under each run's output_dir -- these helpers just display them, plus a
    few comparison plots (PPI vs no-PPI, blind-target) the pipeline itself doesn't produce.
"""

import os
import re
import subprocess
import sys
import time
import yaml
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
    display(fig)


def plot_ppi_vs_noppi_scatter(df, score_label, highlight_list=None, width=560, height=560):
    """Display what build_ppi_vs_noppi_scatter() builds. Use the build_ variant when the
    figure is needed as a value, e.g. as one variant behind a show_figure_dropdown()."""
    fig = build_ppi_vs_noppi_scatter(df, score_label, highlight_list=highlight_list,
                                     width=width, height=height)
    if fig is not None:
        display(fig)


def build_ppi_vs_noppi_scatter(df, score_label, highlight_list=None, width=560, height=560):
    """x = no-PPI score, y = PPI score, one point per residue, hover shows unires+unipos+gene;
    optionally highlights a list of unipos values (e.g. the top-N by |delta|). Returns the
    figure rather than displaying it."""
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
        plot_bgcolor='white', paper_bgcolor='white', width=width, height=height,
    )
    return fig


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


def _molstar_color_payload(chain_values, vmax=2.0, highlight_top_n=None):
    """The color_data payload pdbe-molstar's visual.select() expects -- same structure
    color_molstar() sets on the widget, built here so the widget-free viewer colors
    identically."""
    flat = [(chain, unipos, value)
            for chain, pos_values in chain_values.items()
            for unipos, value in pos_values.items()]
    top_positions = set()
    if highlight_top_n:
        ranked = sorted(flat, key=lambda x: abs(x[2]), reverse=True)[:highlight_top_n]
        top_positions = {(chain, unipos) for chain, unipos, _ in ranked}

    query = []
    for chain, unipos, value in flat:
        entry = {'struct_asym_id': chain, 'residue_number': int(unipos),
                 'color': value_to_rgb(value, vmax=vmax)}
        if (chain, unipos) in top_positions:
            entry['representation'] = 'spacefill'
            entry['representationColor'] = value_to_rgb(value, vmax=vmax)
        query.append(entry)

    return {'data': query,
            'nonSelectedColor': {'r': 220, 'g': 220, 'b': 220},
            'keepColors': False, 'keepRepresentations': False}


def molstar_viewer_html(pdb_path, views, vmax=2.0, highlight_top_n=None,
                         title='BE3D structure view', height=460, caption=''):
    """
    Build a self-contained HTML page showing `pdb_path` in a pdbe-molstar viewer, with the
    view selector INSIDE the viewer's own toolbar: `views` is {label: {chain: {unipos: value}}}
    and every label's color payload is embedded, so switching is a plain JS call
    (viewerInstance.visual.select) with no kernel round-trip.

    Deliberately uses no ipywidgets/anywidget layer. The widget version (PDBeMolstar +
    Dropdown.observe) needs a live kernel<->frontend comm, which is exactly what fails in
    VS Code ("Failed to load model class 'AnyModel' from module 'anywidget'") and what made
    the selection not take effect. Everything here runs in the browser: the plugin is loaded
    from a CDN and the PDB is inlined, so the page needs nothing from the kernel once
    rendered.
    """
    import base64
    import json as _json

    # Pinned to a version that exists (checked against jsDelivr's package metadata: 3.12.0 is
    # current, and the built files are pdbe-molstar-plugin.min.js / pdbe-molstar-light.min.css).
    # An earlier revision of this page pointed at a made-up 3.1.3 with the non-min filename,
    # which 404'd and left the viewer blank. Candidates are tried in order at runtime so a
    # future rename or a jsDelivr outage falls through to unpkg / the unpinned latest rather
    # than silently showing nothing.
    MOLSTAR_VERSION = '3.12.0'
    script_urls = [
        f'https://cdn.jsdelivr.net/npm/pdbe-molstar@{MOLSTAR_VERSION}/build/pdbe-molstar-plugin.min.js',
        f'https://unpkg.com/pdbe-molstar@{MOLSTAR_VERSION}/build/pdbe-molstar-plugin.min.js',
        f'https://cdn.jsdelivr.net/npm/pdbe-molstar@{MOLSTAR_VERSION}/build/pdbe-molstar-plugin.js',
        'https://cdn.jsdelivr.net/npm/pdbe-molstar/build/pdbe-molstar-plugin.min.js',
    ]

    with open(pdb_path) as f:
        pdb_text = f.read()
    pdb_b64 = base64.b64encode(pdb_text.encode()).decode()

    payloads = {label: _molstar_color_payload(cv, vmax=vmax, highlight_top_n=highlight_top_n)
                for label, cv in views.items()}
    labels = list(payloads)
    options = '\n'.join(
        f'<option value="{i}"{" selected" if i == 0 else ""}>{label}</option>'
        for i, label in enumerate(labels))

    return f"""<!DOCTYPE html>
<html><head><meta charset="utf-8">
<!-- Stylesheet is cosmetic here (the plugin's own controls are hidden), so a miss is harmless. -->
<link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/pdbe-molstar@{MOLSTAR_VERSION}/build/pdbe-molstar-light.min.css">
<style>
  html,body{{margin:0;padding:0;font-family:-apple-system,Segoe UI,Roboto,sans-serif;background:#fff}}
  .wrap{{border:1px solid #d9d9d6;border-radius:6px;overflow:hidden}}
  .bar{{display:flex;align-items:center;gap:10px;padding:7px 10px;background:#f6f6f4;
        border-bottom:1px solid #e2e2df;font-size:12.5px;color:#333;flex-wrap:wrap}}
  .bar label{{font-weight:600}}
  .bar select{{font-size:12.5px;padding:3px 6px;border:1px solid #c4c4c0;border-radius:4px;background:#fff}}
  .bar .cap{{color:#666}}
  #viewer{{position:relative;width:100%;height:{height}px}}
  #status{{padding:6px 10px;font-size:12px;color:#a33;background:#fff4f4;display:none}}
</style></head><body>
<div class="wrap">
  <div class="bar">
    <label for="viewSel">{title}</label>
    <select id="viewSel">{options}</select>
    <span class="cap">{caption}</span>
  </div>
  <div id="viewer"></div>
  <div id="status"></div>
</div>
<script>
  const PAYLOADS = {_json.dumps([payloads[l] for l in labels])};
  const statusEl = document.getElementById('status');
  function fail(msg) {{ statusEl.style.display = 'block'; statusEl.textContent = msg; }}

  // The plugin is loaded by trying a list of CDN URLs in order rather than one hardcoded
  // script tag: a single wrong path (a version that does not exist, or min/non-min naming
  // changing between releases) is otherwise indistinguishable from having no network, and
  // leaves nothing but an empty box. Each candidate is attempted until the global shows up.
  const SCRIPT_URLS = {_json.dumps(script_urls)};

  function loadPlugin(i) {{
    if (typeof PDBeMolstarPlugin !== 'undefined') {{ start(); return; }}
    if (i >= SCRIPT_URLS.length) {{
      fail('Could not load the pdbe-molstar plugin from any CDN candidate: ' + SCRIPT_URLS.join(', ')
           + ' -- check network access (a proxy or offline runtime would do this).');
      return;
    }}
    const s = document.createElement('script');
    s.src = SCRIPT_URLS[i];
    s.onload = () => {{
      if (typeof PDBeMolstarPlugin !== 'undefined') {{ start(); }}
      else {{ console.warn('loaded but no global: ' + SCRIPT_URLS[i]); loadPlugin(i + 1); }}
    }};
    s.onerror = () => {{ console.warn('failed to load: ' + SCRIPT_URLS[i]); loadPlugin(i + 1); }};
    document.head.appendChild(s);
  }}

  function start() {{
    let ready = false;
    let viewerInstance = null;

    const apply = () => {{
      if (!ready) return;
      const idx = parseInt(document.getElementById('viewSel').value, 10);
      try {{ viewerInstance.visual.select(PAYLOADS[idx]); }}
      catch (e) {{ fail('Coloring failed (visual.select): ' + e.message); console.error(e); }}
    }};

    try {{
      const pdbText = atob("{pdb_b64}");
      const blobUrl = URL.createObjectURL(new Blob([pdbText], {{type: 'text/plain'}}));
      viewerInstance = new PDBeMolstarPlugin();
      viewerInstance.render(document.getElementById('viewer'), {{
        customData: {{url: blobUrl, format: 'pdb', binary: false}},
        hideWater: true, visualStyle: 'cartoon', sequencePanel: false,
        hideControls: true, bgColor: {{r: 255, g: 255, b: 255}},
      }});

      if (viewerInstance.events && viewerInstance.events.loadComplete) {{
        viewerInstance.events.loadComplete.subscribe(() => {{ ready = true; apply(); }});
      }} else {{
        fail('Plugin loaded but events.loadComplete is missing -- the pdbe-molstar API differs from what this page expects.');
      }}
      document.getElementById('viewSel').addEventListener('change', apply);

      // If the structure never finishes loading, say so instead of leaving an empty box.
      setTimeout(() => {{
        if (!ready) fail('The structure did not finish loading within 20s -- see the browser console for the underlying error.');
      }}, 20000);
    }} catch (e) {{
      fail('Viewer setup failed: ' + e.message);
      console.error(e);
    }}
  }}

  loadPlugin(0);
</script>
</body></html>
"""


def show_molstar_viewer(pdb_path, views, vmax=2.0, highlight_top_n=None,
                         title='View', height=460, caption='', save_to=None):
    """
    Show molstar_viewer_html() inline in the notebook, in an <iframe srcdoc=...>.

    The iframe is what makes this portable: notebook frontends restrict scripts in ordinary
    HTML output (VS Code especially), but an iframe carrying its own self-contained document
    renders and runs the same way in Colab, Jupyter and VS Code -- and it needs no widget
    comm, so the in-viewer selector keeps working even where ipymolstar will not load.

    save_to : optional path to also write the same page as a standalone .html file, handy for
    sharing or for opening full-screen in a browser.
    """
    import html as _html

    from IPython.display import HTML, display as _display

    page = molstar_viewer_html(pdb_path, views, vmax=vmax, highlight_top_n=highlight_top_n,
                               title=title, height=height, caption=caption)
    if save_to:
        with open(save_to, 'w') as f:
            f.write(page)
        print(f'[html] standalone viewer: {save_to}')

    _display(HTML(
        f'<iframe srcdoc="{_html.escape(page, quote=True)}" '
        f'style="width:100%;height:{height + 90}px;border:0" '
        f'sandbox="allow-scripts allow-same-origin allow-downloads"></iframe>'
    ))


def save_molstar_html(pdb_path, chain_values, out_path, vmax=2.0, highlight_top_n=None,
                       title='BE3D structure view', height=600):
    """Single-view standalone HTML file (kept for callers that want one file per coloring);
    molstar_viewer_html/show_molstar_viewer take several views plus an in-viewer selector."""
    page = molstar_viewer_html(pdb_path, {title: chain_values}, vmax=vmax,
                                highlight_top_n=highlight_top_n, title=title, height=height)
    with open(out_path, 'w') as f:
        f.write(page)
    return out_path


def chain_values_from_df(df, value_col, chain_col='chain', pos_col='unipos'):
    values = pd.to_numeric(df[value_col].replace('-', pd.NA), errors='coerce')
    out = {}
    for chain, sub_idx in df.groupby(chain_col).groups.items():
        sub = df.loc[sub_idx]
        out[chain] = {int(p): float(v) for p, v in zip(sub[pos_col], values.loc[sub_idx]) if pd.notna(v)}
    return out


YAML_FIELD_HELP = {
    'input_gene': 'Gene symbol(s); comma-separated when there is more than one (e.g. PPI mode)',
    'input_uniprot': 'UniProt accession(s), one per gene, same order as input_gene',
    'input_chain': 'PDB chain ID(s), one per gene, same order as input_gene',
    'screen_dir': 'Directory containing the screen data file(s)',
    'screens': 'Screen data filename(s), comma-separated',
    'output_dir': 'Directory the pipeline writes its outputs to',
    'user_pdb': 'Path to a user-supplied PDB structure (blank = fetch automatically)',
    'user_fasta': 'Path to a user-supplied FASTA sequence (blank = use the UniProt sequence)',
    'user_dssp': 'Path to a user-supplied DSSP file (blank = compute automatically)',
    'nRandom': 'Number of random permutations for the null distribution (higher = slower, more precise p-values)',
    'structure_radius': 'Angstrom radius used to build the structural neighbor graph (LFC3D)',
    'clustering_radius': 'Angstrom radius used for spatial clustering / dendrograms',
    'function_for_lfc': 'Aggregation function across guides for LFC (e.g. mean)',
    'function_for_lfc3d': 'Aggregation function across structural neighbors for LFC3D (e.g. mean)',
    'function_for_meta': 'Aggregation function across screens for meta-analysis (e.g. SUM)',
    'score_type': "mode: ppi_diff merge granularity -- 'LFC3D' (per-screen) or 'Meta_LFC3D' (meta-aggregated)",
    'skip_existing': 'Skip re-running a pipeline leg if its output already exists',
    'mut_col': 'Screen column holding each guide\'s mutation type (e.g. Missense, Nonsense)',
    'val_col': 'Screen column holding each guide\'s score (e.g. sgRNA_score, LFC)',
    'gene_col': 'Screen column holding the gene symbol each guide targets',
    'edits_col': 'Screen column holding the (semicolon-joined) per-edit mutation categories',
    'mutation_category': "Raw mut_col values grouped into categories, as 'category: [raw values]'",
    'mut_categories': "Raw mut_col values grouped into categories (this partner's own list)",
    'mutation_priority': 'Most-to-least-deleterious order for collapsing a multi-edit guide to one category',
    'mut_delimiter': "Delimiter separating multiple edits within edits_col (e.g. ';')",
    'mut_list_col': 'Optional column listing individual edits separately (blank = derive from edits_col)',
    'gRNA_col': 'Optional column holding the guide/sgRNA identifier',
}


def _get_path(config, path):
    obj = config
    for part in path.split('.'):
        m = re.match(r'^(\w+)\[(\d+)\]$', part)
        obj = obj[m.group(1)][int(m.group(2))] if m else obj[part]
    return obj


def _set_path(config, path, value):
    parts = path.split('.')
    obj = config
    for part in parts[:-1]:
        m = re.match(r'^(\w+)\[(\d+)\]$', part)
        obj = obj[m.group(1)][int(m.group(2))] if m else obj[part]
    m = re.match(r'^(\w+)\[(\d+)\]$', parts[-1])
    if m:
        obj[m.group(1)][int(m.group(2))] = value
    else:
        obj[parts[-1]] = value


def _flatten_leaf_paths(obj, prefix=''):
    """
    Every leaf setting in a nested yaml config, as dotted/bracketed paths (e.g.
    'database.mut_col', 'partners[0].mut_categories') -- a dict recurses key by key, a
    list of dicts recurses index by index, anything else (scalar, or a plain list of
    scalars) is a leaf.
    """
    if isinstance(obj, dict):
        paths = []
        for k, v in obj.items():
            paths.extend(_flatten_leaf_paths(v, f'{prefix}.{k}' if prefix else k))
        return paths
    if isinstance(obj, list) and obj and isinstance(obj[0], dict):
        paths = []
        for i, item in enumerate(obj):
            paths.extend(_flatten_leaf_paths(item, f'{prefix}[{i}]'))
        return paths
    return [prefix]


def edit_yaml_widgets(yaml_path, key_groups, exclude=('mode',)):
    """
    Per-field ipywidgets form over a pipeline yaml config, grouped into labeled sections
    -- each edit is written straight back to yaml_path as soon as it changes, so the next
    run_be3d()/run_be3d_if_needed() call downstream picks it up.

    key_groups is a list of (title, keys) pairs. `title` is shown as a header above that
    group, or omitted for None (used for the top, unlabeled "important fields" group).
    `keys` is a list of dotted/bracketed paths into the config (e.g. 'database.mut_col',
    'partners[0].mut_categories'); pass None to mean "every leaf setting not already shown
    in an earlier group" -- that auto-fills a trailing catch-all group (e.g. "Advanced /
    other settings") so every field in the yaml ends up exposed somewhere. A dict-valued
    key (e.g. the whole 'mutation_category') is edited as a small YAML block rather than
    exploded field by field.
    """
    import ipywidgets as widgets

    with open(yaml_path) as f:
        config = yaml.safe_load(f)

    def _excluded(key):
        head = key.split('.')[0].split('[')[0]
        return head in exclude

    shown = set()
    resolved_groups = []
    for title, keys in key_groups:
        if keys is None:
            keys = [p for p in _flatten_leaf_paths(config)
                    if p not in shown
                    and not any(p == s or p.startswith(s + '.') or p.startswith(s + '[') for s in shown)]
        keys = [k for k in keys if not _excluded(k)]
        resolved_groups.append((title, keys))
        shown.update(keys)

    value_widgets = {}
    section_boxes = []
    for title, keys in resolved_groups:
        rows = []
        for key in keys:
            try:
                val = _get_path(config, key)
            except (KeyError, IndexError, TypeError):
                continue
            screen_dir_val = config.get('screen_dir')
            leaf_name = key.split('.')[-1].split('[')[0]
            if leaf_name == 'screens' and isinstance(val, str) and screen_dir_val and os.path.isdir(screen_dir_val):
                # Pick from the .tsv files actually present in screen_dir instead of typing
                # filenames by hand -- a typo here fails silently downstream (parse_be_data
                # just won't find the file), so a fixed option list is safer than free text.
                available = sorted(f for f in os.listdir(screen_dir_val) if f.endswith('.tsv'))
                current_screens = [s.strip() for s in val.split(',') if s.strip()]
                options = sorted(set(available) | set(current_screens))
                w = widgets.SelectMultiple(
                    value=tuple(s for s in current_screens if s in options),
                    options=options, description=key,
                    rows=min(6, max(3, len(options))),
                )
            elif isinstance(val, dict):
                # A whole nested mapping (mutation_category, conservation, ppi_chain_gene_dict,
                # ...) edited as one YAML block rather than exploded field by field.
                w = widgets.Textarea(
                    value=yaml.safe_dump(val, sort_keys=False), description=key,
                    layout=widgets.Layout(width='550px', height='120px'),
                )
            elif isinstance(val, bool):
                w = widgets.Checkbox(value=val, description=key, indent=False)
            elif isinstance(val, int):
                w = widgets.IntText(value=val, description=key)
            elif isinstance(val, float):
                w = widgets.FloatText(value=val, description=key)
            elif isinstance(val, list):
                w = widgets.Text(value=', '.join(str(v) for v in val), description=key)
            else:
                w = widgets.Text(value='' if val is None else str(val), description=key)
            w.style = {'description_width': '150px'}
            if not isinstance(w, widgets.Textarea):
                w.layout = widgets.Layout(width='550px')
            value_widgets[key] = w
            help_text = YAML_FIELD_HELP.get(leaf_name, '')
            rows.append(widgets.VBox([
                w,
                widgets.HTML(f"<div style='color:#777;font-size:12px;margin:0 0 8px 154px'>{help_text}</div>"),
            ]))
        if not rows:
            continue
        if title:
            section_boxes.append(widgets.HTML(f"<h4 style='margin:12px 0 4px'>{title}</h4>"))
        section_boxes.append(widgets.VBox(rows))

    def save(change=None):
        with open(yaml_path) as f:
            current = yaml.safe_load(f)
        for key, w in value_widgets.items():
            orig = _get_path(config, key)
            if isinstance(w, widgets.SelectMultiple):
                v = ', '.join(w.value)
            elif isinstance(w, widgets.Textarea):
                v = yaml.safe_load(w.value)
            else:
                v = w.value
                if isinstance(orig, list):
                    v = [x.strip() for x in v.split(',') if x.strip()]
                elif orig is None and v == '':
                    v = None
            _set_path(current, key, v)
        with open(yaml_path, 'w') as f:
            yaml.safe_dump(current, f, sort_keys=False)

    for w in value_widgets.values():
        w.observe(save, names='value')

    display(widgets.VBox(section_boxes))


# The 14 subfolders every be3d_local.py run creates (15 when conservation.run is set,
# adding 'conservation'); the pipeline has no per-step progress signal of its own, so
# counting how many of these have appeared under output_dir is a coarse but honest
# substitute for a bare "please wait" while a run that can take several minutes works.
BASE_OUTPUT_FOLDERS = [
    'screendata', 'screendata_rand', 'screendata_sequence', 'screendata_sequence_rand',
    'LFC', 'LFC3D', 'sequence_structure', 'cluster_LFC', 'cluster_LFC3D', 'cluster_union',
    'characterization', 'meta-aggregate', 'hypothesis_qa', 'g2p_visualization',
]


def expected_output_folders(config):
    folders = list(BASE_OUTPUT_FOLDERS)
    if (config.get('conservation') or {}).get('run'):
        folders.append('conservation')
    return folders


def pthr_fragment(value):
    """
    '0.05' -> '05', '0.001' -> '001' -- the column-name/file-name suffix the pipeline
    itself derives from a pthr value (see be3d_local.py's single_screen_pthr_str /
    multi_screen_pthr_str), reproduced here so notebook cells can pass the pthr_str a
    plot function actually needs instead of relying on its hardcoded '05' default,
    which silently reads the wrong (or a missing) column/file whenever the yaml's own
    threshold isn't 0.05.
    """
    return str(value).split('.')[1]


def run_be3d_with_progress(script, yaml_path):
    """
    Runs be3d_local.py as a subprocess and shows a live progress bar while it's
    working, tracking how many of its known output subfolders (BASE_OUTPUT_FOLDERS)
    have appeared under the yaml's output_dir so far -- run_be3d() is not "just
    BE-QA", it's the whole pipeline running (often several minutes), and a bare
    print gives no sense of whether it's stuck or a third of the way through.
    """
    import ipywidgets as widgets

    with open(yaml_path) as f:
        config = yaml.safe_load(f)
    output_dir = config['output_dir']
    folders = expected_output_folders(config)
    total = len(folders)

    bar = widgets.IntProgress(value=0, min=0, max=total, bar_style='info',
                               layout=widgets.Layout(width='400px'))
    label = widgets.HTML(value=f'starting -- 0/{total} output folders created')
    display(widgets.HBox([bar, label]))

    def count_done():
        return sum(1 for f in folders if os.path.isdir(os.path.join(output_dir, f)))

    proc = subprocess.Popen([sys.executable, script, yaml_path],
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    while proc.poll() is None:
        done = count_done()
        bar.value = done
        label.value = f'running -- {done}/{total} output folders created'
        time.sleep(2)

    output = proc.stdout.read() if proc.stdout else ''
    done = count_done()
    bar.value = done
    if proc.returncode == 0:
        bar.bar_style = 'success'
        label.value = f'done -- {done}/{total} output folders created'
    else:
        bar.bar_style = 'danger'
        label.value = f'FAILED (exit {proc.returncode}) -- {done}/{total} output folders created'
    print(output)
    if proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, [script, yaml_path], output=output)
