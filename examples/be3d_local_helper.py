"""
File: be3d_local_helper.py
Description:
    Small display/plotting helpers for BE3D_local.ipynb. The pipeline (be3d_local.py)
    already writes most plots (residue dot-plots, LFC-vs-LFC3D scatter, dendrograms) as
    SVG/PNG files under each run's output_dir -- these helpers just display them, plus a
    few comparison plots (PPI vs no-PPI, blind-target) the pipeline itself doesn't produce.
"""

import os
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
}


def edit_yaml_widgets(yaml_path, editable_keys):
    """
    Small per-field ipywidgets form over a subset of a pipeline yaml config's top-level
    keys -- each edit is written straight back to yaml_path as soon as it changes, so the
    next run_be3d()/run_be3d_if_needed() call downstream picks it up. Only scalar/list
    top-level keys are exposed here; nested settings (pthr, database, mutation_category,
    conservation, qa, partners, ...) are left at the yaml's defaults -- edit the file
    directly if those need to change.
    """
    import ipywidgets as widgets

    with open(yaml_path) as f:
        config = yaml.safe_load(f)

    value_widgets = {}
    rows = []
    for key in editable_keys:
        if key not in config:
            continue
        val = config[key]
        screen_dir_val = config.get('screen_dir')
        if key == 'screens' and isinstance(val, str) and screen_dir_val and os.path.isdir(screen_dir_val):
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
        w.layout = widgets.Layout(width='550px')
        value_widgets[key] = w
        help_text = YAML_FIELD_HELP.get(key, '')
        rows.append(widgets.VBox([
            w,
            widgets.HTML(f"<div style='color:#777;font-size:12px;margin:0 0 8px 154px'>{help_text}</div>"),
        ]))

    def save(change=None):
        with open(yaml_path) as f:
            current = yaml.safe_load(f)
        for key, w in value_widgets.items():
            orig = config.get(key)
            if isinstance(w, widgets.SelectMultiple):
                v = ', '.join(w.value)
            else:
                v = w.value
                if isinstance(orig, list):
                    v = [x.strip() for x in v.split(',') if x.strip()]
                elif orig is None and v == '':
                    v = None
            current[key] = v
        with open(yaml_path, 'w') as f:
            yaml.safe_dump(current, f, sort_keys=False)

    for w in value_widgets.values():
        w.observe(save, names='value')

    display(widgets.VBox(rows))
