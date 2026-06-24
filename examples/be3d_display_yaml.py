"""
BEClust3D YAML Configuration Viewer
-------------------------------------
Upload a BEClust3D .yaml config file and this cell will display
all parameters with human-readable descriptions, grouped by category.

Usage in Google Colab:
  1. Run this cell to define the viewer.
  2. Call:  display_yaml_config("your_config.yaml")
"""

import yaml
from IPython.display import display, HTML

# ── Parameter descriptions ────────────────────────────────────────────────────

PARAM_DESCRIPTIONS = {
    # ── Paths & Directories ──────────────────────────────────────────────────
    "beclust3d_path": {
        "label": "BEClust3D Package Path",
        "desc": "Path to the root directory of the BEClust3D package installation.",
        "group": "Paths & Directories",
    },
    "screen_dir": {
        "label": "Screen Directory",
        "desc": "Directory where the screen input TSV/TXT files are located.",
        "group": "Paths & Directories",
    },
    "output_dir": {
        "label": "Output Directory",
        "desc": "Top-level output directory; maps to the workdir parameter in all pipeline functions.",
        "group": "Paths & Directories",
    },
    "user_pdb": {
        "label": "User PDB File",
        "desc": "Path to a user-supplied PDB file. If set, skips the AlphaFold structure query (used in sequence_structural_features).",
        "group": "Paths & Directories",
    },
    "user_fasta": {
        "label": "User FASTA File",
        "desc": "Path to a user-supplied FASTA file. If set, skips the UniProt sequence query (used in sequence_structural_features).",
        "group": "Paths & Directories",
    },
    "user_dssp": {
        "label": "User DSSP File",
        "desc": "Path to a user-supplied DSSP file. If set, skips running DSSP locally (used in sequence_structural_features).",
        "group": "Paths & Directories",
    },
    "muscle_path": {
        "label": "MUSCLE Executable Path",
        "desc": "Path to the local MUSCLE executable. Used when mode='run' in conservation.",
        "group": "Paths & Directories",
    },

    # ── Gene / Protein Identity ──────────────────────────────────────────────
    "input_gene": {
        "label": "Gene Name",
        "desc": "Gene name for this run (e.g. 'DNMT3A', 'MEN1'). Used in sequence_structural_features, parse_be_data, calculate_lfc3d, and related functions.",
        "group": "Gene / Protein Identity",
    },
    "input_uniprot": {
        "label": "UniProt Accession",
        "desc": "UniProt accession ID for input_gene. Used in sequence_structural_features and conservation.",
        "group": "Gene / Protein Identity",
    },
    "input_chain": {
        "label": "PDB Chain",
        "desc": "Chain ID of input_gene in the PDB structure. Used in sequence_structural_features and calculate_lfc3d.",
        "group": "Gene / Protein Identity",
    },

    # ── Screen Inputs ────────────────────────────────────────────────────────
    "screens": {
        "label": "Screen File(s)",
        "desc": "Comma-separated list of screen input files (TSV/TXT). Each file represents one screen and is used as the screen_name identifier throughout the pipeline.",
        "group": "Screen Inputs",
    },

    # ── Statistical Parameters ───────────────────────────────────────────────
    "nRandom": {
        "label": "Number of Randomisations",
        "desc": "Number of randomizations to perform. Used in randomize_data, randomize_sequence, calculate_lfc3d, and average_split_meta.",
        "group": "Statistical Parameters",
    },
    "pthr": {
        "label": "P-value Thresholds",
        "desc": "P-value thresholds for significance labeling. Used in znorm_score, znorm_meta, and average_split_bin_plots.",
        "group": "Statistical Parameters",
        "children": {
            "single_screen": "P-value threshold for significance labeling in single-screen analysis (used in znorm_score, average_split_bin_plots).",
            "multi_screen": "P-value threshold for significance labeling in multi-screen / meta-aggregate analysis (used in znorm_meta, average_split_bin_plots).",
        },
    },

    # ── Structure Parameters ─────────────────────────────────────────────────
    "structure_radius": {
        "label": "Structure Radius (Å)",
        "desc": "Neighbor count radius in Angstroms for structural feature calculation (used in sequence_structural_features).",
        "group": "Structure Parameters",
    },
    "clustering_radius": {
        "label": "Clustering Radius (Å)",
        "desc": "Clustering radius in Angstroms for spatial agglomerative clustering of significant residues (used in clustering).",
        "group": "Structure Parameters",
    },
    "atom_level_naa": {
        "label": "Atom-level NAA",
        "desc": "If true, counts structural neighbors at atom level rather than residue level (used in sequence_structural_features).",
        "group": "Structure Parameters",
    },

    # ── Aggregation Functions ────────────────────────────────────────────────
    "function_for_lfc": {
        "label": "LFC Aggregation Function",
        "desc": "Aggregation function applied per residue for LFC scores (used in prioritize_by_sequence and calculate_lfc3d).",
        "group": "Aggregation Functions",
    },
    "function_for_lfc3d": {
        "label": "3-D LFC Aggregation Function",
        "desc": "Aggregation function for LFC3D neighborhood scores (used in calculate_lfc3d).",
        "group": "Aggregation Functions",
    },
    "function_for_meta": {
        "label": "Meta-analysis Aggregation",
        "desc": "Aggregation function for combining scores across screens in meta-analysis (used in average_split_meta and znorm_meta).",
        "group": "Aggregation Functions",
    },

    # ── Database / Column Mappings ───────────────────────────────────────────
    "database": {
        "label": "Screen Column Mappings",
        "desc": "Maps BEClust3D's expected column names to the actual column names in your screen input files.",
        "group": "Database / Column Mappings",
        "children": {
            "mut_list_col": "Mutation list column in input DataFrames (used in parse_be_data).",
            "mut_col": "Mutation category column in input DataFrames (used in parse_be_data, plot_rawdata, hypothesis_test).",
            "val_col": "Numeric measurement / score column in input DataFrames (used in parse_be_data, plot_rawdata, hypothesis_test).",
            "gene_col": "Gene identifier column in input DataFrames (used in parse_be_data, plot_rawdata, hypothesis_test).",
            "edits_col": "Amino acid edits column in input DataFrames, e.g. 'M1V,Q2Q' (used in parse_be_data).",
            "mut_delimiter": "Delimiter used within edits_col to separate multiple edits (used in parse_be_data).",
            "gRNA_col": "Guide RNA identifier column in input DataFrames (used in parse_be_data).",
        },
    },

    # ── Mutation Categories ──────────────────────────────────────────────────
    "mutation_category": {
        "label": "Mutation Category Labels",
        "desc": "Maps BEClust3D's internal category names to the labels used in your screen files. Used in parse_be_data.",
        "group": "Mutation Categories",
        "children": {
            "missense": "Label(s) in mut_col that correspond to missense mutations.",
            "silent": "Label(s) in mut_col that correspond to silent / synonymous mutations.",
            "nonsense": "Label(s) in mut_col that correspond to nonsense (stop-gain) mutations.",
            "no_mutation": "Label(s) in mut_col that correspond to no-edit control guides.",
            "splice": "Label(s) in mut_col that correspond to splice site mutations.",
            "intron": "Label(s) in mut_col that correspond to intronic mutations.",
        },
    },

    # ── Quality Assurance ────────────────────────────────────────────────────
    "qa": {
        "label": "Quality Assurance (QA) Settings",
        "desc": "Runs Mann-Whitney U and Kolmogorov-Smirnov tests to verify that case mutation categories are significantly different from controls (used in hypothesis_test).",
        "group": "Quality Assurance",
        "children": {
            "qa_passed_only": "If true, retains only samples that pass QA for downstream analysis.",
            "qa_only": "If true, runs QA checks only and halts before the main analysis.",
            "cases": "Mutation categories treated as the case group in QA hypothesis testing.",
            "controls": "Mutation categories treated as the control group in QA hypothesis testing.",
        },
    },

    # ── Conservation Analysis ────────────────────────────────────────────────
    "conservation": {
        "label": "Conservation Analysis Settings",
        "desc": "Aligns two protein sequences and generates per-residue conservation scores (used in conservation and parse_be_data).",
        "group": "Conservation Analysis",
        "children": {
            "run": "If true, runs conservation alignment and filters; enables the conservation workflow.",
            "v_score_threshold": "Minimum conservation score to retain residues (-1, 1, 2, or 3) (used in parse_be_data).",
            "alt_gene_name": "Alternate gene name for conservation alignment, e.g. mouse ortholog (used in conservation).",
            "alt_uniprot_id": "UniProt accession ID for the alternate gene (used in conservation).",
            "alt_screen_start": "Prefix label for screens belonging to the alternate gene; used to route screens in conservation filtering.",
        },
    },

    # ── Protein-Protein Interaction ──────────────────────────────────────────
    "ppi_chain_gene_dict": {
        "label": "PPI Chain → Gene Dictionary",
        "desc": "Interacting gene to chain ID mapping for protein-protein interface analysis, e.g. {'GENE1': 'B'}. Null if not used (used in calculate_lfc3d).",
        "group": "Protein-Protein Interaction",
    },
    "ppi_gene_edits_dict": {
        "label": "PPI Gene → Edits Dictionary",
        "desc": "Interacting gene to edits dict mapping for PPI analysis. Null if not used (used in calculate_lfc3d).",
        "group": "Protein-Protein Interaction",
    },
    "priority_on_alternative": {
        "label": "Priority on Alternative Species",
        "desc": "If true, uses the alternate gene sequence as the primary reference in conservation mapping.",
        "group": "Conservation Analysis",
    },
}

# All group headers use the same gray palette
GROUP_HEADER_BG  = "#5a5a5a"
GROUP_HEADER_FG  = "#ffffff"
GROUP_ROW_BG     = "#ffffff"

# ── HTML rendering ────────────────────────────────────────────────────────────

def _val_html(v):
    """Render a YAML value as a tidy string."""
    if v is None:
        return '<span style="color:#aaa;font-style:italic;">null</span>'
    if isinstance(v, bool):
        color = "#27ae60" if v else "#e74c3c"
        return f'<span style="color:{color};font-weight:600;">{"true" if v else "false"}</span>'
    if isinstance(v, list):
        items = ", ".join(f'<code style="background:#f0f0f0;padding:1px 4px;border-radius:3px;">{x}</code>' for x in v)
        return items
    return f'<code style="background:#f0f0f0;padding:1px 4px;border-radius:3px;">{v}</code>'


def _build_html(config: dict) -> str:
    # Group params
    groups: dict[str, list] = {}
    unknown: list = []

    for key, value in config.items():
        info = PARAM_DESCRIPTIONS.get(key)
        if info:
            g = info["group"]
            groups.setdefault(g, []).append((key, value, info))
        else:
            unknown.append((key, value))

    # Collect child keys that are rendered inside parents so we don't double-render
    child_keys = set()
    for key, info in PARAM_DESCRIPTIONS.items():
        if "children" in info and key in config and isinstance(config[key], dict):
            child_keys.update(config[key].keys())

    html_parts = ["""
<style>
  .yaml-viewer { font-family: 'Segoe UI', sans-serif; max-width: 900px; margin: 0 auto; }
  .yaml-group  { margin-bottom: 14px; border-radius: 6px; overflow: hidden;
                  border: 1px solid #ddd; }
  .yaml-group-header { padding: 5px 12px; font-size: 11px; font-weight: 700;
                        letter-spacing: .07em; text-transform: uppercase; }
  .yaml-row    { display: grid; grid-template-columns: 200px 1fr 1fr;
                  gap: 0; padding: 4px 12px; border-top: 1px solid #eee;
                  background: #fff; align-items: baseline; }
  .yaml-row:hover { background: #f7f7f7; }
  .yaml-key    { font-family: monospace; font-size: 12px; color: #222; font-weight: 600; }
  .yaml-label  { font-size: 11px; color: #888; }
  .yaml-desc   { font-size: 12px; color: #555; padding: 0 8px; line-height: 1.4; }
  .yaml-val    { font-size: 12px; color: #333; line-height: 1.4; }
  .yaml-child  { padding: 3px 12px 3px 28px; border-top: 1px solid #f0f0f0;
                  background: #fafafa; display: grid;
                  grid-template-columns: 188px 1fr 1fr; gap: 0; align-items: baseline; }
  .yaml-child:hover { background: #f2f2f2; }
  .child-key   { font-family: monospace; font-size: 11px; color: #666; }
  .child-desc  { font-size: 11px; color: #888; padding: 0 8px; line-height: 1.4; }
  .child-val   { font-size: 11px; color: #333; }
  h2.yaml-title { font-family: 'Segoe UI', sans-serif; color: #222; font-size: 16px;
                   border-bottom: 2px solid #999; padding-bottom: 5px;
                   margin-bottom: 14px; }
</style>
<div class="yaml-viewer">
<h2 class="yaml-title">⚙️ BEClust3D Configuration</h2>
"""]

    for group_name, items in groups.items():
        html_parts.append(f'<div class="yaml-group">')
        html_parts.append(
            f'<div class="yaml-group-header" style="background:{GROUP_HEADER_BG};color:{GROUP_HEADER_FG};">'
            f'{group_name}</div>'
        )

        for key, value, info in items:
            label = info["label"]
            desc  = info["desc"]
            children = info.get("children", {})

            if children and isinstance(value, dict):
                # Parent row (no value cell)
                html_parts.append(
                    f'<div class="yaml-row">'
                    f'<div><div class="yaml-key">{key}</div>'
                    f'<div class="yaml-label">{label}</div></div>'
                    f'<div class="yaml-desc">{desc}</div>'
                    f'<div class="yaml-val"></div>'
                    f'</div>'
                )
                for child_key, child_desc in children.items():
                    child_val = value.get(child_key, "—")
                    html_parts.append(
                        f'<div class="yaml-child">'
                        f'<div class="child-key">↳ {child_key}</div>'
                        f'<div class="child-desc">{child_desc}</div>'
                        f'<div class="child-val">{_val_html(child_val)}</div>'
                        f'</div>'
                    )
            else:
                html_parts.append(
                    f'<div class="yaml-row">'
                    f'<div><div class="yaml-key">{key}</div>'
                    f'<div class="yaml-label">{label}</div></div>'
                    f'<div class="yaml-desc">{desc}</div>'
                    f'<div class="yaml-val">{_val_html(value)}</div>'
                    f'</div>'
                )

        html_parts.append('</div>')  # close group

    # Unknown params
    if unknown:
        html_parts.append('<div class="yaml-group">')
        html_parts.append(
            f'<div class="yaml-group-header" style="background:#888;color:#fff;">'
            'Other / Unrecognised Parameters</div>'
        )
        for key, value in unknown:
            html_parts.append(
                f'<div class="yaml-row">'
                f'<div class="yaml-key">{key}</div>'
                f'<div class="yaml-desc"><em style="color:#aaa;">No description available</em></div>'
                f'<div class="yaml-val">{_val_html(value)}</div>'
                f'</div>'
            )
        html_parts.append('</div>')

    html_parts.append('</div>')
    return "\n".join(html_parts)


# ── Main function ─────────────────────────────────────────────────────────

def display_yaml_config(yaml_path: str):
    """
    Load a BEClust3D YAML config file and display all parameters
    with descriptions in a formatted HTML table.

    Parameters
    ----------
    yaml_path : str
        Path to the .yaml config file (e.g. "dnmt3a.yaml").
    """
    with open(yaml_path, "r") as fh:
        config = yaml.safe_load(fh)

    html = _build_html(config)
    display(HTML(html))
