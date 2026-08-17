"""
BEClust3D YAML Configuration Editor
------------------------------------
Interactive ipywidgets form for creating/editing a BEClust3D .yaml config file
without hand-editing YAML text. Fields are grouped into collapsible sections;
the repeatable "screens" list (one row per input screen file) is rendered as
an add/remove table since it's the most common thing users need to change.

Usage in Google Colab / Jupyter:
    from be3d_config_editor import edit_yaml_config
    editor = edit_yaml_config("DNMT3A/DNMT3A-Lue2022.yaml")
    display(editor)
    # ... fill in fields, click "Save config" ...
    # the yaml file on disk is overwritten in place when saved
"""

import yaml
import ipywidgets as widgets
from IPython.display import display


# ── Field schema ──────────────────────────────────────────────────────────────
# Each entry describes how to render + read back one top-level (or nested)
# config key. `kind` drives which widget is built.
#   text      -> Text (string, or None if left blank and optional)
#   int       -> IntText
#   float     -> FloatText
#   bool      -> Checkbox
#   choice    -> Dropdown, needs "options"
#   list      -> Textarea, comma-separated -> list[str]

FIELD_GROUPS = [
    {
        "title": "Paths & Directories",
        "fields": [
            ("beclust3d_path", "text", "BEClust3D Package Path",
             "Path to the root directory of the BEClust3D package installation."),
            ("screen_dir", "text", "Screen Directory",
             "Directory where the screen input TSV/TXT files are located."),
            ("output_dir", "text", "Output Directory",
             "Top-level output directory; maps to workdir in all pipeline functions."),
            ("muscle_path", "text", "MUSCLE Executable Path",
             "Path to the local MUSCLE executable (used in conservation)."),
            ("user_pdb", "text", "User PDB File (optional)",
             "Path to a user-supplied PDB file. Leave blank to query AlphaFold."),
            ("user_fasta", "text", "User FASTA File (optional)",
             "Path to a user-supplied FASTA file. Leave blank to query UniProt."),
            ("user_dssp", "text", "User DSSP File (optional)",
             "Path to a user-supplied DSSP file. Leave blank to run DSSP locally."),
        ],
    },
    {
        "title": "Gene / Protein Identity",
        "fields": [
            ("input_gene", "text", "Gene Name", "e.g. DNMT3A, MEN1"),
            ("input_uniprot", "text", "UniProt Accession", "e.g. Q9Y6K1"),
            ("input_chain", "text", "PDB Chain", "Chain ID of input_gene in the PDB structure."),
        ],
    },
    {
        "title": "Screen → Column Mappings",
        "fields": [
            ("database.mut_col", "text", "Mutation category column", "e.g. Mut_type"),
            ("database.val_col", "text", "Score column", "e.g. sgRNA_score"),
            ("database.gene_col", "text", "Gene column", "e.g. Gene"),
            ("database.edits_col", "text", "Amino-acid edits column", "e.g. Mutation_list"),
            ("database.mut_delimiter", "text", "Edits delimiter", 'Delimiter within edits_col, e.g. ","'),
            ("database.mut_list_col", "text", "Mutation list column (optional)", ""),
            ("database.gRNA_col", "text", "Guide RNA column (optional)", ""),
        ],
    },
    {
        "title": "Mutation Category Labels",
        "fields": [
            ("mutation_category.missense", "list", "Missense label(s)", ""),
            ("mutation_category.silent", "list", "Silent label(s)", ""),
            ("mutation_category.nonsense", "list", "Nonsense label(s)", ""),
            ("mutation_category.no_mutation", "list", "No-mutation / control label(s)", ""),
            ("mutation_category.splice", "list", "Splice label(s)", ""),
            ("mutation_category.intron", "list", "Intron label(s)", ""),
        ],
    },
    {
        "title": "Quality Assurance",
        "fields": [
            ("qa.cases", "list", "Case categories", "Mutation categories treated as the case group."),
            ("qa.controls", "list", "Control categories", "Mutation categories treated as the control group."),
            ("qa.qa_passed_only", "bool", "QA-passed samples only", ""),
            ("qa.qa_only", "bool", "Run QA only (halt before main analysis)", ""),
        ],
    },
    {
        "title": "Statistical Parameters",
        "fields": [
            ("nRandom", "int", "Number of randomizations", ""),
            ("pthr.single_screen", "float", "P-value threshold (single screen)", ""),
            ("pthr.multi_screen", "float", "P-value threshold (multi screen)", ""),
        ],
    },
    {
        "title": "Structure Parameters",
        "fields": [
            ("structure_radius", "float", "Structure radius (Å)", ""),
            ("clustering_radius", "float", "Clustering radius (Å)", ""),
            ("atom_level_naa", "bool", "Atom-level neighbor counting", ""),
        ],
    },
    {
        "title": "Aggregation Functions",
        "fields": [
            ("function_for_lfc", "choice", "LFC aggregation", "", ["max", "mean", "median", "min", "sum"]),
            ("function_for_lfc3d", "choice", "LFC3D aggregation", "", ["max", "mean", "median", "min", "sum"]),
            ("function_for_meta", "choice", "Meta aggregation", "", ["SUM", "mean", "median", "max", "min"]),
        ],
    },
    {
        "title": "Conservation Analysis",
        "fields": [
            ("conservation.run", "bool", "Run conservation analysis", ""),
            ("conservation.v_score_threshold", "int", "Conservation score threshold", ""),
            ("conservation.alt_gene_name", "text", "Alternate gene name (optional)", "e.g. mouse ortholog"),
            ("conservation.alt_uniprot_id", "text", "Alternate UniProt accession (optional)", ""),
            ("conservation.alt_screen_start", "text", "Alternate screen prefix (optional)", ""),
            ("priority_on_alternative", "bool", "Prioritize alternate species sequence", ""),
        ],
    },
    {
        "title": "Protein-Protein Interaction (PPI)",
        "fields": [
            ("ppi_chain_gene_dict", "dict", "Chain → interacting gene",
             "One 'chain: gene_identifier' pair per line, e.g. 'C: HDAC1'. Leave blank if this isn't a PPI run."),
            ("ppi_gene_edits_dict", "dict", "Gene identifier → BEClust3D output dir",
             "One 'gene_identifier: /content/output_dir' pair per line, e.g. 'HDAC1: /content/HDAC1'."),
        ],
    },
]

# Keys handled specially outside FIELD_GROUPS: "screens" (the add/remove table).


# ── Helpers to get / set nested dotted keys ──────────────────────────────────

def _get_nested(config, dotted_key, default=None):
    keys = dotted_key.split(".")
    value = config
    for k in keys:
        if not isinstance(value, dict) or k not in value:
            return default
        value = value[k]
    return value


def _set_nested(config, dotted_key, value):
    keys = dotted_key.split(".")
    d = config
    for k in keys[:-1]:
        d = d.setdefault(k, {})
    d[keys[-1]] = value


def _list_to_text(v):
    if v is None:
        return ""
    if isinstance(v, list):
        return ", ".join(str(x) for x in v)
    return str(v)


def _text_to_list(s):
    items = [x.strip() for x in s.split(",")]
    return [x for x in items if x]


def _dict_to_text(v):
    if not v:
        return ""
    return "\n".join(f"{k}: {val}" for k, val in v.items())


def _text_to_dict(s):
    result = {}
    for line in s.splitlines():
        line = line.strip()
        if not line or ":" not in line:
            continue
        key, _, value = line.partition(":")
        key, value = key.strip(), value.strip()
        if key:
            result[key] = value
    return result or None


# ── Widget builders ───────────────────────────────────────────────────────────

def _make_widget(kind, current_value, options=None):
    layout = widgets.Layout(width="320px")
    if kind == "text":
        return widgets.Text(value="" if current_value is None else str(current_value), layout=layout)
    if kind == "int":
        return widgets.IntText(value=int(current_value) if current_value is not None else 0, layout=layout)
    if kind == "float":
        return widgets.FloatText(value=float(current_value) if current_value is not None else 0.0, layout=layout)
    if kind == "bool":
        return widgets.Checkbox(value=bool(current_value), indent=False)
    if kind == "choice":
        opts = options or []
        value = current_value if current_value in opts else (opts[0] if opts else None)
        return widgets.Dropdown(options=opts, value=value, layout=layout)
    if kind == "list":
        return widgets.Text(value=_list_to_text(current_value), layout=layout,
                             placeholder="comma-separated, e.g. Nonsense, Splice")
    if kind == "dict":
        return widgets.Textarea(value=_dict_to_text(current_value), layout=widgets.Layout(width="320px", height="60px"),
                                 placeholder="one key: value per line")
    raise ValueError(f"Unknown widget kind: {kind}")


def _read_widget(kind, widget):
    if kind == "text":
        v = widget.value.strip()
        return v if v else None
    if kind == "int":
        return int(widget.value)
    if kind == "float":
        return float(widget.value)
    if kind == "bool":
        return bool(widget.value)
    if kind == "choice":
        return widget.value
    if kind == "list":
        return _text_to_list(widget.value)
    if kind == "dict":
        return _text_to_dict(widget.value)
    raise ValueError(f"Unknown widget kind: {kind}")


# ── Screens table (add/remove rows) ──────────────────────────────────────────

class _ScreensTable:
    """Editable list of screen input filenames, one row per screen."""

    def __init__(self, initial_screens):
        self.rows_box = widgets.VBox([])
        self.add_button = widgets.Button(description="+ Add screen", button_style="info",
                                          layout=widgets.Layout(width="140px"))
        self.add_button.on_click(lambda _btn: self._add_row())
        self.widget = widgets.VBox([
            widgets.HTML("<b>Screen input files</b> "
                         "<span style='color:#888;'>(one row per screen TSV/TXT to analyze)</span>"),
            self.rows_box,
            self.add_button,
        ])
        if not initial_screens:
            initial_screens = [""]
        for s in initial_screens:
            self._add_row(s)

    def _add_row(self, value=""):
        text = widgets.Text(value=value, placeholder="path/to/screen.tsv",
                             layout=widgets.Layout(width="420px"))
        remove_btn = widgets.Button(description="✕", layout=widgets.Layout(width="32px"),
                                     button_style="danger")
        row = widgets.HBox([text, remove_btn])

        def _remove(_btn, row=row):
            rows = list(self.rows_box.children)
            if row in rows:
                rows.remove(row)
                self.rows_box.children = tuple(rows)

        remove_btn.on_click(_remove)
        self.rows_box.children = tuple(self.rows_box.children) + (row,)

    def get_values(self):
        values = []
        for row in self.rows_box.children:
            text_widget = row.children[0]
            v = text_widget.value.strip()
            if v:
                values.append(v)
        return values


# ── Main editor ───────────────────────────────────────────────────────────────

class YamlConfigEditor:
    """
    Builds an interactive ipywidgets form for a BEClust3D yaml config,
    seeded from an existing file (or a blank template if it doesn't exist yet).
    Call `.widget` (or just display the editor object) to render it, and
    read back `.config` at any time, or click "Save config" to write to disk.
    """

    def __init__(self, yaml_path):
        self.yaml_path = yaml_path
        try:
            with open(yaml_path, "r") as fh:
                self.config = yaml.safe_load(fh) or {}
        except FileNotFoundError:
            self.config = {}

        self._field_widgets = {}  # dotted_key -> (kind, widget)
        self._screens_table = _ScreensTable(self._screens_as_list())

        # Flat, always-visible section list rather than a collapsible Accordion:
        # Colab's custom widget manager has long-standing rendering gaps for
        # Accordion/Tab container widgets specifically (independent of the
        # ipywidgets version), while plain VBox/HBox/HTML render reliably there.
        section_titles = ["Screens"] + [g["title"] for g in FIELD_GROUPS]
        sections = [self._build_screens_section()]
        for group in FIELD_GROUPS:
            sections.append(self._build_group_section(group))

        self.status_label = widgets.HTML("")
        self.save_button = widgets.Button(description="\U0001F4BE Save config", button_style="success",
                                           layout=widgets.Layout(width="160px"))
        self.save_button.on_click(lambda _btn: self.save())

        section_blocks = []
        for title, section in zip(section_titles, sections):
            section_blocks.append(widgets.HTML(f"<h4 style='margin:12px 0 4px;'>{title}</h4>"))
            section_blocks.append(section)

        self.widget = widgets.VBox([
            widgets.HTML("<h3>⚙️ BEClust3D Configuration</h3>"),
            *section_blocks,
            widgets.HBox([self.save_button, self.status_label]),
        ])

    def _screens_as_list(self):
        screens = self.config.get("screens")
        if screens is None:
            return []
        if isinstance(screens, str):
            return [s.strip() for s in screens.split(",") if s.strip()]
        return list(screens)

    def _build_screens_section(self):
        return widgets.VBox([self._screens_table.widget])

    def _build_group_section(self, group):
        rows = []
        for key, kind, label, desc, *rest in group["fields"]:
            options = rest[0] if rest else None
            current = _get_nested(self.config, key)
            widget = _make_widget(kind, current, options)
            self._field_widgets[key] = (kind, widget)
            label_box = widgets.VBox([
                widgets.HTML(f"<b>{label}</b>"),
                widgets.HTML(f"<span style='color:#888;font-size:11px;'>{desc}</span>") if desc else widgets.HTML(""),
            ], layout=widgets.Layout(width="360px"))
            rows.append(widgets.HBox([label_box, widget]))
        return widgets.VBox(rows)

    def to_dict(self):
        """Read current widget values back into a plain dict, preserving unknown keys."""
        config = dict(self.config)  # keep any keys we don't render (e.g. ppi_* dicts)
        config["screens"] = self._screens_table.get_values()
        for key, (kind, widget) in self._field_widgets.items():
            _set_nested(config, key, _read_widget(kind, widget))
        return config

    def save(self):
        self.config = self.to_dict()
        with open(self.yaml_path, "w") as fh:
            yaml.dump(self.config, fh, sort_keys=False, default_flow_style=False)
        self.status_label.value = f"<span style='color:#27ae60;'>Saved to {self.yaml_path}</span>"

    def _ipython_display_(self):
        display(self.widget)


def edit_yaml_config(yaml_path):
    """
    Build and display an interactive editor for a BEClust3D yaml config.

    Parameters
    ----------
    yaml_path : str
        Path to the .yaml config file to load (or create on first Save).

    Returns
    -------
    YamlConfigEditor
        Call `.config` / `.to_dict()` to read the current values in-memory,
        or use the "Save config" button to write them to `yaml_path`.
    """
    return YamlConfigEditor(yaml_path)
