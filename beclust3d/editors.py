"""
File: editors.py
Author: Amit Shenoy
Date: 2026-09-03
Description:
    Editor-profile abstraction: a declarative registry of what each base /
    prime editor can install, so "is this residue reachable?" can be answered
    per-editor rather than being hard-wired to CBE (C->T) / ABE (A->G) codon
    logic.

    Base editing has fanned out well beyond the canonical transitions:
      - C->G / C->A glycosylase editors (CGBE / GBE)
      - adenine transversion editors (AYBE, A->C / A->T)
      - dual C+A editors (simultaneous C->T + A->G)
      - PAM-flexible editors (SpG / SpRY) that change *which* base is
        addressable rather than *what* substitution is installed
      - prime editing (arbitrary substitutions + small indels)
      - organelle / dsDNA editors with no PAM and no sgRNA (DdCBE, TALED)

    This module generalizes the CBE/ABE reachability report (PR #19) into a
    lookup parameterized over an editor profile. It is fully standalone,
    deterministic and offline; it does not depend on PR #19 being merged and
    does not modify the LFC3D / clustering pipeline.

    Selected citations (installable-substitution chemistry per editor):
      - CBE / ABE canonical transitions: Komor et al. 2016
        (10.1038/nature17946); Gaudelli et al. 2017 (10.1038/nature24644).
      - CGBE (C->G): Kurt et al. 2021, Nat Biotechnol 39:41-46
        (10.1038/s41587-020-0609-x).
      - GBE (C->A in E. coli, C->G in mammalian cells): Zhao et al. 2021,
        Nat Biotechnol 39:35-40 (10.1038/s41587-020-0592-2, PMID 32690970).
      - AYBE (A->C / A->T): Tong et al. 2023, Nat Biotechnol 41:1080-1084
        (10.1038/s41587-022-01595-6). ACBE / high-accuracy ACBE-Q (A->C):
        Chen et al. 2023 (10.1038/s41587-023-01821-9).
      - Dual "all-transition" (C->T + A->G): Grunewald et al. 2020 SPACE
        (10.1038/s41587-020-0535-y); Sakata et al. 2020 Target-ACEmax
        (10.1038/s41587-020-0509-0); Zhang et al. 2020 A&C-BEmax
        (10.1038/s41587-020-0527-y, PMID 32483363). TadA dual CABE-Ts:
        Lam et al. 2023 (10.1038/s41587-022-01611-9).
      - PAM flexibility (reachability multiplier): Walton et al. 2020 SpRY
        (10.1126/science.aba8853).
      - Prime editing (arbitrary substitutions + small indels): Anzalone
        et al. 2019, Nature 576:149 (10.1038/s41586-019-1711-4, PMID 31634902).
      - DdCBE (mtDNA C->T, split-DddA + paired TALEs, no guide RNA / no PAM,
        strong 5'-TC context): Mok et al. 2020, Nature 583:631
        (10.1038/s41586-020-2477-4, PMID 32641830).
      - TALED (mtDNA A->G, TALE-guided, no PAM): Cho et al. 2022, Cell
        185:1764 (10.1016/j.cell.2022.03.039, PMID 35472301).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, FrozenSet, Optional, Tuple, Union

# ----------------------------------------------------------------------------
# Standard genetic code (sense-strand codon -> single-letter amino acid, '*' =
# stop). Included in-module so amino-acid reachability is fully self-contained.
# ----------------------------------------------------------------------------
STANDARD_CODON_TABLE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# Convenience alias used as the default codon_table keyword argument.
STANDARD = STANDARD_CODON_TABLE

# Watson-Crick complement (for modeling substitutions on the antisense strand).
_COMPLEMENT: Dict[str, str] = {"A": "T", "T": "A", "C": "G", "G": "C"}

# Sentinel marking an editor that can install any substitution (prime editing).
ARBITRARY = "arbitrary"

# A single nucleotide substitution on the target (sense) strand, e.g. ('C','T').
Substitution = Tuple[str, str]


@dataclass(frozen=True)
class EditorProfile:
    """Declarative description of what a base / prime editor can install.

    Attributes:
        name: canonical editor name (e.g. 'CBE', 'ABE', 'PRIME').
        substitutions: either the sentinel ``ARBITRARY`` (prime editing, any
            substitution) or a frozenset of ``(from_base, to_base)`` nucleotide
            substitutions the editor installs on the *target (sense) strand*.
            Reachability additionally considers the reverse-complement of each
            substitution to model editing of the opposite (antisense) strand;
            see :func:`amino_acid_reachable`.
        window: short human note on the edit window (base editors) or the
            targeting/positioning constraint (prime / TALE editors).
        targeting: PAM / guide / TALE targeting note.
        notes: short metadata + source citation for the installable chemistry.
    """

    name: str
    substitutions: Union[FrozenSet[Substitution], str]
    window: str = ""
    targeting: str = ""
    notes: str = ""

    def is_arbitrary(self) -> bool:
        """True if this editor can install any substitution (prime editing)."""
        return self.substitutions == ARBITRARY


def _profile(name, subs, window, targeting, notes) -> EditorProfile:
    if subs == ARBITRARY:
        substitutions: Union[FrozenSet[Substitution], str] = ARBITRARY
    else:
        substitutions = frozenset((f.upper(), t.upper()) for f, t in subs)
    return EditorProfile(
        name=name,
        substitutions=substitutions,
        window=window,
        targeting=targeting,
        notes=notes,
    )


# ----------------------------------------------------------------------------
# Built-in editor registry. Keys are upper-cased so lookup is case-insensitive.
# ----------------------------------------------------------------------------
BUILTIN_EDITORS: Dict[str, EditorProfile] = {
    p.name.upper(): p
    for p in (
        _profile(
            "CBE",
            {("C", "T")},
            window="~4-8 nt activity window on the protospacer",
            targeting="NGG PAM (SpCas9); broadened by SpG / SpRY",
            notes="Cytosine base editor, C->T (G->A on the opposite strand). "
                  "Komor et al. 2016, Nature 533:420 (10.1038/nature17946).",
        ),
        _profile(
            "ABE",
            {("A", "G")},
            window="~4-8 nt activity window on the protospacer",
            targeting="NGG PAM (SpCas9); broadened by SpG / SpRY",
            notes="Adenine base editor, A->G (T->C on the opposite strand). "
                  "Gaudelli et al. 2017, Nature 551:464 (10.1038/nature24644); "
                  "ABE8e: Richter et al. 2020 (10.1038/s41587-020-0453-z).",
        ),
        _profile(
            "CGBE",
            {("C", "G")},
            window="~4-8 nt window; narrower purity than CBE",
            targeting="NGG PAM (SpCas9); broadened by SpG / SpRY",
            notes="C->G glycosylase base editor (APOBEC(R33A)+eUNG+nCas9). "
                  "Kurt et al. 2021, Nat Biotechnol 39:41 "
                  "(10.1038/s41587-020-0609-x). AP-site processing can create "
                  "indels (Nat Cell Biol 2023, 10.1038/s41556-023-01342-2).",
        ),
        _profile(
            "GBE",
            {("C", "A")},
            window="~4-8 nt window",
            targeting="NGG PAM (SpCas9)",
            notes="Glycosylase base editor: C->A in E. coli, C->G in "
                  "mammalian cells. Zhao et al. 2021, Nat Biotechnol 39:35 "
                  "(10.1038/s41587-020-0592-2, PMID 32690970). Registered here "
                  "with its E. coli C->A chemistry; use CGBE for mammalian C->G.",
        ),
        _profile(
            "AYBE",
            {("A", "C"), ("A", "T")},
            window="~4-8 nt window; A->Y transversions",
            targeting="NGG PAM (SpCas9)",
            notes="Adenine transversion base editor (ABE + MPG glycosylase "
                  "excising inosine), A->C / A->T. Tong et al. 2023, Nat "
                  "Biotechnol 41:1080 (10.1038/s41587-022-01595-6).",
        ),
        _profile(
            "ACBE",
            {("A", "C")},
            window="~4-8 nt window",
            targeting="NGG PAM (SpCas9)",
            notes="Adenine-to-cytosine base editor; high-accuracy ACBE-Q. "
                  "Chen et al. 2023 (10.1038/s41587-023-01821-9).",
        ),
        _profile(
            "DUAL_CBE_ABE",
            {("C", "T"), ("A", "G")},
            window="~4-8 nt window; simultaneous C->T and A->G",
            targeting="NGG PAM (SpCas9); broadened by SpG / SpRY",
            notes="Dual / 'all-transition' editor installing C->T and A->G "
                  "at once. Grunewald et al. 2020 SPACE "
                  "(10.1038/s41587-020-0535-y); Sakata et al. 2020 "
                  "Target-ACEmax (10.1038/s41587-020-0509-0); Zhang et al. "
                  "2020 A&C-BEmax (10.1038/s41587-020-0527-y, PMID 32483363); "
                  "Lam et al. 2023 TadCBE / CABE-Ts (10.1038/s41587-022-01611-9).",
        ),
        _profile(
            "PRIME",
            ARBITRARY,
            window="edit templated by the pegRNA RT template (no fixed window)",
            targeting="nCas9(H840A)-RT; PAM near the nick, broadened by SpG/SpRY",
            notes="Prime editing installs all 12 substitutions plus small "
                  "insertions and deletions, no DSB and no donor. Anzalone "
                  "et al. 2019, Nature 576:149 (10.1038/s41586-019-1711-4, "
                  "PMID 31634902). Reachability treats any AA change as "
                  "installable (substitution-level; indels are out of scope "
                  "for codon-substitution reachability).",
        ),
        _profile(
            "DDCBE",
            {("C", "T")},
            window="strong 5'-TC context; TALE-defined half-sites + spacer",
            targeting="NO PAM / NO guide RNA: paired TALE protein half-sites "
                      "(dsDNA); organelle (mtDNA).",
            notes="DddA-derived cytosine base editor for mtDNA, C->T. Split "
                  "DddA toxin + paired TALEs + UGI. Mok et al. 2020, Nature "
                  "583:631 (10.1038/s41586-020-2477-4, PMID 32641830); context "
                  "broadened by DddA6/DddA11 (Mok et al. 2022, "
                  "10.1038/s41587-022-01256-8, PMID 35379962).",
        ),
        _profile(
            "TALED",
            {("A", "G")},
            window="TALE-defined half-sites + spacer",
            targeting="NO PAM / NO guide RNA: TALE-guided dsDNA deaminase; "
                      "organelle (mtDNA).",
            notes="Transcription-activator-like effector-linked deaminase for "
                  "mitochondrial A->G editing. Cho et al. 2022, Cell 185:1764 "
                  "(10.1016/j.cell.2022.03.039, PMID 35472301); strand-selective "
                  "mitoBEs: Yi et al. 2024 (10.1038/s41587-023-01791-y).",
        ),
    )
}


def get_editor(editor: Union[str, EditorProfile]) -> EditorProfile:
    """Resolve ``editor`` (name or profile) to an :class:`EditorProfile`.

    Lookup is case-insensitive. Raises ``KeyError`` listing the known editors
    if the name is not registered.
    """
    if isinstance(editor, EditorProfile):
        return editor
    if not isinstance(editor, str):
        raise TypeError(
            f"editor must be a str name or EditorProfile, got {type(editor)!r}"
        )
    key = editor.strip().upper()
    if key not in BUILTIN_EDITORS:
        known = ", ".join(sorted(BUILTIN_EDITORS))
        raise KeyError(
            f"Unknown editor {editor!r}. Known editors: {known}."
        )
    return BUILTIN_EDITORS[key]


def installable_substitutions(
    editor: Union[str, EditorProfile]
) -> Union[FrozenSet[Substitution], str]:
    """Return the editor's installable substitution set, or ``ARBITRARY``.

    The returned set contains ``(from_base, to_base)`` substitutions on the
    target (sense) strand. For prime editing the sentinel string
    ``'arbitrary'`` (``ARBITRARY``) is returned.
    """
    return get_editor(editor).substitutions


def substitution_reachable(
    from_base: str, to_base: str, editor: Union[str, EditorProfile]
) -> bool:
    """True if ``editor`` can install the ``from_base -> to_base`` substitution.

    Only the target (sense) strand chemistry is considered here (this is a
    direct nucleotide-level query). For arbitrary editors (prime), any real
    substitution between two of A/C/G/T is reachable.
    """
    f, t = from_base.upper(), to_base.upper()
    if f not in _COMPLEMENT or t not in _COMPLEMENT:
        raise ValueError(
            f"bases must be one of A/C/G/T, got {from_base!r}, {to_base!r}"
        )
    subs = installable_substitutions(editor)
    if subs == ARBITRARY:
        return f != t
    return (f, t) in subs


def _effective_sense_changes(
    profile: EditorProfile,
) -> FrozenSet[Substitution]:
    """Sense-strand single-base changes an editor can produce.

    A deaminase acts on whichever strand carries its substrate base, so an
    editor's ``(f, t)`` substitution also produces ``(comp(f), comp(t))`` on
    the sense strand when the substrate base sits on the antisense strand
    (e.g. CBE C->T also yields G->A on the sense strand). We therefore union
    each declared substitution with its Watson-Crick complement. This is the
    documented "both strands as appropriate" assumption.
    """
    changes = set()
    for f, t in profile.substitutions:  # type: ignore[union-attr]
        changes.add((f, t))
        changes.add((_COMPLEMENT[f], _COMPLEMENT[t]))
    return frozenset(changes)


def _codons_for_aa(aa: str, codon_table: Dict[str, str]):
    aa = aa.upper()
    return [c for c, a in codon_table.items() if a == aa]


def amino_acid_reachable(
    ref_aa: str,
    alt_aa: str,
    editor: Union[str, EditorProfile],
    *,
    codon_table: Dict[str, str] = STANDARD,
) -> bool:
    """True if ``editor`` can convert amino acid ``ref_aa`` to ``alt_aa``.

    Generalizes PR #19's CBE/ABE reachability to any editor profile. For every
    standard codon encoding ``ref_aa``, apply each single-nucleotide
    substitution the editor allows (on either strand; see
    :func:`_effective_sense_changes`) at each of the three positions, translate
    the resulting codon, and return ``True`` if any yields ``alt_aa``.

    For prime / ``ARBITRARY`` editors any amino-acid change (``ref_aa`` !=
    ``alt_aa``) is reachable at the substitution level.

    Args:
        ref_aa, alt_aa: single-letter amino acids ('*' = stop is accepted).
        editor: editor name or :class:`EditorProfile`.
        codon_table: sense-strand codon -> amino-acid map (default STANDARD).
    """
    ref, alt = ref_aa.upper(), alt_aa.upper()
    profile = get_editor(editor)

    if profile.is_arbitrary():
        return ref != alt

    ref_codons = _codons_for_aa(ref, codon_table)
    if not ref_codons:
        raise ValueError(f"unknown reference amino acid {ref_aa!r}")

    changes = _effective_sense_changes(profile)
    if not changes:
        return False

    for codon in ref_codons:
        for i, base in enumerate(codon):
            for f, t in changes:
                if base == f:
                    mutated = codon[:i] + t + codon[i + 1:]
                    if codon_table.get(mutated) == alt:
                        return True
    return False


__all__ = [
    "EditorProfile",
    "BUILTIN_EDITORS",
    "STANDARD_CODON_TABLE",
    "STANDARD",
    "ARBITRARY",
    "get_editor",
    "installable_substitutions",
    "substitution_reachable",
    "amino_acid_reachable",
]
