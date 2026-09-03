"""
File: tests/test_editor_profiles.py
Description:
    Deterministic, offline tests for the editor-profile abstraction
    (beclust3d/editors.py): the registry contents, nucleotide-substitution
    reachability, amino-acid reachability generalized over editor profiles,
    and get_editor error behavior.
"""

import pytest

from beclust3d.editors import (
    ARBITRARY,
    BUILTIN_EDITORS,
    EditorProfile,
    amino_acid_reachable,
    get_editor,
    installable_substitutions,
    substitution_reachable,
)


# ---------------------------------------------------------------------------
# (a) registry contains the expected editors with the right nucleotide sets
# ---------------------------------------------------------------------------
def test_registry_has_expected_editors():
    for name in ("CBE", "ABE", "CGBE", "AYBE", "DUAL_CBE_ABE", "PRIME",
                 "DDCBE", "TALED"):
        assert name in BUILTIN_EDITORS
        assert isinstance(BUILTIN_EDITORS[name], EditorProfile)


def test_registry_nucleotide_sets():
    assert installable_substitutions("CBE") == frozenset({("C", "T")})
    assert installable_substitutions("ABE") == frozenset({("A", "G")})
    assert installable_substitutions("CGBE") == frozenset({("C", "G")})
    assert installable_substitutions("AYBE") == frozenset({("A", "C"), ("A", "T")})
    assert installable_substitutions("DUAL_CBE_ABE") == frozenset(
        {("C", "T"), ("A", "G")}
    )
    # organelle / dsDNA editors (documented as PAM-free, TALE-guided)
    assert installable_substitutions("DDCBE") == frozenset({("C", "T")})
    assert installable_substitutions("TALED") == frozenset({("A", "G")})
    # prime editing installs anything
    assert installable_substitutions("PRIME") == ARBITRARY


def test_lookup_is_case_insensitive():
    assert get_editor("cbe") is BUILTIN_EDITORS["CBE"]
    assert get_editor("  Abe ") is BUILTIN_EDITORS["ABE"]


# ---------------------------------------------------------------------------
# (b) substitution_reachable at the nucleotide level
# ---------------------------------------------------------------------------
def test_substitution_reachable_cbe_abe():
    assert substitution_reachable("C", "T", "CBE") is True
    assert substitution_reachable("A", "G", "CBE") is False
    assert substitution_reachable("A", "G", "ABE") is True
    assert substitution_reachable("C", "T", "ABE") is False


def test_substitution_reachable_transversion_editors():
    assert substitution_reachable("C", "G", "CGBE") is True
    assert substitution_reachable("A", "C", "AYBE") is True
    assert substitution_reachable("A", "T", "AYBE") is True
    assert substitution_reachable("A", "G", "AYBE") is False  # not a transition


def test_substitution_reachable_prime_is_any_real_substitution():
    for f, t in (("A", "C"), ("G", "T"), ("C", "A"), ("T", "G")):
        assert substitution_reachable(f, t, "PRIME") is True
    assert substitution_reachable("A", "A", "PRIME") is False  # no change


def test_substitution_reachable_rejects_bad_bases():
    with pytest.raises(ValueError):
        substitution_reachable("C", "U", "CBE")  # RNA base not allowed here


# ---------------------------------------------------------------------------
# (c) amino_acid_reachable, the key generalization of PR #19
# ---------------------------------------------------------------------------
def test_prime_reaches_any_aa_change():
    # arbitrary editor: every real amino-acid change is reachable
    assert amino_acid_reachable("A", "W", "PRIME") is True
    assert amino_acid_reachable("R", "H", "PRIME") is True
    assert amino_acid_reachable("M", "V", "PRIME") is True
    # ...but a no-op change is not a "change"
    assert amino_acid_reachable("R", "R", "PRIME") is False


def test_cbe_reachable_aa_change_hand_verified():
    # HAND-VERIFIED: Arg (R) codon CGT -> His (H) codon CAT via a middle
    # G->A change. CBE deaminates a C on the *antisense* strand, which reads
    # out as G->A on the sense strand (the reverse complement of C->T).
    #   CGT (R) --[pos 2 G->A]--> CAT (H)
    # So R->H is CBE-reachable.
    assert amino_acid_reachable("R", "H", "CBE") is True


def test_cbe_cannot_do_transition_only_change():
    # HAND-VERIFIED: Met (M) has the single codon ATG; Val (V) codon GTG needs
    # a position-1 A->G transition. CBE installs only C->T (and G->A on the
    # opposite strand), never A->G, so M->V is NOT CBE-reachable.
    assert amino_acid_reachable("M", "V", "CBE") is False


def test_abe_reaches_met_to_val_but_not_the_cbe_case():
    # ABE installs A->G: ATG (M) --[pos 1 A->G]--> GTG (V) is reachable.
    assert amino_acid_reachable("M", "V", "ABE") is True
    # The R->H case above needed a G->A (=antisense C->T) change, which ABE
    # cannot install, and no other Arg codon reaches His under A->G/T->C.
    assert amino_acid_reachable("R", "H", "ABE") is False


def test_amino_acid_reachable_rejects_unknown_ref_aa():
    with pytest.raises(ValueError):
        amino_acid_reachable("Z", "H", "CBE")


# ---------------------------------------------------------------------------
# (d) get_editor unknown name raises with a listing of known editors
# ---------------------------------------------------------------------------
def test_get_editor_unknown_raises_with_listing():
    with pytest.raises(KeyError) as exc:
        get_editor("NOPE")
    msg = str(exc.value)
    assert "NOPE" in msg
    for known in ("CBE", "ABE", "PRIME"):
        assert known in msg


def test_get_editor_accepts_profile_passthrough():
    prof = BUILTIN_EDITORS["CBE"]
    assert get_editor(prof) is prof
