#!/usr/bin/env python3
"""
Formal projective covers in blocks of BGG category O.

The easiest interface is the sl2 mode, which follows the notation in the
notes directly:

    ./category_o_projectives.py sl2 2

This computes the projective cover P(2) of the irreducible highest-weight
module L(2), and prints Verma sections as M(μ).

For sl_n, use fundamental-weight coordinates:

    ./category_o_projectives.py sl3 1,0
    ./category_o_projectives.py sl4 0,1,0

For a finite Weyl group W and a regular integral block, simples, Vermas,
and indecomposable projectives are indexed by W.  This program computes the
Verma flag of the projective cover P(w) from BGG reciprocity and
Kazhdan-Lusztig polynomials.

Two common indexing conventions are supported:

    dominant:     [Delta(x) : L(w)] = P_{x,w}(1), x <= w
    antidominant: [Delta(x) : L(w)] = P_{w,x}(1), w <= x

The output is therefore a formal description of the projective cover by its
standard filtration.  It does not construct category O modules inside Sage.
"""

from __future__ import annotations

import argparse
import ast
import csv
import json
import os
import re
import sys
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path
from typing import Any, Iterable, Sequence


def _configure_dot_sage() -> None:
    """Keep Sage cache files inside this project when run via plain Python."""

    if "DOT_SAGE" in os.environ:
        return

    project_root = Path(__file__).resolve().parent
    venv_dir = project_root / ".venv"
    if venv_dir.exists():
        dot_sage = venv_dir / ".sage"
        dot_sage.mkdir(parents=True, exist_ok=True)
        os.environ["DOT_SAGE"] = str(dot_sage)


_configure_dot_sage()

try:
    from sage.all import PolynomialRing, WeylGroup, ZZ
    from sage.combinat.kazhdan_lusztig import KazhdanLusztigPolynomial
except Exception as exc:  # pragma: no cover - exercised only outside Sage.
    raise SystemExit(
        "This program needs SageMath. Run it with one of:\n"
        "  .venv/bin/python category_o_projectives.py A2 --projective e\n"
        "  DOT_SAGE=.venv/.sage .venv/bin/python category_o_projectives.py A2 --projective e\n\n"
        f"Import error: {exc}"
    ) from exc


BASE_WEIGHT_CHOICES = ("dominant", "antidominant")


@dataclass(frozen=True)
class Sl2Term:
    """One term labelled by an actual sl2 highest weight."""

    weight: Fraction
    multiplicity: int = 1

    def as_dict(self) -> dict[str, Any]:
        return {
            "weight": format_weight(self.weight),
            "multiplicity": self.multiplicity,
        }


@dataclass(frozen=True)
class Sl2ProjectiveCover:
    """Projective cover data in an sl2 block, using highest-weight labels."""

    highest_weight: Fraction
    block_kind: str
    block_weights: tuple[Fraction, ...]
    verma_flag: tuple[Sl2Term, ...]
    composition_factors: tuple[Sl2Term, ...] = ()
    structure: tuple[str, ...] = ()

    def as_dict(self) -> dict[str, Any]:
        payload = {
            "projective": f"P({format_weight(self.highest_weight)})",
            "top_simple": f"L({format_weight(self.highest_weight)})",
            "highest_weight": format_weight(self.highest_weight),
            "block_kind": self.block_kind,
            "block_weights": [format_weight(weight) for weight in self.block_weights],
            "verma_flag": [
                {
                    "verma": f"M({format_weight(term.weight)})",
                    **term.as_dict(),
                }
                for term in self.verma_flag
            ],
        }
        if self.composition_factors:
            payload["composition_factors"] = [
                {
                    "simple": f"L({format_weight(term.weight)})",
                    **term.as_dict(),
                }
                for term in self.composition_factors
            ]
        if self.structure:
            payload["structure"] = list(self.structure)
        return payload


@dataclass(frozen=True)
class SlnTerm:
    """One Verma or simple term labelled by sl_n fundamental coordinates."""

    weight: tuple[Fraction, ...]
    reduced_word: tuple[int, ...]
    length: int
    multiplicity: int
    kl_polynomial: str = "1"

    def as_verma_dict(self) -> dict[str, Any]:
        return {
            "verma": f"M({format_weight_vector(self.weight)})",
            "weight": format_weight_vector(self.weight),
            "reduced_word": list(self.reduced_word),
            "length": self.length,
            "kl_polynomial": self.kl_polynomial,
            "multiplicity": self.multiplicity,
        }

    def as_simple_dict(self) -> dict[str, Any]:
        return {
            "simple": f"L({format_weight_vector(self.weight)})",
            "weight": format_weight_vector(self.weight),
            "reduced_word": list(self.reduced_word),
            "length": self.length,
            "multiplicity": self.multiplicity,
        }


@dataclass(frozen=True)
class SlnProjectiveCover:
    """Projective cover data in a regular integral sl_n block."""

    n: int
    highest_weight: tuple[Fraction, ...]
    dominant_weight: tuple[Fraction, ...]
    block_weights: tuple[tuple[Fraction, ...], ...]
    projective_word: tuple[int, ...]
    projective_length: int
    verma_flag: tuple[SlnTerm, ...]
    composition_factors: tuple[SlnTerm, ...] = ()

    def as_dict(self) -> dict[str, Any]:
        payload = {
            "lie_algebra": f"sl{self.n}",
            "projective": f"P({format_weight_vector(self.highest_weight)})",
            "top_simple": f"L({format_weight_vector(self.highest_weight)})",
            "highest_weight": format_weight_vector(self.highest_weight),
            "dominant_block_weight": format_weight_vector(self.dominant_weight),
            "block_weights": [
                format_weight_vector(weight) for weight in self.block_weights
            ],
            "projective_reduced_word": list(self.projective_word),
            "projective_length": self.projective_length,
            "verma_flag": [term.as_verma_dict() for term in self.verma_flag],
        }
        if self.composition_factors:
            payload["composition_factors"] = [
                term.as_simple_dict() for term in self.composition_factors
            ]
        return payload


@dataclass(frozen=True)
class VermaFlagTerm:
    """One standard object appearing in a projective's Verma flag."""

    verma: str
    reduced_word: tuple[int, ...]
    length: int
    kl_polynomial: str
    multiplicity: int

    def as_dict(self) -> dict[str, Any]:
        return {
            "verma": self.verma,
            "reduced_word": list(self.reduced_word),
            "length": self.length,
            "kl_polynomial": self.kl_polynomial,
            "multiplicity": self.multiplicity,
        }


@dataclass(frozen=True)
class CompositionTerm:
    """One simple object appearing as a composition factor."""

    simple: str
    reduced_word: tuple[int, ...]
    length: int
    multiplicity: int

    def as_dict(self) -> dict[str, Any]:
        return {
            "simple": self.simple,
            "reduced_word": list(self.reduced_word),
            "length": self.length,
            "multiplicity": self.multiplicity,
        }


@dataclass(frozen=True)
class ProjectiveCover:
    """Formal data for the projective cover P(w)."""

    projective: str
    reduced_word: tuple[int, ...]
    length: int
    verma_flag: tuple[VermaFlagTerm, ...]
    composition_factors: tuple[CompositionTerm, ...] = ()

    def as_dict(self) -> dict[str, Any]:
        payload = {
            "projective": self.projective,
            "top_simple": f"L({self.projective})",
            "reduced_word": list(self.reduced_word),
            "length": self.length,
            "verma_flag": [term.as_dict() for term in self.verma_flag],
        }
        if self.composition_factors:
            payload["composition_factors"] = [
                term.as_dict() for term in self.composition_factors
            ]
        return payload


def parse_cartan_type(value: str) -> list[Any]:
    """
    Parse a compact Cartan type such as A2, G2, E8, A,2, or ['A', 2].
    """

    raw = value.strip()
    if not raw:
        raise argparse.ArgumentTypeError("Cartan type cannot be empty.")

    if raw[0] in "[(":
        try:
            parsed = ast.literal_eval(raw)
        except (SyntaxError, ValueError) as exc:
            raise argparse.ArgumentTypeError(
                f"Could not parse Cartan type literal {value!r}: {exc}"
            ) from exc
        if not isinstance(parsed, (list, tuple)) or len(parsed) < 2:
            raise argparse.ArgumentTypeError(
                "Cartan type literal should look like ['A', 2]."
            )
        return [parsed[0], *parsed[1:]]

    match = re.fullmatch(r"\s*([A-Za-z]+)\s*[, _-]?\s*(\d+)\s*", raw)
    if not match:
        raise argparse.ArgumentTypeError(
            "Use a finite Cartan type like A2, B3, G2, E8, or ['A', 2]."
        )
    letter, rank = match.groups()
    return [letter.upper(), int(rank)]


def format_cartan_type(cartan_type: Sequence[Any]) -> str:
    return "".join(str(part) for part in cartan_type)


def format_weight(weight: Fraction) -> str:
    if weight.denominator == 1:
        return str(weight.numerator)
    return f"{weight.numerator}/{weight.denominator}"


def parse_highest_weight(value: str) -> Fraction:
    try:
        return Fraction(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"Could not parse highest weight λ = {value!r}; try 0, 2, -2, or 1/2."
        ) from exc


def looks_like_highest_weight(value: str) -> bool:
    try:
        Fraction(value)
    except ValueError:
        return False
    return True


def parse_weight_vector(value: str, rank: int) -> tuple[Fraction, ...]:
    raw = value.strip()
    if raw == "0":
        return tuple(Fraction(0) for _ in range(rank))

    if raw.startswith("[") or raw.startswith("("):
        try:
            parsed = ast.literal_eval(raw)
        except (SyntaxError, ValueError) as exc:
            raise argparse.ArgumentTypeError(
                f"Could not parse weight vector {value!r}: {exc}"
            ) from exc
        if not isinstance(parsed, (list, tuple)):
            raise argparse.ArgumentTypeError(
                f"Weight vector {value!r} should be a list or tuple."
            )
        entries = [Fraction(entry) for entry in parsed]
    else:
        pieces = [piece for piece in re.split(r"[,\s]+", raw) if piece]
        entries = [Fraction(piece) for piece in pieces]

    if len(entries) != rank:
        raise argparse.ArgumentTypeError(
            f"Expected {rank} fundamental-weight coordinates for sl_{rank + 1}, "
            f"but got {len(entries)} from {value!r}."
        )
    return tuple(entries)


def format_weight_vector(weight: Sequence[Fraction]) -> str:
    return ",".join(format_weight(entry) for entry in weight)


def format_fundamental_weight_basis(rank: int) -> str:
    if rank <= 0:
        return ""
    if rank == 1:
        return "ω_1"
    if rank <= 3:
        return ", ".join(f"ω_{index}" for index in range(1, rank + 1))
    return f"ω_1, ..., ω_{rank}"


def format_sln_rho(rank: int) -> str:
    if rank <= 0:
        return "ρ = 0"
    if rank == 1:
        return "ρ = ω_1"
    if rank <= 3:
        return "ρ = " + " + ".join(f"ω_{index}" for index in range(1, rank + 1))
    return f"ρ = ω_1 + ... + ω_{rank}"


def sln_cartan_type(n: int) -> str:
    return f"A{n - 1}"


def sln_rho_epsilon(n: int) -> tuple[Fraction, ...]:
    return tuple(Fraction(n + 1 - 2 * i, 2) for i in range(1, n + 1))


def sln_dynkin_to_epsilon(weight: Sequence[Fraction]) -> tuple[Fraction, ...]:
    n = len(weight) + 1
    offset = -sum(Fraction(index + 1) * entry for index, entry in enumerate(weight)) / n
    return tuple(offset + sum(weight[index:]) for index in range(n - 1)) + (offset,)


def sln_epsilon_to_dynkin(epsilon_weight: Sequence[Fraction]) -> tuple[Fraction, ...]:
    return tuple(
        epsilon_weight[index] - epsilon_weight[index + 1]
        for index in range(len(epsilon_weight) - 1)
    )


def sln_shifted_epsilon(weight: Sequence[Fraction]) -> tuple[Fraction, ...]:
    return tuple(
        entry + rho_entry
        for entry, rho_entry in zip(sln_dynkin_to_epsilon(weight), sln_rho_epsilon(len(weight) + 1))
    )


def is_integer_fraction(value: Fraction) -> bool:
    return value.denominator == 1


def sln_is_regular_integral(weight: Sequence[Fraction]) -> bool:
    shifted = sln_shifted_epsilon(weight)
    for left_index, left in enumerate(shifted):
        for right in shifted[left_index + 1 :]:
            difference = left - right
            if difference == 0 or not is_integer_fraction(difference):
                return False
    return True


def sln_dominant_representative(
    weight: Sequence[Fraction],
) -> tuple[Fraction, ...]:
    n = len(weight) + 1
    sorted_shifted = tuple(sorted(sln_shifted_epsilon(weight), reverse=True))
    dominant_epsilon = tuple(
        entry - rho_entry
        for entry, rho_entry in zip(sorted_shifted, sln_rho_epsilon(n))
    )
    return sln_epsilon_to_dynkin(dominant_epsilon)


def sln_cartan_entry(rank: int, row: int, column: int) -> int:
    if row == column:
        return 2
    if abs(row - column) == 1:
        return -1
    return 0


def sln_coordinate_formula(simple_reflection: int, rank: int) -> str:
    variables = [f"a_{index}" for index in range(1, rank + 1)]
    reflected_index = simple_reflection - 1
    coordinates: list[str] = []

    for index, variable in enumerate(variables):
        cartan_entry = sln_cartan_entry(rank, reflected_index, index)
        if cartan_entry == 2:
            coordinates.append(f"-a_{simple_reflection} - 2")
        elif cartan_entry == -1:
            adjacent_variables = sorted((index + 1, simple_reflection))
            coordinates.append(
                f"a_{adjacent_variables[0]} + a_{adjacent_variables[1]} + 1"
            )
        else:
            coordinates.append(variable)

    return (
        f"s_{simple_reflection} · ({', '.join(variables)}) = "
        f"({', '.join(coordinates)})"
    )


def sln_coordinate_formulas(rank: int) -> tuple[str, ...]:
    return tuple(
        sln_coordinate_formula(simple_reflection, rank)
        for simple_reflection in range(1, rank + 1)
    )


def sln_apply_simple_dot(
    weight: Sequence[Fraction], simple_reflection: int
) -> tuple[Fraction, ...]:
    rank = len(weight)
    row = simple_reflection - 1
    coefficient = weight[row] + 1
    return tuple(
        entry - coefficient * sln_cartan_entry(rank, row, column)
        for column, entry in enumerate(weight)
    )


def sln_dot_action(
    dominant_weight: Sequence[Fraction], element: Any
) -> tuple[Fraction, ...]:
    weight = tuple(dominant_weight)
    for simple_reflection in reversed(element.reduced_word()):
        weight = sln_apply_simple_dot(weight, int(simple_reflection))
    return weight


def sln_find_weyl_element(
    W: Any, dominant_weight: Sequence[Fraction], target_weight: Sequence[Fraction]
) -> Any:
    target = tuple(target_weight)
    for element in sorted_weyl_elements(W):
        if sln_dot_action(dominant_weight, element) == target:
            return element
    raise ValueError(
        f"Could not find {format_weight_vector(target)} in the dot orbit of "
        f"{format_weight_vector(dominant_weight)}."
    )


def sl2_dot_partner(weight: Fraction) -> Fraction:
    """The nontrivial sl2 dot-action partner: s · λ = -λ - 2."""

    return -weight - 2


def sl2_block_weights(weight: Fraction) -> tuple[Fraction, ...]:
    if weight.denominator != 1 or weight == -1:
        return (weight,)

    partner = sl2_dot_partner(weight)
    if weight >= 0:
        return (weight, partner)
    return (partner, weight)


def sl2_block_kind(weight: Fraction) -> str:
    if weight.denominator != 1:
        return "nonintegral singleton block"
    if weight == -1:
        return "singular integral singleton block"
    return "regular integral block"


def sl2_projective_cover(
    highest_weight: Fraction, include_composition: bool = False
) -> Sl2ProjectiveCover:
    block_weights = sl2_block_weights(highest_weight)
    block_kind = sl2_block_kind(highest_weight)

    if len(block_weights) == 1:
        label = format_weight(highest_weight)
        return Sl2ProjectiveCover(
            highest_weight=highest_weight,
            block_kind=block_kind,
            block_weights=block_weights,
            verma_flag=(Sl2Term(highest_weight),),
            composition_factors=(Sl2Term(highest_weight),)
            if include_composition
            else (),
            structure=(f"P({label}) = M({label}) = L({label}).",),
        )

    dominant_weight, antidominant_weight = block_weights
    if highest_weight == dominant_weight:
        dominant_label = format_weight(dominant_weight)
        antidominant_label = format_weight(antidominant_weight)
        return Sl2ProjectiveCover(
            highest_weight=highest_weight,
            block_kind=block_kind,
            block_weights=block_weights,
            verma_flag=(Sl2Term(dominant_weight),),
            composition_factors=(
                Sl2Term(dominant_weight),
                Sl2Term(antidominant_weight),
            )
            if include_composition
            else (),
            structure=(
                f"P({dominant_label}) = M({dominant_label}).",
                "0 -> "
                f"M({antidominant_label}) = L({antidominant_label}) -> "
                f"M({dominant_label}) -> L({dominant_label}) -> 0.",
            ),
        )

    dominant_label = format_weight(dominant_weight)
    antidominant_label = format_weight(antidominant_weight)
    return Sl2ProjectiveCover(
        highest_weight=highest_weight,
        block_kind=block_kind,
        block_weights=block_weights,
        verma_flag=(Sl2Term(antidominant_weight), Sl2Term(dominant_weight)),
        composition_factors=(
            Sl2Term(antidominant_weight, 2),
            Sl2Term(dominant_weight),
        )
        if include_composition
        else (),
        structure=(
            "0 -> "
            f"P({dominant_label}) = M({dominant_label}) -> "
            f"P({antidominant_label}) -> "
            f"M({antidominant_label}) = L({antidominant_label}) -> 0.",
        ),
    )


def sl2_projective_covers(
    highest_weight: Fraction,
    include_composition: bool = False,
    all_in_block: bool = False,
) -> tuple[Sl2ProjectiveCover, ...]:
    weights = sl2_block_weights(highest_weight) if all_in_block else (highest_weight,)
    return tuple(
        sl2_projective_cover(weight, include_composition=include_composition)
        for weight in weights
    )


def sort_terms_for_projective(
    terms: Iterable[SlnTerm], projective_word: tuple[int, ...]
) -> tuple[SlnTerm, ...]:
    return tuple(
        sorted(
            terms,
            key=lambda term: (
                term.reduced_word != projective_word,
                -term.length,
                term.reduced_word,
            ),
        )
    )


def sln_projective_cover(
    n: int,
    W: Any,
    KL: Any,
    dominant_weight: tuple[Fraction, ...],
    block_weights: tuple[tuple[Fraction, ...], ...],
    projective_element: Any,
    include_composition: bool = False,
) -> SlnProjectiveCover:
    projective_word = word_tuple(projective_element)
    verma_terms: list[SlnTerm] = []
    for verma_element in sorted_weyl_elements(W):
        polynomial, multiplicity = standard_multiplicity(
            KL, verma_element, projective_element, "dominant"
        )
        if multiplicity == 0:
            continue
        verma_terms.append(
            SlnTerm(
                weight=sln_dot_action(dominant_weight, verma_element),
                reduced_word=word_tuple(verma_element),
                length=int(verma_element.length()),
                kl_polynomial=polynomial,
                multiplicity=multiplicity,
            )
        )

    simple_terms: list[SlnTerm] = []
    if include_composition:
        multiplicities: dict[Any, int] = {element: 0 for element in W}
        for verma_element in sorted_weyl_elements(W):
            _, projective_standard_multiplicity = standard_multiplicity(
                KL, verma_element, projective_element, "dominant"
            )
            if projective_standard_multiplicity == 0:
                continue
            for simple_element in sorted_weyl_elements(W):
                _, verma_composition_multiplicity = standard_multiplicity(
                    KL, verma_element, simple_element, "dominant"
                )
                multiplicities[simple_element] += (
                    projective_standard_multiplicity
                    * verma_composition_multiplicity
                )
        for simple_element, multiplicity in multiplicities.items():
            if multiplicity == 0:
                continue
            simple_terms.append(
                SlnTerm(
                    weight=sln_dot_action(dominant_weight, simple_element),
                    reduced_word=word_tuple(simple_element),
                    length=int(simple_element.length()),
                    multiplicity=multiplicity,
                )
            )

    return SlnProjectiveCover(
        n=n,
        highest_weight=sln_dot_action(dominant_weight, projective_element),
        dominant_weight=dominant_weight,
        block_weights=block_weights,
        projective_word=projective_word,
        projective_length=int(projective_element.length()),
        verma_flag=sort_terms_for_projective(verma_terms, projective_word),
        composition_factors=sort_terms_for_projective(simple_terms, projective_word),
    )


def sln_projective_covers(
    n: int,
    highest_weight: tuple[Fraction, ...],
    include_composition: bool = False,
    all_in_block: bool = False,
) -> tuple[SlnProjectiveCover, ...]:
    if n < 2:
        raise ValueError("sl_n mode requires n >= 2.")
    if len(highest_weight) != n - 1:
        raise ValueError(
            f"sl_{n} weights need {n - 1} fundamental-weight coordinates."
        )
    if not sln_is_regular_integral(highest_weight):
        raise ValueError(
            "sln mode currently supports regular integral blocks only. "
            "In fundamental-weight coordinates, this includes dominant "
            "integral weights such as 0, 1,0, or 0,1,0. Singular and "
            "nonintegral blocks need parabolic/integral-subgroup KL data."
        )

    W = WeylGroup(["A", n - 1], prefix="s")
    KL = make_kl(W)
    dominant_weight = sln_dominant_representative(highest_weight)
    block_elements = sorted_weyl_elements(W)
    block_weights = tuple(
        sln_dot_action(dominant_weight, element) for element in block_elements
    )
    projective_element = sln_find_weyl_element(W, dominant_weight, highest_weight)
    elements = block_elements if all_in_block else [projective_element]
    return tuple(
        sln_projective_cover(
            n,
            W,
            KL,
            dominant_weight,
            block_weights,
            element,
            include_composition=include_composition,
        )
        for element in elements
    )


def parse_word(text: str, W: Any) -> Any:
    """
    Parse a Weyl group element from a reduced word.

    Examples accepted by the CLI:
      e, id, identity, []       identity element
      w0, longest              longest element
      s1, 1, 1,2,1, s1*s2*s1  products of simple reflections
    """

    raw = text.strip()
    lowered = raw.lower()
    if lowered in {"e", "id", "identity", "[]", "()"}:
        return W.one()
    if lowered in {"w0", "longest"}:
        return W.long_element()

    indices = [int(token) for token in re.findall(r"\d+", raw)]
    if not indices:
        raise argparse.ArgumentTypeError(
            f"Could not parse Weyl group word {text!r}; try e, w0, or 1,2,1."
        )

    try:
        element = W.from_reduced_word(indices)
    except Exception as exc:
        raise argparse.ArgumentTypeError(
            f"Invalid reduced word {text!r} for this Weyl group: {exc}"
        ) from exc
    if int(element.length()) != len(indices):
        raise argparse.ArgumentTypeError(
            f"The word {text!r} is not reduced; use {word_label(element)!r} instead."
        )
    return element


def word_tuple(element: Any) -> tuple[int, ...]:
    return tuple(int(i) for i in element.reduced_word())


def word_label(element: Any) -> str:
    word = word_tuple(element)
    if not word:
        return "e"
    return "*".join(f"s{i}" for i in word)


def sorted_weyl_elements(W: Any) -> list[Any]:
    return sorted(list(W), key=lambda element: (element.length(), word_tuple(element)))


def make_kl(W: Any) -> Any:
    q_ring = PolynomialRing(ZZ, "q")
    return KazhdanLusztigPolynomial(W, q_ring.gen())


def kl_at_one(KL: Any, lower: Any, upper: Any) -> tuple[str, int]:
    polynomial = KL.P(lower, upper)
    return str(polynomial), int(polynomial(1))


def standard_multiplicity(
    KL: Any, verma_element: Any, simple_element: Any, base_weight: str
) -> tuple[str, int]:
    """Return [Delta(verma_element) : L(simple_element)]."""

    if base_weight == "dominant":
        if not verma_element.bruhat_le(simple_element):
            return "0", 0
        return kl_at_one(KL, verma_element, simple_element)

    if base_weight == "antidominant":
        if not simple_element.bruhat_le(verma_element):
            return "0", 0
        return kl_at_one(KL, simple_element, verma_element)

    raise ValueError(f"Unknown base weight convention: {base_weight!r}")


def convention_text(base_weight: str) -> str:
    if base_weight == "dominant":
        return "[P(w): Delta(x)] = [Delta(x): L(w)] = P_{x,w}(1), x <= w"
    if base_weight == "antidominant":
        return "[P(w): Delta(x)] = [Delta(x): L(w)] = P_{w,x}(1), w <= x"
    raise ValueError(f"Unknown base weight convention: {base_weight!r}")


def kl_header(base_weight: str) -> str:
    return "P_{x,w}(q)" if base_weight == "dominant" else "P_{w,x}(q)"


def verma_flag(
    W: Any, KL: Any, projective_element: Any, base_weight: str
) -> tuple[VermaFlagTerm, ...]:
    terms: list[VermaFlagTerm] = []
    for verma_element in sorted_weyl_elements(W):
        polynomial, multiplicity = standard_multiplicity(
            KL, verma_element, projective_element, base_weight
        )
        if multiplicity == 0:
            continue
        terms.append(
            VermaFlagTerm(
                verma=word_label(verma_element),
                reduced_word=word_tuple(verma_element),
                length=int(verma_element.length()),
                kl_polynomial=polynomial,
                multiplicity=multiplicity,
            )
        )
    return tuple(terms)


def composition_factors(
    W: Any, KL: Any, projective_element: Any, base_weight: str
) -> tuple[CompositionTerm, ...]:
    """
    Compute composition multiplicities from the Verma flag.

    This reports total composition multiplicities, not Loewy layers.
    """

    multiplicities: dict[Any, int] = {element: 0 for element in W}
    for verma_element in sorted_weyl_elements(W):
        _, projective_standard_multiplicity = standard_multiplicity(
            KL, verma_element, projective_element, base_weight
        )
        if projective_standard_multiplicity == 0:
            continue
        for simple_element in sorted_weyl_elements(W):
            _, verma_composition_multiplicity = standard_multiplicity(
                KL, verma_element, simple_element, base_weight
            )
            multiplicities[simple_element] += (
                projective_standard_multiplicity * verma_composition_multiplicity
            )

    return tuple(
        CompositionTerm(
            simple=word_label(element),
            reduced_word=word_tuple(element),
            length=int(element.length()),
            multiplicity=multiplicity,
        )
        for element, multiplicity in sorted(
            multiplicities.items(), key=lambda item: (item[0].length(), word_tuple(item[0]))
        )
        if multiplicity
    )


def projective_cover(
    W: Any,
    KL: Any,
    projective_element: Any,
    base_weight: str,
    include_composition: bool = False,
) -> ProjectiveCover:
    return ProjectiveCover(
        projective=word_label(projective_element),
        reduced_word=word_tuple(projective_element),
        length=int(projective_element.length()),
        verma_flag=verma_flag(W, KL, projective_element, base_weight),
        composition_factors=composition_factors(W, KL, projective_element, base_weight)
        if include_composition
        else (),
    )


def regular_block_projectives(
    cartan_type: Sequence[Any],
    projective_words: Iterable[str] | None = None,
    base_weight: str = "dominant",
    include_composition: bool = False,
) -> tuple[ProjectiveCover, ...]:
    if base_weight not in BASE_WEIGHT_CHOICES:
        raise ValueError(
            f"base_weight must be one of {', '.join(BASE_WEIGHT_CHOICES)}."
        )
    W = WeylGroup(list(cartan_type), prefix="s")
    KL = make_kl(W)
    if projective_words is None:
        elements = sorted_weyl_elements(W)
    else:
        elements = [parse_word(word, W) for word in projective_words]
    return tuple(
        projective_cover(W, KL, element, base_weight, include_composition)
        for element in elements
    )


def text_table(headers: Sequence[str], rows: Sequence[Sequence[Any]]) -> str:
    if not rows:
        return "  (none)"
    string_rows = [[str(item) for item in row] for row in rows]
    widths = [
        max(len(str(header)), *(len(row[index]) for row in string_rows))
        for index, header in enumerate(headers)
    ]
    header_line = "  " + "  ".join(
        str(header).ljust(widths[index]) for index, header in enumerate(headers)
    )
    rule = "  " + "  ".join("-" * width for width in widths)
    body = [
        "  "
        + "  ".join(row[index].ljust(widths[index]) for index in range(len(headers)))
        for row in string_rows
    ]
    return "\n".join([header_line, rule, *body])


def render_text(
    cartan_type: Sequence[Any], covers: Sequence[ProjectiveCover], base_weight: str
) -> str:
    lines = [
        f"Regular integral block of category O for Cartan type {format_cartan_type(cartan_type)}",
        "Dot action: w · λ = w(λ + ρ) - ρ",
        f"Base weight convention: {base_weight}",
        f"Convention: {convention_text(base_weight)}.",
        "",
    ]
    for index, cover in enumerate(covers):
        if index:
            lines.append("")
        lines.append(
            f"P({cover.projective}), projective cover of L({cover.projective}) "
            f"(length {cover.length})"
        )
        lines.append("Verma flag:")
        lines.append(
            text_table(
                ["Delta(x)", "word", "len", kl_header(base_weight), "mult"],
                [
                    (
                        f"Delta({term.verma})",
                        list(term.reduced_word) or "e",
                        term.length,
                        term.kl_polynomial,
                        term.multiplicity,
                    )
                    for term in cover.verma_flag
                ],
            )
        )
        if cover.composition_factors:
            lines.append("Composition factors:")
            lines.append(
                text_table(
                    ["L(y)", "word", "len", "mult"],
                    [
                        (
                            f"L({term.simple})",
                            list(term.reduced_word) or "e",
                            term.length,
                            term.multiplicity,
                        )
                        for term in cover.composition_factors
                    ],
                )
            )
    return "\n".join(lines)


def render_json(
    cartan_type: Sequence[Any], covers: Sequence[ProjectiveCover], base_weight: str
) -> str:
    return json.dumps(
        {
            "cartan_type": list(cartan_type),
            "block": "regular_integral",
            "base_weight": base_weight,
            "dot_action": "w · λ = w(λ + ρ) - ρ",
            "convention": convention_text(base_weight),
            "projective_covers": [cover.as_dict() for cover in covers],
        },
        indent=2,
        ensure_ascii=False,
    )


def render_csv(covers: Sequence[ProjectiveCover]) -> str:
    from io import StringIO

    buffer = StringIO()
    writer = csv.writer(buffer)
    writer.writerow(
        [
            "projective",
            "projective_reduced_word",
            "verma",
            "verma_reduced_word",
            "verma_length",
            "kl_polynomial",
            "multiplicity",
        ]
    )
    for cover in covers:
        for term in cover.verma_flag:
            writer.writerow(
                [
                    cover.projective,
                    " ".join(map(str, cover.reduced_word)),
                    term.verma,
                    " ".join(map(str, term.reduced_word)),
                    term.length,
                    term.kl_polynomial,
                    term.multiplicity,
                ]
            )
    return buffer.getvalue().rstrip()


def render_sl2_text(
    input_weight: Fraction, covers: Sequence[Sl2ProjectiveCover]
) -> str:
    first_cover = covers[0]
    block_weights = first_cover.block_weights
    lines = [
        "Block of category O for sl2",
        f"Input simple: L({format_weight(input_weight)})",
        "Dot action: s · λ = -λ - 2",
        f"Block type: {first_cover.block_kind}",
        "Simple objects: "
        + ", ".join(f"L({format_weight(weight)})" for weight in block_weights),
        "",
    ]

    for index, cover in enumerate(covers):
        if index:
            lines.append("")
        label = format_weight(cover.highest_weight)
        lines.append(f"P({label}), projective cover of L({label})")
        lines.append("Verma flag:")
        lines.append(
            text_table(
                ["Verma section", "mult"],
                [
                    (f"M({format_weight(term.weight)})", term.multiplicity)
                    for term in cover.verma_flag
                ],
            )
        )
        if cover.composition_factors:
            lines.append("Composition factors:")
            lines.append(
                text_table(
                    ["Simple", "mult"],
                    [
                        (f"L({format_weight(term.weight)})", term.multiplicity)
                        for term in cover.composition_factors
                    ],
                )
            )
        if cover.structure:
            lines.append("Structure:")
            lines.extend(f"  {statement}" for statement in cover.structure)
    return "\n".join(lines)


def render_sl2_json(
    input_weight: Fraction, covers: Sequence[Sl2ProjectiveCover]
) -> str:
    first_cover = covers[0]
    return json.dumps(
        {
            "lie_algebra": "sl2",
            "input_simple": f"L({format_weight(input_weight)})",
            "input_highest_weight": format_weight(input_weight),
            "dot_action": "s · λ = -λ - 2",
            "block_kind": first_cover.block_kind,
            "block_weights": [
                format_weight(weight) for weight in first_cover.block_weights
            ],
            "projective_covers": [cover.as_dict() for cover in covers],
        },
        indent=2,
        ensure_ascii=False,
    )


def render_sl2_csv(covers: Sequence[Sl2ProjectiveCover]) -> str:
    from io import StringIO

    buffer = StringIO()
    writer = csv.writer(buffer)
    writer.writerow(["projective", "verma", "multiplicity"])
    for cover in covers:
        projective = f"P({format_weight(cover.highest_weight)})"
        for term in cover.verma_flag:
            writer.writerow(
                [
                    projective,
                    f"M({format_weight(term.weight)})",
                    term.multiplicity,
                ]
            )
    return buffer.getvalue().rstrip()


def render_sln_text(
    n: int,
    input_weight: tuple[Fraction, ...],
    covers: Sequence[SlnProjectiveCover],
) -> str:
    first_cover = covers[0]
    block_weights = first_cover.block_weights
    basis = format_fundamental_weight_basis(n - 1)
    lines = [
        f"Regular integral block of category O for sl_{n}",
        f"Cartan type: {sln_cartan_type(n)}",
        f"Basis: fundamental weights {basis}",
        format_sln_rho(n - 1),
        "Dot action: w · λ = w(λ + ρ) - ρ",
        "Simple-reflection formulas:",
        *[f"  {formula}" for formula in sln_coordinate_formulas(n - 1)],
        f"Input simple: L({format_weight_vector(input_weight)})",
        f"Dominant block representative: {format_weight_vector(first_cover.dominant_weight)}",
        f"Block size: {len(block_weights)} simples",
    ]
    if len(block_weights) <= 24:
        lines.append(
            "Simple objects: "
            + ", ".join(f"L({format_weight_vector(weight)})" for weight in block_weights)
        )
    lines.append("")

    for index, cover in enumerate(covers):
        if index:
            lines.append("")
        label = format_weight_vector(cover.highest_weight)
        lines.append(f"P({label}), projective cover of L({label})")
        lines.append(f"Projective Weyl word: {list(cover.projective_word) or 'e'}")
        lines.append("Verma flag:")
        lines.append(
            text_table(
                ["Verma section", "word", "len", "KL", "mult"],
                [
                    (
                        f"M({format_weight_vector(term.weight)})",
                        list(term.reduced_word) or "e",
                        term.length,
                        term.kl_polynomial,
                        term.multiplicity,
                    )
                    for term in cover.verma_flag
                ],
            )
        )
        if cover.composition_factors:
            lines.append("Composition factors:")
            lines.append(
                text_table(
                    ["Simple", "word", "len", "mult"],
                    [
                        (
                            f"L({format_weight_vector(term.weight)})",
                            list(term.reduced_word) or "e",
                            term.length,
                            term.multiplicity,
                        )
                        for term in cover.composition_factors
                    ],
                )
            )
    return "\n".join(lines)


def render_sln_json(
    n: int,
    input_weight: tuple[Fraction, ...],
    covers: Sequence[SlnProjectiveCover],
) -> str:
    first_cover = covers[0]
    return json.dumps(
        {
            "lie_algebra": f"sl{n}",
            "cartan_type": sln_cartan_type(n),
            "basis": f"fundamental weights {format_fundamental_weight_basis(n - 1)}",
            "rho": format_sln_rho(n - 1),
            "dot_action": "w · λ = w(λ + ρ) - ρ",
            "simple_reflection_formulas": list(sln_coordinate_formulas(n - 1)),
            "input_simple": f"L({format_weight_vector(input_weight)})",
            "input_highest_weight": format_weight_vector(input_weight),
            "block": "regular_integral",
            "dominant_block_representative": format_weight_vector(
                first_cover.dominant_weight
            ),
            "block_weights": [
                format_weight_vector(weight) for weight in first_cover.block_weights
            ],
            "projective_covers": [cover.as_dict() for cover in covers],
        },
        indent=2,
        ensure_ascii=False,
    )


def render_sln_csv(covers: Sequence[SlnProjectiveCover]) -> str:
    from io import StringIO

    buffer = StringIO()
    writer = csv.writer(buffer)
    writer.writerow(
        [
            "projective",
            "projective_reduced_word",
            "verma",
            "verma_reduced_word",
            "verma_length",
            "kl_polynomial",
            "multiplicity",
        ]
    )
    for cover in covers:
        projective = f"P({format_weight_vector(cover.highest_weight)})"
        for term in cover.verma_flag:
            writer.writerow(
                [
                    projective,
                    " ".join(map(str, cover.projective_word)),
                    f"M({format_weight_vector(term.weight)})",
                    " ".join(map(str, term.reduced_word)),
                    term.length,
                    term.kl_polynomial,
                    term.multiplicity,
                ]
            )
    return buffer.getvalue().rstrip()


def build_weyl_parser(prog: str = "category_o_projectives.py weyl") -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Compute formal projective covers in a regular integral block of "
            "BGG category O using Sage's Kazhdan-Lusztig polynomials."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  .venv/bin/python category_o_projectives.py weyl A1 --projective e
  .venv/bin/python category_o_projectives.py weyl A2 --projective s1
  .venv/bin/python category_o_projectives.py weyl A3 --projective 1,2,1 --composition-factors
  .venv/bin/python category_o_projectives.py weyl G2 --all --format json

Notes:
  The word e means the identity element. A bare number like 1 means s1.
  The word w0 or longest means the longest Weyl group element.
  This venv's Sage launcher is useful for interactive Sage, but the Python
  launcher is the reliable way to pass CLI options to this script.
""",
    )
    parser.add_argument("cartan_type", type=parse_cartan_type, help="Finite type, e.g. A2.")
    projective_group = parser.add_mutually_exclusive_group()
    projective_group.add_argument(
        "--projective",
        "-p",
        action="append",
        default=None,
        help=(
            "Projective label w as a reduced word. May be passed multiple times. "
            "Defaults to e."
        ),
    )
    projective_group.add_argument(
        "--all", action="store_true", help="Compute projective covers for every w in W."
    )
    parser.add_argument(
        "--composition-factors",
        action="store_true",
        help="Also compute total simple composition multiplicities of P(w).",
    )
    parser.add_argument(
        "--base-weight",
        choices=BASE_WEIGHT_CHOICES,
        default="dominant",
        help=(
            "Indexing convention for the regular integral block. "
            "dominant uses [Delta(x):L(w)] = P_{x,w}(1); "
            "antidominant uses [Delta(x):L(w)] = P_{w,x}(1)."
        ),
    )
    parser.add_argument(
        "--format",
        choices=("text", "json", "csv"),
        default="text",
        help="Output format.",
    )
    return parser


def build_sl2_parser(prog: str = "category_o_projectives.py sl2") -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Compute projective covers in sl2 category O using actual "
            "highest-weight labels, as in P(λ), L(λ), and M(λ)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  .venv/bin/python category_o_projectives.py sl2 0
  .venv/bin/python category_o_projectives.py sl2 2 --composition-factors
  .venv/bin/python category_o_projectives.py sl2 -2
  .venv/bin/python category_o_projectives.py sl2 0 --all-in-block

Notes:
  For sl2, the dot action is s · λ = -λ - 2.
  If n >= 0, the regular integral block has simples L(n) and L(-n-2).
""",
    )
    parser.add_argument(
        "highest_weight",
        type=parse_highest_weight,
        help="Highest weight λ in L(λ), e.g. 0, 2, -2, or 1/2.",
    )
    parser.add_argument(
        "--all-in-block",
        action="store_true",
        help="Compute the projective covers of every simple in the same sl2 block.",
    )
    parser.add_argument(
        "--composition-factors",
        action="store_true",
        help="Also show total simple composition multiplicities.",
    )
    parser.add_argument(
        "--format",
        choices=("text", "json", "csv"),
        default="text",
        help="Output format.",
    )
    return parser


def build_sln_parser(
    fixed_n: int | None = None, prog: str = "category_o_projectives.py sln"
) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Compute projective covers in regular integral sl_n blocks using "
            "fundamental-weight coordinates."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  .venv/bin/python category_o_projectives.py sln 3 0
  .venv/bin/python category_o_projectives.py sln 3 1,0 --composition-factors
  .venv/bin/python category_o_projectives.py sl3 1,0 --all-in-block
  .venv/bin/python category_o_projectives.py sl4 0,1,0

Notes:
  Coordinates are in the fundamental-weight basis.
  For example, sl3 1,0 means ω_1, and sl4 0,1,0 means ω_2.
  The dot action is w · λ = w(λ + ρ) - ρ.
  The output also prints the coordinate formula for each simple reflection.
  This mode currently supports regular integral blocks.
""",
    )
    if fixed_n is None:
        parser.add_argument("n", type=int, help="The n in sl_n.")
    parser.add_argument(
        "highest_weight",
        help=(
            "Highest weight in fundamental coordinates. Use 0 for the zero "
            "weight, or comma-separated labels like 1,0 or 0,1,0."
        ),
    )
    parser.add_argument(
        "--all-in-block",
        action="store_true",
        help="Compute the projective covers of every simple in the same block.",
    )
    parser.add_argument(
        "--composition-factors",
        action="store_true",
        help="Also show total simple composition multiplicities.",
    )
    parser.add_argument(
        "--format",
        choices=("text", "json", "csv"),
        default="text",
        help="Output format.",
    )
    return parser


def top_level_help() -> str:
    return """usage:
  category_o_projectives.py sl2 <highest-weight> [options]
  category_o_projectives.py sln <n> <highest-weight-vector> [options]
  category_o_projectives.py sl3 <highest-weight-vector> [options]
  category_o_projectives.py <highest-weight> [options]
  category_o_projectives.py weyl <cartan-type> [options]

Friendly sl2 examples:
  ./category_o_projectives.py sl2 0
  ./category_o_projectives.py sl2 2 --composition-factors
  ./category_o_projectives.py sl2 -2
  ./category_o_projectives.py 0 --all-in-block

Friendly sl_n examples:
  ./category_o_projectives.py sln 3 0
  ./category_o_projectives.py sl3 1,0 --composition-factors
  ./category_o_projectives.py sl4 0,1,0

General Weyl/Kazhdan-Lusztig examples:
  ./category_o_projectives.py weyl A2 --projective s1
  ./category_o_projectives.py weyl A3 --projective 1,2,1
  ./category_o_projectives.py weyl G2 --all --format json

Use one of these for detailed options:
  ./category_o_projectives.py sl2 --help
  ./category_o_projectives.py sln --help
  ./category_o_projectives.py weyl --help
"""


def run_weyl(argv: Sequence[str], prog: str = "category_o_projectives.py weyl") -> int:
    parser = build_weyl_parser(prog=prog)
    args = parser.parse_args(argv)

    projective_words = None if args.all else (args.projective or ["e"])
    try:
        covers = regular_block_projectives(
            args.cartan_type,
            projective_words=projective_words,
            base_weight=args.base_weight,
            include_composition=args.composition_factors,
        )
    except Exception as exc:
        parser.exit(2, f"error: {exc}\n")

    if args.format == "json":
        print(render_json(args.cartan_type, covers, args.base_weight))
    elif args.format == "csv":
        print(render_csv(covers))
    else:
        print(render_text(args.cartan_type, covers, args.base_weight))
    return 0


def run_sl2(argv: Sequence[str], prog: str = "category_o_projectives.py sl2") -> int:
    parser = build_sl2_parser(prog=prog)
    args = parser.parse_args(argv)
    covers = sl2_projective_covers(
        args.highest_weight,
        include_composition=args.composition_factors,
        all_in_block=args.all_in_block,
    )

    if args.format == "json":
        print(render_sl2_json(args.highest_weight, covers))
    elif args.format == "csv":
        print(render_sl2_csv(covers))
    else:
        print(render_sl2_text(args.highest_weight, covers))
    return 0


def run_sln(
    argv: Sequence[str],
    fixed_n: int | None = None,
    prog: str = "category_o_projectives.py sln",
) -> int:
    parser = build_sln_parser(fixed_n=fixed_n, prog=prog)
    args = parser.parse_args(argv)
    n = fixed_n if fixed_n is not None else args.n
    try:
        highest_weight = parse_weight_vector(args.highest_weight, n - 1)
        covers = sln_projective_covers(
            n,
            highest_weight,
            include_composition=args.composition_factors,
            all_in_block=args.all_in_block,
        )
    except Exception as exc:
        parser.exit(2, f"error: {exc}\n")

    if args.format == "json":
        print(render_sln_json(n, highest_weight, covers))
    elif args.format == "csv":
        print(render_sln_csv(covers))
    else:
        print(render_sln_text(n, highest_weight, covers))
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] in {"-h", "--help"}:
        print(top_level_help())
        return 0

    command, *rest = argv
    if command == "sl2":
        return run_sl2(rest)
    if command == "sln":
        return run_sln(rest)
    sln_match = re.fullmatch(r"sl(\d+)", command)
    if sln_match:
        n = int(sln_match.group(1))
        if n == 2:
            return run_sl2(rest, prog=f"category_o_projectives.py {command}")
        return run_sln(rest, fixed_n=n, prog=f"category_o_projectives.py {command}")
    if command == "weyl":
        return run_weyl(rest)
    if looks_like_highest_weight(command):
        return run_sl2(argv, prog="category_o_projectives.py")

    return run_weyl(argv, prog="category_o_projectives.py")


if __name__ == "__main__":
    raise SystemExit(main())
