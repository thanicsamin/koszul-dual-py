#!/usr/bin/env python3
"""
Koszul-dual computations for Lie algebra representations.

This script is intentionally explicit: it computes the finite graded algebras
appearing in the notes by their homogeneous bases and multiplication tables.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from dataclasses import dataclass
from fractions import Fraction
from typing import Sequence

import category_o_projectives as category_o


@dataclass(frozen=True)
class BasisElement:
    name: str
    degree: int
    source: str | None = None
    target: str | None = None

    def as_dict(self) -> dict[str, object]:
        payload: dict[str, object] = {"name": self.name, "degree": self.degree}
        if self.source is not None:
            payload["source"] = self.source
        if self.target is not None:
            payload["target"] = self.target
        return payload


class FiniteGradedAlgebra:
    def __init__(
        self,
        name: str,
        basis: Sequence[BasisElement],
        products: dict[tuple[str, str], str],
        zero_products: set[tuple[str, str]] | None = None,
        notes: Sequence[str] = (),
    ) -> None:
        self.name = name
        self.basis = tuple(basis)
        self.products = dict(products)
        self.zero_products = set(zero_products or set())
        self.notes = tuple(notes)
        self._basis_by_name = {element.name: element for element in self.basis}

    def degree_pieces(self) -> dict[int, list[str]]:
        pieces: dict[int, list[str]] = {}
        for element in self.basis:
            pieces.setdefault(element.degree, []).append(element.name)
        return dict(sorted(pieces.items()))

    def hilbert_series(self) -> str:
        terms: list[str] = []
        for degree, names in self.degree_pieces().items():
            coefficient = len(names)
            if degree == 0:
                terms.append(str(coefficient))
            elif coefficient == 1:
                terms.append("t" if degree == 1 else f"t^{degree}")
            else:
                terms.append(f"{coefficient}t" if degree == 1 else f"{coefficient}t^{degree}")
        return " + ".join(terms)

    def named_products(self) -> list[tuple[str, str, str]]:
        return [
            (left, right, result)
            for (left, right), result in sorted(self.products.items())
        ]

    def named_zero_products(self) -> list[tuple[str, str]]:
        return sorted(self.zero_products)

    def as_dict(self) -> dict[str, object]:
        return {
            "name": self.name,
            "basis": [element.as_dict() for element in self.basis],
            "degree_pieces": {
                str(degree): names for degree, names in self.degree_pieces().items()
            },
            "hilbert_series": self.hilbert_series(),
            "nonzero_products": [
                {"left": left, "right": right, "result": result}
                for left, right, result in self.named_products()
            ],
            "declared_zero_products": [
                {"left": left, "right": right}
                for left, right in self.named_zero_products()
            ],
            "notes": list(self.notes),
        }


@dataclass(frozen=True)
class DualityEntry:
    projective_label: str
    projective_word: tuple[int, ...]
    dual_simple_label: str
    dual_word: tuple[int, ...]

    def as_dict(self) -> dict[str, object]:
        return {
            "projective": self.projective_label,
            "projective_reduced_word": list(self.projective_word),
            "dual_simple": self.dual_simple_label,
            "dual_reduced_word": list(self.dual_word),
        }


@dataclass(frozen=True)
class RegularBlockDuality:
    input_label: str
    cartan_type: str
    langlands_dual_type: str
    convention: str
    entries: tuple[DualityEntry, ...]
    notes: tuple[str, ...] = ()

    def as_dict(self) -> dict[str, object]:
        return {
            "input": self.input_label,
            "cartan_type": self.cartan_type,
            "langlands_dual_type": self.langlands_dual_type,
            "convention": self.convention,
            "entries": [entry.as_dict() for entry in self.entries],
            "notes": list(self.notes),
        }


@dataclass(frozen=True)
class Sl2FunctorObject:
    module_kind: str
    weight: Fraction
    shift: int = 0
    dual_verma: bool = False

    def module_label(self) -> str:
        label = f"{self.module_kind}({category_o.format_weight(self.weight)})"
        if self.dual_verma:
            label += "ᵛ"
        return label

    def shifted_label(self) -> str:
        return format_shifted(self.module_label(), self.shift)


@dataclass(frozen=True)
class SlnFunctorObject:
    module_kind: str
    weight: tuple[Fraction, ...]
    shift: int = 0

    def module_label(self) -> str:
        return f"{self.module_kind}({category_o.format_weight_vector(self.weight)})"

    def shifted_label(self) -> str:
        return format_shifted(self.module_label(), self.shift)


def principal_block_algebra() -> FiniteGradedAlgebra:
    """The algebra A = End_O(P(-2) ⊕ P(0)) from main.tex."""

    basis = [
        BasisElement("e_0", 0, "0", "0"),
        BasisElement("e_-2", 0, "-2", "-2"),
        BasisElement("i", 1, "0", "-2"),
        BasisElement("p", 1, "-2", "0"),
        BasisElement("ip", 2, "-2", "-2"),
    ]
    products = {
        ("e_0", "e_0"): "e_0",
        ("e_-2", "e_-2"): "e_-2",
        ("e_-2", "i"): "i",
        ("i", "e_0"): "i",
        ("e_0", "p"): "p",
        ("p", "e_-2"): "p",
        ("i", "p"): "ip",
        ("e_-2", "ip"): "ip",
        ("ip", "e_-2"): "ip",
    }
    zero_products = {
        ("e_0", "e_-2"),
        ("e_-2", "e_0"),
        ("p", "i"),
        ("ip", "i"),
        ("p", "ip"),
        ("ip", "ip"),
    }
    return FiniteGradedAlgebra(
        "A = End_O(P(-2) ⊕ P(0))",
        basis,
        products,
        zero_products,
        notes=(
            "This is the principal block algebra for category O of sl2.",
            "The relation p i = 0 and the nonzero product i p span the socle of P(-2).",
            "Degrees: e_0,e_-2 in degree 0; i,p in degree 1; ip in degree 2.",
        ),
    )


def ext_algebra_of_principal_block() -> FiniteGradedAlgebra:
    """B = Ext_A^*(L,L), equivalently Ext_O^*(L,L), from main.tex."""

    basis = [
        BasisElement("ε_0", 0, "0", "0"),
        BasisElement("ε_-2", 0, "-2", "-2"),
        BasisElement("α", 1, "0", "-2"),
        BasisElement("β", 1, "-2", "0"),
        BasisElement("γ", 2, "0", "0"),
    ]
    products = {
        ("ε_0", "ε_0"): "ε_0",
        ("ε_-2", "ε_-2"): "ε_-2",
        ("ε_-2", "α"): "α",
        ("α", "ε_0"): "α",
        ("ε_0", "β"): "β",
        ("β", "ε_-2"): "β",
        ("β", "α"): "γ",
        ("ε_0", "γ"): "γ",
        ("γ", "ε_0"): "γ",
    }
    zero_products = {
        ("ε_0", "ε_-2"),
        ("ε_-2", "ε_0"),
        ("α", "β"),
        ("α", "γ"),
        ("γ", "β"),
        ("γ", "γ"),
    }
    return FiniteGradedAlgebra(
        "B = Ext_A^*(L,L)",
        basis,
        products,
        zero_products,
        notes=(
            "Here L = L(0) ⊕ L(-2).",
            "α ∈ Ext^1(L(0), L(-2)), β ∈ Ext^1(L(-2), L(0)), γ ∈ Ext^2(L(0), L(0)).",
            "Yoneda products: β α = γ and α β = 0.",
        ),
    )


def koszul_dual_of_ext_algebra() -> FiniteGradedAlgebra:
    """E(B) = Ext_B^*(B_0,B_0), from the self-duality section."""

    basis = [
        BasisElement("η_0", 0, "0", "0"),
        BasisElement("η_-2", 0, "-2", "-2"),
        BasisElement("u", 1, "0", "-2"),
        BasisElement("v", 1, "-2", "0"),
        BasisElement("w", 2, "-2", "-2"),
    ]
    products = {
        ("η_0", "η_0"): "η_0",
        ("η_-2", "η_-2"): "η_-2",
        ("η_-2", "u"): "u",
        ("u", "η_0"): "u",
        ("η_0", "v"): "v",
        ("v", "η_-2"): "v",
        ("u", "v"): "w",
        ("η_-2", "w"): "w",
        ("w", "η_-2"): "w",
    }
    zero_products = {
        ("η_0", "η_-2"),
        ("η_-2", "η_0"),
        ("v", "u"),
        ("w", "u"),
        ("v", "w"),
        ("w", "w"),
    }
    return FiniteGradedAlgebra(
        "E(B) = Ext_B^*(B_0,B_0)",
        basis,
        products,
        zero_products,
        notes=(
            "u ∈ Ext^1(S(0), S(-2)), v ∈ Ext^1(S(-2), S(0)), w ∈ Ext^2(S(-2), S(-2)).",
            "Yoneda products: u v = w and v u = 0.",
            "This has the same multiplication table as B after swapping the two vertices.",
        ),
    )


def vertex_swap_isomorphism() -> dict[str, str]:
    return {
        "ε_0": "η_-2",
        "ε_-2": "η_0",
        "α": "v",
        "β": "u",
        "γ": "w",
    }


def principal_block_resolutions() -> dict[str, list[str]]:
    return {
        "right_A_modules": [
            "0 -> Q(0) -> Q(-2) -> L(-2) -> 0",
            "0 -> Q(0) -> Q(-2) -> Q(0) -> L(0) -> 0",
        ],
        "left_A_koszul_resolutions": [
            "0 -> P(0)<1> --·p--> P(-2) -> L(-2) -> 0",
            "0 -> P(0)<2> --·p--> P(-2)<1> --·i--> P(0) -> L(0) -> 0",
        ],
        "left_B_koszul_resolutions": [
            "0 -> R(-2)<1> --·α--> R(0) -> S(0) -> 0",
            "0 -> R(-2)<2> --·α--> R(0)<1> --·β--> R(-2) -> S(-2) -> 0",
        ],
    }


def principal_block_duality_functor_images() -> dict[str, str]:
    return {
        "K(M(0)[0])": "L(-2)[0]",
        "K(M(-2)[0]) = K(L(-2)[0])": "P(0)[0]",
        "K(L(0)[0])": "P(-2)[0]",
        "K(P(-2)[0])": "L(0)[0]",
    }


def format_shifted(module_label: str, shift: int) -> str:
    return f"{module_label}[{shift}]"


def sl2_weight_label(weight: Fraction) -> str:
    return category_o.format_weight(weight)


def sl2_module_label(module_kind: str, weight: Fraction, dual_verma: bool = False) -> str:
    label = f"{module_kind}({sl2_weight_label(weight)})"
    if dual_verma:
        label += "ᵛ"
    return label


def sl2_shifted_module_label(
    module_kind: str, weight: Fraction, shift: int, dual_verma: bool = False
) -> str:
    return format_shifted(sl2_module_label(module_kind, weight, dual_verma), shift)


def parse_sl2_shifted_object(label: str) -> Sl2FunctorObject:
    normalized = label.strip().replace(" ", "")
    if normalized.startswith("K(") and normalized.endswith(")"):
        normalized = normalized[2:-1]

    shift = 0
    shift_match = re.fullmatch(r"(.+)\[([+-]?\d+)\]", normalized)
    if shift_match:
        normalized = shift_match.group(1)
        shift = int(shift_match.group(2))
    elif "[" in normalized or "]" in normalized:
        raise ValueError(
            f"could not parse cohomological shift in {label!r}; use an integer shift like [0] or [-1]."
        )

    normalized = normalized.replace("^vee", "ᵛ")
    normalized = normalized.replace("^v", "ᵛ")
    normalized = normalized.replace("^∨", "ᵛ")
    normalized = normalized.replace("∨", "ᵛ")
    normalized = normalized.replace("^*", "*")

    dual_verma = False
    for suffix in ("ᵛ", "vee", "dual", "v", "*"):
        if normalized.endswith(suffix):
            normalized = normalized[: -len(suffix)]
            dual_verma = True
            break

    if not normalized or normalized[0] not in {"M", "L", "P"}:
        allowed = "M(0)[0], M0[0], M(0)v[0], L(0)[0], P(-2)[0]"
        raise ValueError(f"unknown sl2 object {label!r}; try {allowed}.")

    module_kind = normalized[0]
    weight_text = normalized[1:]
    if weight_text.startswith("(") and weight_text.endswith(")"):
        weight_text = weight_text[1:-1]
    if not weight_text:
        raise ValueError(f"missing highest weight in {label!r}.")
    if dual_verma and module_kind != "M":
        raise ValueError(
            f"the dual-Verma suffix is only supported for M(...), not {module_kind}(...)."
        )

    try:
        weight = category_o.parse_highest_weight(weight_text)
    except Exception as exc:
        raise ValueError(str(exc)) from exc

    if category_o.sl2_block_kind(weight) != "regular integral block":
        raise ValueError(
            "sl2 currently supports regular integral sl2 blocks only; "
            f"{sl2_weight_label(weight)} is in a {category_o.sl2_block_kind(weight)}."
        )

    return Sl2FunctorObject(module_kind, weight, shift, dual_verma)


def normalize_sl2_object_label(label: str) -> str:
    return parse_sl2_shifted_object(label).module_label()


def sl2_block_pair(weight: Fraction) -> tuple[Fraction, Fraction]:
    dominant, antidominant = category_o.sl2_block_weights(weight)
    return dominant, antidominant


def sl2_opposite_block_weight(weight: Fraction) -> Fraction:
    dominant, antidominant = sl2_block_pair(weight)
    return antidominant if weight == dominant else dominant


def sl2_functor_input_label(spec: Sl2FunctorObject) -> str:
    base = f"K({spec.shifted_label()})"
    dominant, antidominant = sl2_block_pair(spec.weight)
    if spec.module_kind == "M" and not spec.dual_verma and spec.weight == antidominant:
        return f"{base} = K({sl2_shifted_module_label('L', spec.weight, spec.shift)})"
    if spec.module_kind == "M" and spec.dual_verma and spec.weight == antidominant:
        return f"{base} = K({sl2_shifted_module_label('L', spec.weight, spec.shift)})"
    if spec.module_kind == "P" and spec.weight == dominant:
        return f"{base} = K({sl2_shifted_module_label('M', spec.weight, spec.shift)})"
    return base


def sl2_concentrated_functor_image(spec: Sl2FunctorObject) -> dict[str, object]:
    dominant, antidominant = sl2_block_pair(spec.weight)
    target_shift = -spec.shift
    input_label = sl2_functor_input_label(spec)

    if spec.module_kind == "M" and not spec.dual_verma:
        if spec.weight == dominant:
            return {
                "kind": "concentrated",
                "input": input_label,
                "result": sl2_shifted_module_label("L", antidominant, target_shift),
            }
        return {
            "kind": "concentrated",
            "input": input_label,
            "result": sl2_shifted_module_label("P", dominant, target_shift),
        }

    if spec.module_kind == "L":
        target_weight = sl2_opposite_block_weight(spec.weight)
        return {
            "kind": "concentrated",
            "input": input_label,
            "result": sl2_shifted_module_label("P", target_weight, target_shift),
        }

    if spec.module_kind == "P":
        target_weight = sl2_opposite_block_weight(spec.weight)
        return {
            "kind": "concentrated",
            "input": input_label,
            "result": sl2_shifted_module_label("L", target_weight, target_shift),
        }

    # The antidominant dual Verma is already the simple L(antidominant).
    return {
        "kind": "concentrated",
        "input": input_label,
        "result": sl2_shifted_module_label("P", dominant, target_shift),
    }


def sl2_functor_image(label: str) -> dict[str, object]:
    spec = parse_sl2_shifted_object(label)
    dominant, antidominant = sl2_block_pair(spec.weight)
    if not spec.dual_verma or spec.weight == antidominant:
        return sl2_concentrated_functor_image(spec)

    target_shift = -spec.shift
    graded_terms = [
        (target_shift, sl2_module_label("L", dominant)),
        (target_shift + 1, sl2_module_label("L", dominant)),
        (target_shift + 2, sl2_module_label("L", antidominant)),
    ]
    graded_object = " ⊕ ".join(
        format_shifted(label, degree) for degree, label in graded_terms
    )
    return {
        "kind": "cohomology",
        "input": sl2_functor_input_label(spec),
        "result": graded_object,
        "graded_object": graded_object,
        "projective_resolution": (
            f"0 -> P({sl2_weight_label(dominant)}) -> P({sl2_weight_label(antidominant)}) "
            f"-> P({sl2_weight_label(antidominant)}) -> "
            f"M({sl2_weight_label(dominant)})ᵛ -> 0"
        ),
        "raw_ext_groups": [
            f"Ext^0(M({sl2_weight_label(dominant)})ᵛ, L({sl2_weight_label(antidominant)})) ≅ C",
            f"Ext^1(M({sl2_weight_label(dominant)})ᵛ, L({sl2_weight_label(antidominant)})) ≅ C",
            f"Ext^2(M({sl2_weight_label(dominant)})ᵛ, L({sl2_weight_label(dominant)})) ≅ C",
        ],
        "target_cohomology_after_vertex_swap": [
            f"H^{degree} ≅ {label}" for degree, label in graded_terms
        ],
    }


def sl2_functor_images(labels: Sequence[str] | None = None) -> list[dict[str, object]]:
    requested = labels or ["M(0)", "M(-2)", "L(0)", "P(-2)", "M(0)ᵛ"]
    return [sl2_functor_image(label) for label in requested]


def format_direct_sum(labels: Sequence[str]) -> str:
    return " ⊕ ".join(labels)


def sl2_l_block_from_image(image: dict[str, object]) -> str | None:
    input_text = str(image["input"])
    match = re.search(r"[MLP]\((-?\d+(?:/\d+)?)\)", input_text)
    if not match:
        return None
    weight = category_o.parse_highest_weight(match.group(1))
    return format_direct_sum(
        f"L({sl2_weight_label(block_weight)})"
        for block_weight in category_o.sl2_block_weights(weight)
    )


def sl2_l_block_description(images: Sequence[dict[str, object]]) -> str:
    descriptions = []
    for image in images:
        description = sl2_l_block_from_image(image)
        if description is not None and description not in descriptions:
            descriptions.append(description)
    if len(descriptions) == 1:
        return descriptions[0]
    if descriptions:
        return "; ".join(descriptions)
    return "the sum of simples in the relevant block"


def langlands_dual_cartan_type(cartan_type: Sequence[object]) -> list[object]:
    letter = str(cartan_type[0]).upper()
    rest = list(cartan_type[1:])
    if letter == "B":
        return ["C", *rest]
    if letter == "C":
        return ["B", *rest]
    return [letter, *rest]


def dual_element(w0: object, element: object, dual_map: str) -> object:
    if dual_map == "left":
        return w0 * element
    if dual_map == "right":
        return element * w0
    raise ValueError(f"unknown dual map convention: {dual_map}")


def dual_map_convention(dual_map: str) -> str:
    if dual_map == "left":
        return "P(w · λ) ↦ L((w0 w) · λᵛ)"
    if dual_map == "right":
        return "P(w · λ) ↦ L((w w0) · λᵛ)"
    raise ValueError(f"unknown dual map convention: {dual_map}")


def regular_block_duality_by_cartan_type(
    cartan_type: Sequence[object], dual_map: str = "left"
) -> RegularBlockDuality:
    W = category_o.WeylGroup(list(cartan_type), prefix="s")
    w0 = W.long_element()
    entries: list[DualityEntry] = []
    for element in category_o.sorted_weyl_elements(W):
        target = dual_element(w0, element, dual_map)
        source_label = category_o.word_label(element)
        target_label = category_o.word_label(target)
        entries.append(
            DualityEntry(
                projective_label=f"P({source_label})",
                projective_word=category_o.word_tuple(element),
                dual_simple_label=f"L({target_label})",
                dual_word=category_o.word_tuple(target),
            )
        )

    dual_type = langlands_dual_cartan_type(cartan_type)
    return RegularBlockDuality(
        input_label=f"regular integral block of type {category_o.format_cartan_type(cartan_type)}",
        cartan_type=category_o.format_cartan_type(cartan_type),
        langlands_dual_type=category_o.format_cartan_type(dual_type),
        convention=dual_map_convention(dual_map),
        entries=tuple(entries),
        notes=(
            "BGS Koszul duality identifies a regular integral block with the regular block for the Langlands-dual Lie algebra.",
            "For types A, D, E, F4, and G2 the Cartan type is self-dual; types B and C are interchanged.",
            "This reports the dual block and projective-to-simple correspondence, not the full higher-rank Ext multiplication table.",
        ),
    )


def regular_sln_duality(
    n: int,
    highest_weight: tuple[object, ...],
    dual_map: str = "left",
) -> RegularBlockDuality:
    if n < 2:
        raise ValueError("sl_n mode requires n >= 2.")
    if not category_o.sln_is_regular_integral(highest_weight):
        raise ValueError(
            "regular sl_n Koszul duality mode currently supports regular integral blocks only."
        )

    W = category_o.WeylGroup(["A", n - 1], prefix="s")
    w0 = W.long_element()
    dominant_weight = category_o.sln_dominant_representative(highest_weight)
    entries: list[DualityEntry] = []
    for element in category_o.sorted_weyl_elements(W):
        target = dual_element(w0, element, dual_map)
        source_weight = category_o.sln_dot_action(dominant_weight, element)
        target_weight = category_o.sln_dot_action(dominant_weight, target)
        entries.append(
            DualityEntry(
                projective_label=f"P({category_o.format_weight_vector(source_weight)})",
                projective_word=category_o.word_tuple(element),
                dual_simple_label=f"L({category_o.format_weight_vector(target_weight)})",
                dual_word=category_o.word_tuple(target),
            )
        )

    return RegularBlockDuality(
        input_label=(
            f"regular integral block of sl_{n} containing "
            f"L({category_o.format_weight_vector(highest_weight)})"
        ),
        cartan_type=f"A{n - 1}",
        langlands_dual_type=f"A{n - 1}",
        convention=dual_map_convention(dual_map),
        entries=tuple(entries),
        notes=(
            f"sl_{n} is Langlands self-dual, so the dual Cartan type is A{n - 1}.",
            f"Dominant representative: {category_o.format_weight_vector(dominant_weight)}.",
            "Entries are printed in fundamental-weight coordinates.",
            "This reports the dual block and projective-to-simple correspondence, not the full higher-rank Ext multiplication table.",
        ),
    )


def is_sln_object_text(value: str) -> bool:
    normalized = value.strip().replace(" ", "")
    if normalized.startswith("K(") and normalized.endswith(")"):
        normalized = normalized[2:-1]
    return bool(normalized) and normalized[0] in {"M", "L", "P"}


def parse_sln_shifted_object(label: str, n: int) -> SlnFunctorObject:
    rank = n - 1
    normalized = label.strip().replace(" ", "")
    if normalized.startswith("K(") and normalized.endswith(")"):
        normalized = normalized[2:-1]

    shift = 0
    shift_match = re.fullmatch(r"(.+)\[([+-]?\d+)\]", normalized)
    if shift_match:
        normalized = shift_match.group(1)
        shift = int(shift_match.group(2))
    elif "[" in normalized or "]" in normalized:
        raise ValueError(
            f"could not parse cohomological shift in {label!r}; use an integer shift like [0] or [-1]."
        )

    if not normalized or normalized[0] not in {"M", "L", "P"}:
        raise ValueError(f"unknown sl_{n} object {label!r}; try P(0)[0] or L(0)[0].")

    module_kind = normalized[0]
    if module_kind == "M":
        raise ValueError(
            "higher-rank sl_n object mode currently supports P(...) and L(...). "
            "Use a bare weight, such as 0 or 0,1, to print the whole block correspondence."
        )

    weight_text = normalized[1:]
    if weight_text.startswith("(") and weight_text.endswith(")"):
        weight_text = weight_text[1:-1]
    if not weight_text:
        raise ValueError(f"missing highest weight in {label!r}.")

    try:
        weight = category_o.parse_weight_vector(weight_text, rank)
    except Exception as exc:
        raise ValueError(str(exc)) from exc
    if not category_o.sln_is_regular_integral(weight):
        raise ValueError(
            f"sl_{n} object mode currently supports regular integral blocks only; "
            f"{category_o.format_weight_vector(weight)} is not regular integral."
        )
    return SlnFunctorObject(module_kind, weight, shift)


def sln_shifted_module_label(
    module_kind: str, weight: Sequence[Fraction], shift: int
) -> str:
    return format_shifted(
        f"{module_kind}({category_o.format_weight_vector(weight)})", shift
    )


def sln_functor_image(
    n: int, label: str, dual_map: str = "left"
) -> dict[str, object]:
    spec = parse_sln_shifted_object(label, n)
    W = category_o.WeylGroup(["A", n - 1], prefix="s")
    w0 = W.long_element()
    dominant_weight = category_o.sln_dominant_representative(spec.weight)
    element = category_o.sln_find_weyl_element(W, dominant_weight, spec.weight)
    target = dual_element(w0, element, dual_map)
    target_weight = category_o.sln_dot_action(dominant_weight, target)
    target_kind = "L" if spec.module_kind == "P" else "P"
    return {
        "kind": "concentrated",
        "input": f"K({spec.shifted_label()})",
        "result": sln_shifted_module_label(target_kind, target_weight, -spec.shift),
        "source_weight": category_o.format_weight_vector(spec.weight),
        "source_reduced_word": list(category_o.word_tuple(element)),
        "target_weight": category_o.format_weight_vector(target_weight),
        "target_reduced_word": list(category_o.word_tuple(target)),
        "dominant_weight": category_o.format_weight_vector(dominant_weight),
        "cartan_type": f"A{n - 1}",
        "convention": dual_map_convention(dual_map),
    }


def sln_functor_images(
    n: int, labels: Sequence[str], dual_map: str = "left"
) -> list[dict[str, object]]:
    return [sln_functor_image(n, label, dual_map=dual_map) for label in labels]


def render_sln_functor_images(
    n: int, images: Sequence[dict[str, object]], format_name: str
) -> str:
    l_block = sln_l_block_description(n, images)
    if format_name == "json":
        return json.dumps(
            {
                "functor": "K = RHom_O(-, L_block)",
                "lie_algebra": f"sl_{n}",
                "L_block": l_block,
                "images": images,
            },
            indent=2,
            ensure_ascii=False,
        )

    lines = [
        f"Koszul duality functor for regular integral blocks of category O for sl_{n}",
        "K = RHom_O(-, L_block), where L_block is the sum of simples in the relevant block",
        f"Here, L_block = {l_block}",
        "",
    ]
    lines.extend(f"{image['input']} = {image['result']}" for image in images)
    return "\n".join(lines)


def sln_l_block_description(n: int, images: Sequence[dict[str, object]]) -> str:
    descriptions: list[str] = []
    W = category_o.WeylGroup(["A", n - 1], prefix="s")
    for image in images:
        dominant_text = image.get("dominant_weight")
        if not isinstance(dominant_text, str):
            continue
        dominant_weight = category_o.parse_weight_vector(dominant_text, n - 1)
        labels = [
            f"L({category_o.format_weight_vector(category_o.sln_dot_action(dominant_weight, element))})"
            for element in category_o.sorted_weyl_elements(W)
        ]
        description = format_direct_sum(labels)
        if description not in descriptions:
            descriptions.append(description)
    if len(descriptions) == 1:
        return descriptions[0]
    if descriptions:
        return "; ".join(descriptions)
    return "the sum of simples in the relevant block"


def zero_sln_weight(n: int) -> tuple[Fraction, ...]:
    return tuple(Fraction(0) for _ in range(n - 1))


def render_sln_inputs(
    n: int,
    inputs: Sequence[str],
    dual_map: str,
    format_name: str,
) -> str:
    if n < 2:
        raise ValueError("sl_n mode requires n >= 2.")
    if not inputs:
        return render_regular_block_duality(
            regular_sln_duality(n, zero_sln_weight(n), dual_map=dual_map),
            format_name,
        )

    object_flags = [is_sln_object_text(input_text) for input_text in inputs]
    if any(object_flags):
        if not all(object_flags):
            raise ValueError(
                "Do not mix bare block weights with object inputs. "
                "Use either sl3 0,1 or sl3 'P(0,1)[0]' 'L(-2,2)[1]'."
            )
        return render_sln_functor_images(
            n, sln_functor_images(n, inputs, dual_map=dual_map), format_name
        )

    if len(inputs) > 1:
        raise ValueError(
            "Expected one bare highest weight for a block table, or object inputs like P(0,1)[0]."
        )

    highest_weight = category_o.parse_weight_vector(inputs[0], n - 1)
    return render_regular_block_duality(
        regular_sln_duality(n, highest_weight, dual_map=dual_map),
        format_name,
    )


def render_algebra(algebra: FiniteGradedAlgebra) -> str:
    lines = [
        algebra.name,
        f"Hilbert series: {algebra.hilbert_series()}",
        "Basis by degree:",
    ]
    for degree, names in algebra.degree_pieces().items():
        lines.append(f"  degree {degree}: {', '.join(names)}")
    lines.append("Nonzero products:")
    for left, right, result in algebra.named_products():
        lines.append(f"  {left} {right} = {result}")
    if algebra.zero_products:
        lines.append("Key zero products:")
        for left, right in algebra.named_zero_products():
            lines.append(f"  {left} {right} = 0")
    if algebra.notes:
        lines.append("Notes:")
        lines.extend(f"  {note}" for note in algebra.notes)
    return "\n".join(lines)


def render_principal_sl2(format_name: str) -> str:
    payload = {
        "model": "principal block of category O for sl2",
        "algebra_A": principal_block_algebra().as_dict(),
        "koszul_dual_B": ext_algebra_of_principal_block().as_dict(),
        "double_dual_EB": koszul_dual_of_ext_algebra().as_dict(),
        "self_duality_isomorphism_B_to_EB": vertex_swap_isomorphism(),
        "projective_resolutions": principal_block_resolutions(),
        "koszul_duality_functor_images": principal_block_duality_functor_images(),
    }
    if format_name == "json":
        return json.dumps(payload, indent=2, ensure_ascii=False)

    sections = [
        "Koszul dual computation for the principal block of category O for sl2",
        "",
        render_algebra(principal_block_algebra()),
        "",
        render_algebra(ext_algebra_of_principal_block()),
        "",
        render_algebra(koszul_dual_of_ext_algebra()),
        "",
        "Self-duality isomorphism B -> Bꜝ = E(B):",
    ]
    sections.extend(
        f"  {source} ↦ {target}"
        for source, target in vertex_swap_isomorphism().items()
    )
    sections.append("Projective resolutions:")
    for group, resolutions in principal_block_resolutions().items():
        sections.append(f"  {group}:")
        sections.extend(f"    {resolution}" for resolution in resolutions)
    sections.append("Koszul duality functor images from may26.tex:")
    sections.extend(
        f"  {source} ≅ {target}"
        for source, target in principal_block_duality_functor_images().items()
    )
    return "\n".join(sections)


def render_sl2_functor_images(images: Sequence[dict[str, object]], format_name: str) -> str:
    l_block = sl2_l_block_description(images)
    if format_name == "json":
        return json.dumps(
            {"functor": "K = RHom_O(-, L_block)", "L_block": l_block, "images": images},
            indent=2,
            ensure_ascii=False,
        )

    lines = [
        "Koszul duality functor for regular integral blocks of category O for sl2",
        "K = RHom_O(-, L_block), where L_block is the sum of simples in the relevant block",
        f"Here, L_block = {l_block}",
        "",
    ]
    for image in images:
        lines.append(f"{image['input']} = {image['result']}")
    return "\n".join(lines)


def text_table(headers: Sequence[str], rows: Sequence[Sequence[object]]) -> str:
    if not rows:
        return "  (none)"
    string_rows = [[str(item) for item in row] for row in rows]
    widths = [
        max(len(str(header)), *(len(row[index]) for row in string_rows))
        for index, header in enumerate(headers)
    ]
    lines = [
        "  " + "  ".join(str(header).ljust(widths[index]) for index, header in enumerate(headers)),
        "  " + "  ".join("-" * width for width in widths),
    ]
    lines.extend(
        "  " + "  ".join(row[index].ljust(widths[index]) for index in range(len(headers)))
        for row in string_rows
    )
    return "\n".join(lines)


def render_regular_block_duality(duality: RegularBlockDuality, format_name: str) -> str:
    if format_name == "json":
        return json.dumps(duality.as_dict(), indent=2, ensure_ascii=False)

    rows = [
        (
            entry.projective_label,
            list(entry.projective_word) or "e",
            entry.dual_simple_label,
            list(entry.dual_word) or "e",
        )
        for entry in duality.entries
    ]
    lines = [
        "Koszul duality for a regular integral block",
        f"Input: {duality.input_label}",
        f"Cartan type: {duality.cartan_type}",
        f"Langlands-dual Cartan type: {duality.langlands_dual_type}",
        f"Convention: {duality.convention}",
        "Projective-to-simple correspondence:",
        text_table(["projective", "word", "dual simple", "dual word"], rows),
    ]
    if duality.notes:
        lines.append("Notes:")
        lines.extend(f"  {note}" for note in duality.notes)
    return "\n".join(lines)


def build_sln_parser(
    fixed_n: int | None = None, prog: str = "koszul_duals.py sln"
) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog=prog,
        description=(
            "Find the Koszul-dual regular integral block for sl_n using "
            "fundamental-weight coordinates."
        ),
    )
    if fixed_n is None:
        parser.add_argument("n", type=int, help="The n in sl_n.")
    parser.add_argument(
        "inputs",
        nargs="*",
        help=(
            "Either a highest weight for a block table, e.g. 0 or 1,0, "
            "or object inputs such as P(0,1)[0] and L(-2,2)[1]. "
            "Defaults to the principal regular block."
        ),
    )
    parser.add_argument(
        "--dual-map",
        choices=("left", "right"),
        default="left",
        help="Use w0 w or w w0 for the projective-to-simple vertex correspondence.",
    )
    parser.add_argument(
        "--format", choices=("text", "json"), default="text", help="Output format."
    )
    return parser


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Compute Koszul-dual data for Lie algebra representation examples "
            "from main.tex and may26.tex."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""Examples:
  ./koszul_duals.py sl2-principal
  ./koszul_duals.py sl2
  ./koszul_duals.py sl2 M0 M-2 L0 P-2 M0v
  ./koszul_duals.py sl2-principal --format json
  ./koszul_duals.py regular-block A2
  ./koszul_duals.py regular-block B3
  ./koszul_duals.py sl3
  ./koszul_duals.py sl3 'P(0,1)[0]'
  ./koszul_duals.py sl3 0,1
""",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    sl2_parser = subparsers.add_parser(
        "sl2-principal",
        help="Compute the principal sl2 block algebra, its Ext algebra, and self-duality.",
    )
    sl2_parser.add_argument(
        "--format", choices=("text", "json"), default="text", help="Output format."
    )

    functor_parser = subparsers.add_parser(
        "sl2-functor",
        aliases=("sl2",),
        help="Compute K(M), K(L), K(P), and dominant dual Vermas in regular sl2 blocks.",
    )
    functor_parser.add_argument(
        "objects",
        nargs="*",
        help=(
            "Objects to compute. Use aliases like M0, M-2, L0, P-2, M0v, or quoted M(0)[1]. "
            "Defaults to the list from the notes."
        ),
    )
    functor_parser.add_argument(
        "--format", choices=("text", "json"), default="text", help="Output format."
    )

    regular_parser = subparsers.add_parser(
        "regular-block",
        help="Find the BGS Koszul-dual block for a regular integral finite type block.",
    )
    regular_parser.add_argument("cartan_type", type=category_o.parse_cartan_type)
    regular_parser.add_argument(
        "--dual-map",
        choices=("left", "right"),
        default="left",
        help="Use w0 w or w w0 for the projective-to-simple vertex correspondence.",
    )
    regular_parser.add_argument(
        "--format", choices=("text", "json"), default="text", help="Output format."
    )

    sln_parser = subparsers.add_parser(
        "sln",
        help="Find the Koszul-dual regular integral block for sl_n.",
    )
    sln_parser.add_argument("n", type=int, help="The n in sl_n.")
    sln_parser.add_argument(
        "inputs",
        nargs="*",
        help=(
            "Either a highest weight for a block table, e.g. 0 or 1,0, "
            "or object inputs such as P(0,1)[0] and L(-2,2)[1]. "
            "Defaults to the principal regular block."
        ),
    )
    sln_parser.add_argument(
        "--dual-map",
        choices=("left", "right"),
        default="left",
        help="Use w0 w or w w0 for the projective-to-simple vertex correspondence.",
    )
    sln_parser.add_argument(
        "--format", choices=("text", "json"), default="text", help="Output format."
    )
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == "sl2":
        argv[0] = "sl2-functor"

    if argv:
        sln_match = re.fullmatch(r"sl(\d+)", argv[0])
        if sln_match:
            n = int(sln_match.group(1))
            parser = build_sln_parser(fixed_n=n, prog=f"koszul_duals.py {argv[0]}")
            args = parser.parse_args(argv[1:])
            try:
                output = render_sln_inputs(
                    n, args.inputs, dual_map=args.dual_map, format_name=args.format
                )
            except Exception as exc:
                parser.exit(2, f"error: {exc}\n")
            print(output)
            return 0

    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "sl2-principal":
        print(render_principal_sl2(args.format))
        return 0

    if args.command in {"sl2-functor", "sl2"}:
        try:
            images = sl2_functor_images(args.objects)
        except Exception as exc:
            parser.exit(2, f"error: {exc}\n")
        print(render_sl2_functor_images(images, args.format))
        return 0

    if args.command == "regular-block":
        try:
            duality = regular_block_duality_by_cartan_type(
                args.cartan_type, dual_map=args.dual_map
            )
        except Exception as exc:
            parser.exit(2, f"error: {exc}\n")
        print(render_regular_block_duality(duality, args.format))
        return 0

    if args.command == "sln":
        try:
            output = render_sln_inputs(
                args.n, args.inputs, dual_map=args.dual_map, format_name=args.format
            )
        except Exception as exc:
            parser.exit(2, f"error: {exc}\n")
        print(output)
        return 0

    parser.exit(2, f"error: unknown command {args.command!r}\n")


if __name__ == "__main__":
    raise SystemExit(main())
