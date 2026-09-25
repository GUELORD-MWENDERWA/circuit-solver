"""Circuit description and SPICE-style netlist parser.

Supported elements (node 0 is ground):

    Rname n+ n- value            resistor
    Cname n+ n- value            capacitor
    Lname n+ n- value            inductor
    Vname n+ n- [DC] value [AC mag [phase]]   independent voltage source
    Iname n+ n- [DC] value [AC mag [phase]]   independent current source (flows n+ -> n- inside the source)
    Ename n+ n- nc+ nc- gain     voltage-controlled voltage source

Values accept SPICE suffixes: f p n u m k meg g t (case-insensitive), e.g. 4.7k, 100n, 1meg.
Lines starting with '*' are comments; '.end' stops parsing.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

_SUFFIX = {"t": 1e12, "g": 1e9, "meg": 1e6, "k": 1e3, "m": 1e-3, "u": 1e-6, "n": 1e-9, "p": 1e-12, "f": 1e-15}
_NUM = re.compile(r"^([-+]?\d*\.?\d+(?:e[-+]?\d+)?)(meg|[tgkmunpf])?[a-z]*$", re.IGNORECASE)


def parse_value(text: str) -> float:
    m = _NUM.match(text.strip())
    if not m:
        raise ValueError(f"invalid value {text!r}")
    number, suffix = m.groups()
    return float(number) * (_SUFFIX[suffix.lower()] if suffix else 1.0)


@dataclass
class Element:
    name: str
    kind: str                 # R C L V I E
    nodes: tuple[str, ...]
    value: float = 0.0        # R, C, L value; DC value for sources; gain for E
    ac_mag: float = 0.0
    ac_phase: float = 0.0     # degrees


@dataclass
class Circuit:
    elements: list[Element] = field(default_factory=list)
    title: str = ""

    def add(self, name: str, *nodes, value: float = 0.0, ac: float = 0.0, phase: float = 0.0) -> "Circuit":
        kind = name[0].upper()
        if kind not in "RCLVIE":
            raise ValueError(f"unsupported element {name}")
        expected = 4 if kind == "E" else 2
        if len(nodes) != expected:
            raise ValueError(f"{name} needs {expected} nodes")
        if kind in "RCL" and value <= 0:
            raise ValueError(f"{name} must have a positive value")
        if any(e.name.upper() == name.upper() for e in self.elements):
            raise ValueError(f"duplicate element {name}")
        self.elements.append(Element(name, kind, tuple(str(n) for n in nodes), value, ac, phase))
        return self

    def nodes(self) -> list[str]:
        """Non-ground node names in order of first appearance."""
        seen: list[str] = []
        for e in self.elements:
            for n in e.nodes:
                if n != "0" and n.lower() != "gnd" and n not in seen:
                    seen.append(n)
        return seen

    def branch_elements(self) -> list[Element]:
        """Elements that add a current unknown to the MNA system."""
        return [e for e in self.elements if e.kind in "VEL"]


def _parse_element(c: Circuit, line: str) -> None:
    tok = line.split()
    name, kind = tok[0], tok[0][0].upper()
    if kind == "E":
        if len(tok) != 6:
            raise ValueError(f"{name}: expected 'E n+ n- nc+ nc- gain'")
        c.add(name, *tok[1:5], value=parse_value(tok[5]))
    elif kind in "RCL":
        if len(tok) != 4:
            raise ValueError(f"{name}: expected '{kind} n+ n- value'")
        c.add(name, tok[1], tok[2], value=parse_value(tok[3]))
    elif kind in "VI":
        if len(tok) < 4:
            raise ValueError(f"{name}: missing value")
        rest = [t.lower() for t in tok[3:]]
        dc_val, ac_mag, ac_ph = 0.0, 0.0, 0.0
        j = 0
        while j < len(rest):
            if rest[j] == "dc":
                dc_val = parse_value(rest[j + 1])
                j += 2
            elif rest[j] == "ac":
                ac_mag = parse_value(rest[j + 1])
                j += 2
                if j < len(rest) and _NUM.match(rest[j]):
                    ac_ph = parse_value(rest[j])
                    j += 1
            else:
                dc_val = parse_value(rest[j])
                j += 1
        c.add(name, tok[1], tok[2], value=dc_val, ac=ac_mag, phase=ac_ph)
    else:
        raise ValueError(f"unsupported element {name}")


def parse_netlist(text: str) -> Circuit:
    """Parse a netlist.

    As in SPICE, the first line may be a title. It is treated as a title when it
    does not parse as a valid element line, so short netlists without a title work too.
    """
    c = Circuit()
    lines = text.strip().splitlines()
    for i, raw in enumerate(lines):
        line = raw.split(";")[0].strip()
        if not line or line.startswith("*"):
            continue
        if line.lower().startswith(".end"):
            break
        if line.startswith("."):
            continue  # analysis cards are handled by the CLI
        if i == 0:
            probe = Circuit()
            try:
                _parse_element(probe, line)
            except (ValueError, IndexError):
                c.title = line
                continue
        try:
            _parse_element(c, line)
        except (ValueError, IndexError) as exc:
            raise ValueError(f"line {i + 1}: {exc}") from exc
    return c
