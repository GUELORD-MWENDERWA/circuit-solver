"""Modified nodal analysis.

The unknown vector is [node voltages..., branch currents...]. Voltage sources,
VCVS and inductors add a branch current unknown. Each element "stamps" its
contribution into the matrix G (and C for reactive parts) so that

    DC:        G x = b
    AC:        (G + j w C) x = b_ac
    transient: (G + C/h) x_{n+1} = b + (C/h) x_n        (backward Euler)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .netlist import Circuit


def _is_ground(n: str) -> bool:
    return n == "0" or n.lower() == "gnd"


class _System:
    def __init__(self, circuit: Circuit):
        self.circuit = circuit
        self.node_names = circuit.nodes()
        self.idx = {n: i for i, n in enumerate(self.node_names)}
        branches = circuit.branch_elements()
        self.branch = {e.name: len(self.node_names) + k for k, e in enumerate(branches)}
        self.size = len(self.node_names) + len(branches)
        if self.size == 0:
            raise ValueError("empty circuit")

    def i(self, node: str) -> int | None:
        return None if _is_ground(node) else self.idx[node]

    def build(self, mode: str, omega: float = 0.0):
        """Return (G, C, b) for mode 'dc', 'ac' or 'tran'."""
        n = self.size
        dtype = complex if mode == "ac" else float
        g = np.zeros((n, n), dtype=dtype)
        c = np.zeros((n, n), dtype=dtype)
        b = np.zeros(n, dtype=dtype)

        def stamp2(m, a, bnode, val):
            ia, ib = self.i(a), self.i(bnode)
            if ia is not None:
                m[ia, ia] += val
            if ib is not None:
                m[ib, ib] += val
            if ia is not None and ib is not None:
                m[ia, ib] -= val
                m[ib, ia] -= val

        def source_value(e):
            if mode == "ac":
                return e.ac_mag * np.exp(1j * np.deg2rad(e.ac_phase))
            return e.value

        for e in self.circuit.elements:
            if e.kind == "R":
                stamp2(g, *e.nodes, 1.0 / e.value)
            elif e.kind == "C":
                if mode != "dc":
                    stamp2(c, *e.nodes, e.value)  # open circuit in DC
            elif e.kind in "VL":
                k = self.branch[e.name]
                p, q = self.i(e.nodes[0]), self.i(e.nodes[1])
                if p is not None:
                    g[p, k] += 1
                    g[k, p] += 1
                if q is not None:
                    g[q, k] -= 1
                    g[k, q] -= 1
                if e.kind == "V":
                    b[k] = source_value(e)
                elif mode != "dc":
                    c[k, k] -= e.value  # v = L di/dt ; in DC the inductor is a short (row stays v+ - v- = 0)
            elif e.kind == "I":
                val = source_value(e)
                p, q = self.i(e.nodes[0]), self.i(e.nodes[1])
                if p is not None:
                    b[p] -= val
                if q is not None:
                    b[q] += val
            elif e.kind == "E":
                k = self.branch[e.name]
                p, q, cp, cq = (self.i(x) for x in e.nodes)
                if p is not None:
                    g[p, k] += 1
                    g[k, p] += 1
                if q is not None:
                    g[q, k] -= 1
                    g[k, q] -= 1
                if cp is not None:
                    g[k, cp] -= e.value
                if cq is not None:
                    g[k, cq] += e.value
        return g, c, b

    def solve(self, a, b):
        try:
            return np.linalg.solve(a, b)
        except np.linalg.LinAlgError as exc:
            raise ValueError("singular MNA matrix: check for floating nodes or loops of voltage sources") from exc


@dataclass
class DCResult:
    voltages: dict[str, float]
    currents: dict[str, float]

    def v(self, node: str) -> float:
        return 0.0 if _is_ground(node) else self.voltages[node]

    def power(self, circuit: Circuit) -> dict[str, float]:
        """Power absorbed by each element (negative means delivered)."""
        out = {}
        for e in circuit.elements:
            vd = self.v(e.nodes[0]) - self.v(e.nodes[1])
            if e.kind == "R":
                out[e.name] = vd * vd / e.value
            elif e.kind in "VEL":
                out[e.name] = vd * self.currents[e.name]
            elif e.kind == "I":
                out[e.name] = -vd * e.value
        return out


def dc(circuit: Circuit) -> DCResult:
    """DC operating point. Capacitors are open, inductors are short circuits.

    Branch currents follow the SPICE convention: positive current flows into the
    positive terminal of the element (so a source delivering power has a negative current).
    """
    s = _System(circuit)
    g, _, b = s.build("dc")
    x = s.solve(g, b)
    volts = {n: float(x[i]) for n, i in s.idx.items()}
    amps = {name: float(x[k]) for name, k in s.branch.items()}
    return DCResult(volts, amps)


def ac(circuit: Circuit, frequencies) -> dict[str, np.ndarray]:
    """Small-signal AC sweep. Returns complex node voltages for each frequency."""
    s = _System(circuit)
    freqs = np.atleast_1d(np.asarray(frequencies, dtype=float))
    out = {n: np.zeros(len(freqs), dtype=complex) for n in s.node_names}
    for j, f in enumerate(freqs):
        w = 2 * np.pi * f
        g, c, b = s.build("ac", w)
        x = s.solve(g + 1j * w * c, b)
        for n, i in s.idx.items():
            out[n][j] = x[i]
    out["freq"] = freqs
    return out


def transient(circuit: Circuit, t_stop: float, h: float, initial: dict[str, float] | None = None) -> dict[str, np.ndarray]:
    """Transient analysis with backward Euler (A-stable, first order).

    Sources are constant (their DC value) from t = 0. Unless `initial` gives node
    voltages, all capacitor voltages and inductor currents start at zero.
    """
    s = _System(circuit)
    g, c, b = s.build("tran")
    steps = int(round(t_stop / h))
    x = np.zeros(s.size)
    for n, v in (initial or {}).items():
        x[s.idx[n]] = v
    a = g + c / h
    t = np.linspace(0, steps * h, steps + 1)
    hist = np.zeros((steps + 1, s.size))
    hist[0] = x
    for k in range(1, steps + 1):
        x = s.solve(a, b + c @ x / h)
        hist[k] = x
    result = {n: hist[:, i] for n, i in s.idx.items()}
    result.update({f"I({name})": hist[:, k] for name, k in s.branch.items()})
    result["time"] = t
    return result
