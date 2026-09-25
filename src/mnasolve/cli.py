"""Command-line interface.

    mnasolve circuit.cir                 DC operating point
    mnasolve circuit.cir --ac 10 1e6 20  AC sweep, 20 points per decade
    mnasolve circuit.cir --tran 5m 10u   transient to 5 ms with a 10 us step
"""

from __future__ import annotations

import argparse
import math

import numpy as np

from .analysis import ac, dc, transient
from .netlist import parse_netlist, parse_value


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(prog="mnasolve", description="Modified nodal analysis circuit simulator")
    ap.add_argument("netlist")
    ap.add_argument("--ac", nargs=3, metavar=("FSTART", "FSTOP", "PTS_PER_DEC"))
    ap.add_argument("--tran", nargs=2, metavar=("TSTOP", "STEP"))
    ap.add_argument("--node", action="append", help="node(s) to print for AC/transient")
    a = ap.parse_args(argv)

    with open(a.netlist, encoding="utf-8") as fh:
        circuit = parse_netlist(fh.read())
    if circuit.title:
        print(circuit.title)

    if a.ac:
        f0, f1, ppd = parse_value(a.ac[0]), parse_value(a.ac[1]), int(a.ac[2])
        n = max(2, int(round(math.log10(f1 / f0) * ppd)) + 1)
        res = ac(circuit, np.logspace(math.log10(f0), math.log10(f1), n))
        nodes = a.node or circuit.nodes()
        print("freq_hz " + " ".join(f"|V({n})|_dB phase({n})_deg" for n in nodes))
        for j, f in enumerate(res["freq"]):
            cols = []
            for nd in nodes:
                v = res[nd][j]
                cols.append(f"{20 * math.log10(max(abs(v), 1e-30)):9.3f} {math.degrees(np.angle(v)):9.2f}")
            print(f"{f:10.4g} " + " ".join(cols))
    elif a.tran:
        res = transient(circuit, parse_value(a.tran[0]), parse_value(a.tran[1]))
        nodes = a.node or circuit.nodes()
        print("time_s " + " ".join(f"V({n})" for n in nodes))
        stride = max(1, len(res["time"]) // 50)
        for k in range(0, len(res["time"]), stride):
            print(f"{res['time'][k]:10.4g} " + " ".join(f"{res[n][k]:10.5f}" for n in nodes))
    else:
        op = dc(circuit)
        print("Node voltages")
        for n, v in op.voltages.items():
            print(f"  V({n}) = {v:.6g} V")
        if op.currents:
            print("Branch currents")
            for n, i in op.currents.items():
                print(f"  I({n}) = {i:.6g} A")
        print("Power absorbed")
        for n, p in op.power(circuit).items():
            print(f"  {n}: {p:.6g} W")


if __name__ == "__main__":
    main()
