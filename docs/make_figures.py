"""Regenerate the figures in docs/images by simulating the example netlists.

    pip install -e . matplotlib
    python docs/make_figures.py
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from mnasolve import ac, parse_netlist, transient

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs" / "images"
plt.rcParams.update({"figure.dpi": 150, "axes.grid": True, "grid.alpha": 0.3, "axes.spines.top": False, "axes.spines.right": False})

RC = parse_netlist((ROOT / "examples" / "rc_lowpass.cir").read_text())
RLC = parse_netlist("""Series RLC band-pass, f0 = 5.03 kHz
V1 in 0 DC 1 AC 1
L1 in a 10m
C1 a out 100n
R1 out 0 100
.end""")


def bode() -> None:
    f = np.logspace(1, 6, 400)
    fig, (am, ap) = plt.subplots(2, 1, figsize=(8, 5.5), sharex=True)
    for name, circuit in (("RC low-pass, 1 k / 100 n", RC), ("RLC band-pass, 10 mH / 100 nF / 100 R", RLC)):
        h = ac(circuit, f)["out"]
        am.semilogx(f, 20 * np.log10(np.abs(h)), label=name)
        ap.semilogx(f, np.degrees(np.angle(h)))
    fc = 1 / (2 * np.pi * 1e3 * 100e-9)
    am.axvline(fc, color="gray", ls=":")
    am.annotate(f"fc = {fc:.0f} Hz, -3 dB", (fc, -3), xytext=(fc / 40, -12), arrowprops={"arrowstyle": "->"})
    am.set(ylabel="gain (dB)", ylim=(-60, 5), title="AC sweep computed by mnasolve")
    am.legend()
    ap.set(xlabel="frequency (Hz)", ylabel="phase (degrees)")
    fig.tight_layout()
    fig.savefig(OUT / "bode.png")


def step_response() -> None:
    tau = 1e3 * 100e-9
    r = transient(RC, t_stop=6 * tau, h=tau / 50)
    t = r["time"]
    fig, ax = plt.subplots(figsize=(8, 3.8))
    ax.plot(t * 1e3, r["out"], lw=2, label="mnasolve, backward Euler, h = tau / 50")
    ax.plot(t * 1e3, 5 * (1 - np.exp(-t / tau)), "k--", lw=1, label="analytic 5 (1 - exp(-t / RC))")
    ax.axhline(5 * (1 - np.exp(-1)), color="gray", ls=":")
    ax.set(xlabel="time (ms)", ylabel="V(out)", title="RC charging from a 5 V step, tau = 100 us")
    ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / "rc_step.png")


if __name__ == "__main__":
    OUT.mkdir(exist_ok=True)
    bode()
    step_response()
