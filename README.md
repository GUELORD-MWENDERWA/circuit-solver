# mnasolve: Circuit Simulation with Modified Nodal Analysis

![tests](https://github.com/GUELORD-MWENDERWA/circuit-solver/actions/workflows/tests.yml/badge.svg)
![python](https://img.shields.io/badge/python-3.10%2B-3776AB)
![license](https://img.shields.io/badge/license-MIT-green)

A compact SPICE-like circuit simulator written in Python. It reads a standard netlist, builds the modified nodal analysis (MNA) equations by stamping each component, and solves them for the DC operating point, a small-signal AC frequency sweep, or a time-domain transient.

The project connects two subjects: circuit theory from general electronics (Kirchhoff's laws, Thevenin equivalents, RC and RLC responses) and numerical methods from computer science (sparse linear systems, complex arithmetic, implicit integration of differential equations).

## Capabilities

| Analysis | Method | Output |
| --- | --- | --- |
| DC operating point | Solve `G x = b`; capacitors open, inductors shorted | Node voltages, source and inductor currents, power per element |
| AC sweep | Solve `(G + jωC) x = b` for each frequency | Complex node voltages (magnitude in dB, phase in degrees) |
| Transient | Backward Euler: `(G + C/h) x[n+1] = b + (C/h) x[n]` | Node voltages and branch currents over time |

**Elements:** resistor `R`, capacitor `C`, inductor `L`, independent voltage source `V` and current source `I` (DC and AC values), voltage-controlled voltage source `E` (used to model an op-amp).

**Diagnostics:** floating nodes and loops of voltage sources are reported as a singular MNA matrix with an explicit message.

## Installation

```bash
git clone https://github.com/GUELORD-MWENDERWA/circuit-solver.git
cd circuit-solver
pip install -e ".[dev]"
```

## Netlist format

```spice
First-order RC low-pass, fc = 1.59 kHz
V1 in 0 DC 5 AC 1
R1 in out 1k
C1 out 0 100n
.end
```

Node `0` (or `gnd`) is ground. Values accept SPICE suffixes: `f p n u m k meg g t`. The first line is a title when it is not a valid element line; `*` starts a comment.

## Examples

**DC operating point** of an unbalanced Wheatstone bridge ([`examples/wheatstone.cir`](examples/wheatstone.cir)):

```console
$ mnasolve examples/wheatstone.cir
Node voltages
  V(top) = 10 V
  V(a) = 5.10593 V
  V(b) = 5.12712 V
Branch currents
  I(V1) = -0.00976695 A
Power absorbed
  V1: -0.0976695 W
  ...
```

The 21.2 mV across the galvanometer is exactly the Thevenin prediction; the test suite checks it.

**AC sweep** of the RC low-pass filter:

```console
$ mnasolve examples/rc_lowpass.cir --ac 100 100k 1 --node out
freq_hz |V(out)|_dB phase(out)_deg
       100    -0.017     -3.60
      1000    -1.445    -32.14
     1e+04   -16.072    -80.96
     1e+05   -35.965    -89.09
```

**Op-amp stage** modelled as a VCVS with an open-loop gain of 100 000 ([`examples/opamp_inverting.cir`](examples/opamp_inverting.cir)): the output settles at -1.99978 V for a 0.2 V input, which shows the finite-gain error of an otherwise ideal gain of -10, and the inverting input sits at 20 µV (virtual ground).

**Transient**, from Python:

```python
from mnasolve import parse_netlist, transient

c = parse_netlist("V1 in 0 5\nR1 in out 1k\nC1 out 0 1u")
res = transient(c, t_stop=5e-3, h=1e-6)
res["time"], res["out"]        # capacitor charging curve, tau = 1 ms
```

## How it works

Each unknown is either a node voltage or the current through an element that constrains a voltage (voltage sources, VCVS, inductors). Components add their contribution to the system matrix independently:

- A conductance `1/R` between nodes `a` and `b` adds `+g` on the diagonal entries `(a,a)`, `(b,b)` and `-g` on `(a,b)`, `(b,a)`.
- A voltage source adds a row and a column coupling its current to the two node voltages, and its value to the right-hand side.
- Capacitors and inductors stamp into a second matrix `C`, which becomes `jωC` in AC analysis and `C/h` in the backward Euler transient.

Backward Euler is only first-order accurate but unconditionally stable, which suits stiff circuits with very different time constants.

## Testing

```bash
pytest
```

Every analysis is compared with a closed-form result: divider and superposition in DC, the -3 dB and -45 degree point of an RC filter, unity transfer at series RLC resonance, the exponential charge of RC and RL circuits, and the virtual ground of an op-amp stage.

## Limitations

- Linear elements only: no diodes or transistors yet (they need Newton-Raphson iteration around the MNA solve).
- Dense matrices: fine for circuits of a few hundred nodes.
- Transient sources are constant steps from t = 0.

## Roadmap

- Diode model with Newton-Raphson for non-linear DC analysis
- Pulse and sine sources, trapezoidal integration
- Sparse solver and `.op`, `.ac`, `.tran` cards in the netlist

## License

MIT. See [LICENSE](LICENSE).
