import math

import numpy as np
import pytest

from mnasolve import Circuit, ac, dc, parse_netlist, transient
from mnasolve.netlist import parse_value


@pytest.mark.parametrize("text,value", [("4.7k", 4700), ("1meg", 1e6), ("100n", 1e-7), ("10uF", 1e-5), ("2.2", 2.2), ("1e3", 1000)])
def test_parse_value(text, value):
    assert parse_value(text) == pytest.approx(value)


def test_voltage_divider():
    c = Circuit().add("V1", "in", "0", value=12).add("R1", "in", "out", value=2e3).add("R2", "out", "0", value=1e3)
    op = dc(c)
    assert op.v("out") == pytest.approx(4.0)
    assert op.currents["V1"] == pytest.approx(-4e-3)  # SPICE sign convention
    p = op.power(c)
    assert p["V1"] == pytest.approx(-48e-3) and p["R1"] + p["R2"] == pytest.approx(48e-3)


def test_current_source_and_superposition():
    c = parse_netlist("""
    V1 a 0 10
    R1 a b 1k
    R2 b 0 1k
    I1 0 b 2m
    """)
    # Superposition: 5 V from V1 plus 2 mA * (1k || 1k) = 1 V
    assert dc(c).v("b") == pytest.approx(6.0)


def test_wheatstone_bridge_example():
    with open("examples/wheatstone.cir", encoding="utf-8") as fh:
        op = dc(parse_netlist(fh.read()))
    # Thevenin check: open-circuit voltage 10*(1.1/2.1 - 1/2), source resistance
    # (1k || 1k) + (1k || 1.1k), loaded by the 100 ohm galvanometer.
    vth = 10 * (1.1 / 2.1 - 0.5)
    rth = 500 + 1e3 * 1.1e3 / 2.1e3
    assert op.v("b") - op.v("a") == pytest.approx(vth * 100 / (rth + 100), rel=1e-9)


def test_inductor_is_short_and_capacitor_is_open_in_dc():
    c = Circuit().add("V1", "a", "0", value=5).add("L1", "a", "b", value=1e-3).add("R1", "b", "0", value=100).add("C1", "b", "0", value=1e-6)
    op = dc(c)
    assert op.v("b") == pytest.approx(5.0)
    assert op.currents["L1"] == pytest.approx(0.05)


def test_vcvs_inverting_amplifier():
    with open("examples/opamp_inverting.cir", encoding="utf-8") as fh:
        op = dc(parse_netlist(fh.read()))
    assert op.v("out") == pytest.approx(-2.0, rel=1e-3)
    assert abs(op.v("minus")) < 1e-4  # virtual ground


def test_rc_lowpass_ac_minus_3db_at_cutoff():
    c = parse_netlist("V1 in 0 AC 1\nR1 in out 1k\nC1 out 0 100n")
    fc = 1 / (2 * math.pi * 1e3 * 100e-9)
    res = ac(c, [fc / 100, fc, fc * 100])
    mag = np.abs(res["out"])
    assert mag[0] == pytest.approx(1, abs=1e-3)
    assert 20 * math.log10(mag[1]) == pytest.approx(-3.0103, abs=1e-3)
    assert math.degrees(np.angle(res["out"][1])) == pytest.approx(-45, abs=1e-6)
    assert 20 * math.log10(mag[2]) == pytest.approx(-40, abs=0.01)


def test_series_rlc_resonance():
    c = parse_netlist("V1 in 0 AC 1\nL1 in a 10m\nC1 a b 1u\nR1 b 0 10")
    f0 = 1 / (2 * math.pi * math.sqrt(10e-3 * 1e-6))
    res = ac(c, [f0])
    assert abs(res["b"][0]) == pytest.approx(1.0, rel=1e-6)


def test_rc_transient_matches_analytic():
    c = parse_netlist("V1 in 0 5\nR1 in out 1k\nC1 out 0 1u")
    tau = 1e-3
    res = transient(c, 5 * tau, 1e-6)
    t = res["time"]
    exact = 5 * (1 - np.exp(-t / tau))
    assert np.max(np.abs(res["out"] - exact)) < 5e-3
    assert res["out"][-1] == pytest.approx(5 * (1 - math.exp(-5)), abs=5e-3)


def test_rl_transient():
    c = parse_netlist("V1 a 0 10\nR1 a b 10\nL1 b 0 10m")
    res = transient(c, 5e-3, 1e-6)
    exact = 1.0 * (1 - np.exp(-res["time"] / 1e-3))
    assert np.max(np.abs(res["I(L1)"] - exact)) < 5e-3


def test_title_line_is_detected():
    c = parse_netlist("Inverting stage\nV1 a 0 1\nR1 a 0 1k")
    assert c.title == "Inverting stage" and len(c.elements) == 2
    c = parse_netlist("V1 a 0 1\nR1 a 0 1k")
    assert c.title == "" and len(c.elements) == 2


def test_floating_node_is_reported():
    c = Circuit().add("V1", "a", "0", value=1).add("R1", "b", "c", value=1)
    with pytest.raises(ValueError, match="singular"):
        dc(c)


def test_cli_dc(capsys):
    from mnasolve.cli import main

    main(["examples/wheatstone.cir"])
    out = capsys.readouterr().out
    assert "V(a)" in out and "Power absorbed" in out
