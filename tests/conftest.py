import pytest
import sensitivity_calculator.sensitivity as sens
import sys


@pytest.fixture(scope="session")
def broadband_output():
    i = sens.getInputs("src/sensitivity_calculator/input.yaml")
    angle = 90 - i["observationElevationAngle"]
    coldSpillOverEfficiency = sens.getSpillEfficiency(i)

    calculate = sens.calcByAngle(i["diameter"], i["t"], i["wfe"], i["eta"], i["doe"], i["t_int"], i["pixelYield"], i["szCamNumPoln"], i["eorSpecNumPoln"],
                                 i["t_filter_cold"], i["t_lens_cold"], i["t_uhdpe_window"], coldSpillOverEfficiency, i["singleModedAOmegaLambda2"],
                                 i["spatialPixels"], i["fpi"], i["eqbw"], i["centerFrequency"], i["detectorNEP"],
                                 i["backgroundSubtractionDegradationFactor"], i["sensitivity"], i["hoursPerYear"], i["sensPerBeam"], i["r"], i["signal"])

    outputs = calculate(angle)

    valueDisplay = sens.valDisplayPartial(
        i["outputFreq"], i["centerFrequency"], outputs["wavelength"], i["decimalPlaces"])
    quartileDisplay = sens.quartDisplayPartial(
        i["outputFreq"], i["centerFrequency"], outputs["wavelength"], i["decimalPlaces"])

    return i, coldSpillOverEfficiency, calculate, outputs, valueDisplay, quartileDisplay


@pytest.fixture
def capture_stdout(monkeypatch):
    buffer = {"stdout": "", "write_calls": 0}

    def fake_write(s):
        buffer["stdout"] += s
        buffer["write_calls"] += 1

    monkeypatch.setattr(sys.stdout, 'write', fake_write)
    return buffer
