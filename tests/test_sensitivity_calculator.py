import sensitivity_calculator.sensitivity as sens
import numpy as np


def test_broadband_sensitivity_file(broadband_output):
    i, coldSpillOverEfficiency, calculate, outputs, valueDisplay, quartileDisplay = broadband_output
    sens.outputSensitivityFile(outputs, valueDisplay, quartileDisplay)
    assert True


def test_broadband_power_file(broadband_output):
    i, coldSpillOverEfficiency, calculate, outputs, valueDisplay, quartileDisplay = broadband_output
    sens.outputPowerFile(i, outputs, calculate, quartileDisplay)
    assert True


def test_broadband_spill_efficiency_file(broadband_output):
    i, coldSpillOverEfficiency, calculate, outputs, valueDisplay, quartileDisplay = broadband_output
    sens.outputSpillEfficiencyFile(i, calculate, coldSpillOverEfficiency)
    assert True


def test_broadband_loadings(broadband_output):
    i, coldSpillOverEfficiency, calculate, outputs, valueDisplay, quartileDisplay = broadband_output
    sens.outputLoadings(i, calculate)
    assert True


def test_broadband_noise_curves(broadband_output):
    i, coldSpillOverEfficiency, calculate, outputs, valueDisplay, quartileDisplay = broadband_output
    noiseCurves = sens.getNoiseCurves(i, outputs)
    ell, N_ell_T_full, N_ell_P_full = noiseCurves
    assert np.all(N_ell_T_full > 0)
    assert np.all(N_ell_P_full > 0)


def test_eor_noise_curves(capture_stdout):
    inputs = {'diameter': 5.7, 't': 273, 'wfe': 10.7, 'eta': 0.98, 'doe': 0.8, 'pixelYield': 0.8,
              'eorSpecNumPoln': 2, 't_filter_cold': np.array([1, 1]), 't_lens_cold': np.array([.98, .98]), 't_uhdpe_window': np.array([1, 1]), 'spatialPixels': np.array([3456, 3072]),
              'centerFrequency': np.array([262.5*10**9, 367.5*10**9]), 'detectorNEP': 0,
              'backgroundSubtractionDegradationFactor': 1, 'observationElevationAngle': 45, 'detectorSpacing': np.array([2.75, 2.09]), 'lyotStopAngle': 13.4}
    rfpairs = np.array([(101, 250*10**9), (102, 350*10**9),
                       (103, 275*10**9), (104, 100*10**9)])
    noiseCurves = sens.eorNoiseCurves(inputs, rfpairs)
    assert np.all(noiseCurves[(101, 250*10**9)][1] > 0)
    assert np.all(noiseCurves[(101, 250*10**9)][2] > 0)
    assert np.all(noiseCurves[(102, 350*10**9)][1] > 0)
    assert np.all(noiseCurves[(102, 350*10**9)][2] > 0)
    assert np.all(noiseCurves[(103, 275*10**9)][1] > 0)
    assert np.all(noiseCurves[(103, 275*10**9)][2] > 0)
    assert capture_stdout["stdout"] == "Skipped 104, 100000000000 for not being in any frequency ranges\n"
