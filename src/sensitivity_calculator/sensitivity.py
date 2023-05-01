import os
import yaml
import numpy as np
from functools import partial
import matplotlib.pyplot as plt
from texttable import Texttable
import sensitivity_calculator.noise as noise_file
import scipy.optimize as op
import sensitivity_calculator.pwvCalculator as ACTPWV
from matplotlib import rc
import mapsims
#from sensitivity_calculator.ad_fns import *
import healpy as hp
import warnings
from pathlib import Path
warnings.filterwarnings("ignore")
absolute_path = os.path.dirname(__file__)

fullData = True

# Values taken from https://arxiv.org/pdf/2007.04262.pdf
ccatPWVQ1 = 0.36
ccatPWVQ2 = 0.67
ccatPWVQ3 = 1.28

# These functions work on numpy arrays (componentwise) and help with code clarity
_pi = np.pi
_ln = np.log
_sqrt = np.sqrt
_e = np.e
_h = 6.62606957e-34
_k = 1.3806488e-23
_c = 299792458
_arrayify = np.ndarray.tolist


def getInputs(filePath):
    """Returns a dictionary of variables corresponding to the yaml file at [filePath]. The key of each entry is the same as in the yaml file."""
    # Declare all variables with line breaks repesenting section breaks in the excel file this code is based on
    diameter = None
    t = None
    wfe = None

    eta = None

    doe = None
    t_int = None
    pixelYield = None
    szCamNumPoln = None
    eorSpecNumPoln = None

    t_filter_cold = None
    t_lens_cold = None
    t_uhdpe_window = None
    coldSpillOverEfficiency = None
    singleModedAOmegaLambda2 = None
    spatialPixels = None
    fpi = None

    eqbw = None
    centerFrequency = None
    detectorNEP = None
    backgroundSubtractionDegradationFactor = None
    sensitivity = None
    hoursPerYear = None
    sensPerBeam = None

    r = None
    signal = None

    decimalPlaces = None
    observationElevationAngle = None
    outputFreq = None
    detectorSpacing = None
    lyotStopAngle = None

    # Read in yaml file and assign all input data to their respective variables
    stream = open(filePath, 'r')
    dictionary = yaml.safe_load(stream)
    for key, value in dictionary.items():
        if key == "diameter":
            diameter = value
        if key == "t":
            t = value
        if key == "wfe":
            wfe = value

        if key == "eta":
            eta = value

        if key == "doe":
            doe = value
        if key == "t_int":
            t_int = value
        if key == "pixelYield":
            pixelYield = value
        if key == "sz-camNumPoln":
            szCamNumPoln = value
        if key == "eor-specNumPoln":
            eorSpecNumPoln = value

        if key == "t_filter_cold":
            t_filter_cold = np.array(value)
        if key == "t_lens_cold":
            t_lens_cold = np.array(value)
        if key == "t_uhdpe_window":
            t_uhdpe_window = np.array(value)
        if key == "coldSpillOverEfficiency":
            coldSpillOverEfficiency = np.array(value)
        if key == "singleModedAOmegaLambda2":
            singleModedAOmegaLambda2 = np.array(value)
        if key == "spatialPixels":
            spatialPixels = np.array(value)
        if key == "fpi":
            fpi = np.array(value)

        if key == "eqbw":
            eqbw = np.array(value)
        if key == "centerFrequency":
            centerFrequency = np.array(value)
        if key == "detectorNEP":
            detectorNEP = value
        if key == "backgroundSubtractionDegradationFactor":
            backgroundSubtractionDegradationFactor = value
        if key == "sensitivity":
            sensitivity = value
        if key == "hoursPerYear":
            hoursPerYear = value
        if key == "sensPerBeam":
            sensPerBeam = value

        if key == "r":
            r = np.array(value)
        if key == "signal":
            signal = np.array(value)

        if key == "decimalPlaces":
            decimalPlaces = value
        if key == "observationElevationAngle":
            observationElevationAngle = value
        if key == "outputFreq":
            outputFreq = value
        if key == "detectorSpacing":
            detectorSpacing = np.array(value)
        if key == "lyotStopAngle":
            lyotStopAngle = value

    stream.close()
    return {"diameter": diameter, "t": t, "wfe": wfe, "eta": eta, "doe": doe, "t_int": t_int, "pixelYield": pixelYield, "szCamNumPoln": szCamNumPoln, "eorSpecNumPoln": eorSpecNumPoln, "t_filter_cold": t_filter_cold, "t_lens_cold": t_lens_cold, "t_uhdpe_window": t_uhdpe_window, "coldSpillOverEfficiency": coldSpillOverEfficiency, "singleModedAOmegaLambda2": singleModedAOmegaLambda2, "spatialPixels": spatialPixels, "fpi": fpi, "eqbw": eqbw, "centerFrequency": centerFrequency, "detectorNEP": detectorNEP, "backgroundSubtractionDegradationFactor": backgroundSubtractionDegradationFactor, "sensitivity": sensitivity, "hoursPerYear": hoursPerYear, "sensPerBeam": sensPerBeam, "r": r, "signal": signal, "decimalPlaces": decimalPlaces, "observationElevationAngle": observationElevationAngle, "outputFreq": outputFreq, "detectorSpacing": detectorSpacing, "lyotStopAngle": lyotStopAngle}


def _a_to_CMB(f):
    kb = 1.3806488e-23
    T = 2.725
    v = f*1e9
    x = _h*v/(kb*T)
    return 1./(x**2*np.exp(x)/(np.exp(x)-1)**2)


def _calculate(diameter, t, wfe, eta, doe, t_int, pixelYield, szCamNumPoln, eorSpecNumPoln, t_filter_cold, t_lens_cold, t_uhdpe_window, coldSpillOverEfficiency, singleModedAOmegaLambda2, spatialPixels, eqbw, centerFrequency, detectorNEP, backgroundSubtractionDegradationFactor, r, eqtrans):
    wavelength = _c/centerFrequency*10**6

    aToCMB = np.array([_a_to_CMB(i/1e9) for i in centerFrequency])

    # Telescope
    a = _pi*(diameter/2)**2

    # Detector
    szCamCorForPoln = abs(3-szCamNumPoln)
    eorSpecCorForPoln = abs(3-eorSpecNumPoln)

    # Instrument, throughput
    t_cold = t_filter_cold*t_lens_cold**3
    e_window_warm = 1 - t_uhdpe_window

    # Instrument, beams/area of focal plane
    beam = 1.2*wavelength/diameter/1000000*206265
    print(beam)
    solidAngle = _pi/4/_ln(2)*(beam/206264)**2

    # Weather quartiles/broadband
    e_warm = (t_uhdpe_window[:, None])*((1-eqtrans)
                                        * eta+(1-eta))+(e_window_warm[:, None])
    ruze = 1/(_e**((4*_pi*wfe/(wavelength))**2))
    weather_t_cold = (t_cold)
    occupancy = np.ones((len(wavelength), 3)) * 1 / \
        (_e**(_h*_c/((wavelength)*10**(-6)*_k*t))-1)[:, None]
    acceptedModes = szCamNumPoln * \
        (coldSpillOverEfficiency)*(singleModedAOmegaLambda2)
    powerPerPixel = _h*_c*e_warm*(weather_t_cold[:, None])*(acceptedModes[:, None])*occupancy*(eqbw[:, None])*doe/(
        ((wavelength)*10**(-6))[:, None])  # Differs from excel sheet becauses of correctly using t(cold)
    photonNoiseNEP = _h*_c/((wavelength)[:, None]*10**(-6))*(1/t_int*acceptedModes[:, None]*eqbw[:, None]
                                                             * e_warm*weather_t_cold[:, None]*doe*occupancy*(1+e_warm*weather_t_cold[:, None]*doe*occupancy))**0.5
    totalNEP = _sqrt(photonNoiseNEP**2 + detectorNEP**2)
    coldTerminatedSpillover = (coldSpillOverEfficiency)
    nef = szCamCorForPoln*totalNEP/(a*weather_t_cold[:, None]*doe*eta*ruze[:, None]
                                    * eqtrans*coldTerminatedSpillover[:, None]*backgroundSubtractionDegradationFactor)
    nefd = nef*10**26*1000/(eqbw[:, None])
    net = nefd*10**(-29)/((solidAngle)[:, None])*((wavelength)
                                                  [:, None]*0.000001)**2/(2*1.38*10**(-23))*1000

    # Unlabeled section/final noise equivalent calculations
    arrayNETRJ = (net/((spatialPixels)[:, None]*pixelYield)**0.5)
    arrayNETCMB = arrayNETRJ*aToCMB[:, None]*1000
    arrayNEI = nefd/1000/((solidAngle)[:, None]) / \
        ((spatialPixels)[:, None]*pixelYield)**0.5
    netW8Avg = _sqrt(3/(1/arrayNETCMB[:, 0]**2+1/arrayNETCMB[:, 1]
                        ** 2+1/arrayNETCMB[:, 2]**2).astype(float))
    netW8RJ = netW8Avg/aToCMB
    neiW8 = _sqrt(3/(1/arrayNEI[:, 0]**2+1/arrayNEI[:, 1]
                     ** 2+1/arrayNEI[:, 2]**2).astype(float))
    nefd = arrayNEI*(solidAngle)[:, None] * \
        _sqrt((spatialPixels)[:, None]*pixelYield)*1000

    # EoR Spec
    eorEqBw = _c/(wavelength[:, None]*10**(-6)*r)*_pi/2
    eorEqTrans = eqtrans
    eorE_warm = (t_uhdpe_window)[:, None]*((1-eorEqTrans)
                                           * eta+(1-eta))+(e_window_warm)[:, None]
    eorRuze = 1/_e**((4*_pi*wfe/(wavelength))**2)
    eorT_cold = (t_cold)*0.9  # Always takes from 739 um in the sheet
    eorOccupancy = 1/(_e**(_h*_c/((wavelength)*10**(-6)*_k*t))-1)
    eorAcceptedModes = eorSpecNumPoln*a*(coldSpillOverEfficiency)*(solidAngle)/(
        (wavelength)*10**(-6))**2  # In sheet 1071 um is calculated incorrectly
    eorPhotonNoiseNEP = _h*_c/((wavelength)[:, None]*10**(-6))*(eorAcceptedModes[:, None]*eorEqBw*eorE_warm *
                                                                eorT_cold[:, None]*doe*eorOccupancy[:, None]*(1+eorE_warm*eorT_cold[:, None]*doe*eorOccupancy[:, None]))**0.5
    eorTotalNEP = _sqrt(eorPhotonNoiseNEP**2 + detectorNEP**2)
    eorColdTerminatedSpillover = (coldSpillOverEfficiency)
    eorNEF = eorSpecCorForPoln*eorTotalNEP/(a*eorT_cold[:, None]*doe*eta*eorRuze[:, None]
                                            * eorEqTrans*eorColdTerminatedSpillover[:, None]*backgroundSubtractionDegradationFactor)
    eorNEFD = eorNEF/eorEqBw*10**26*1000
    eorNEI = eorNEFD/(solidAngle)[:, None]/1000
    eorPowerPerPixel = _h*_c/(wavelength[:, None]*10**(-6))*eorE_warm*eorT_cold[:, None]*eorAcceptedModes[:,
                                                                                                          None]*eorOccupancy[:, None]*eorEqBw*doe  # Look at t_cold in excel sheet

    return {"netW8Avg": netW8Avg, "netW8RJ":
            netW8RJ, "neiW8": neiW8, "eorNEFD": eorNEFD, "eorNEI": eorNEI, "powerPerPixel":
            powerPerPixel, "eorPowerPerPixel": eorPowerPerPixel, "wavelength": wavelength, "beam": beam}


def _pwvToTick(pwv):
    return int(pwv/(ccatPWVQ2/20))


def _averageTransSE(filePath, start, end, col=1):
    file = open(filePath + ".out", "r")
    data = file.readlines()
    x = [float(i.split(" ")[0]) for i in data]
    y = [float(i.split(" ")[col]) for i in data]
    file.close()

    transmission = 0
    for (f, trans) in zip(x, y):
        if f < start:
            continue
        if f >= end:
            break
        transmission += trans

    return transmission / ((end - start) * 100)


def _averageTransHelper(filePath, center, width, col=1):
    return _averageTransSE(filePath, (center-width/2)/1e9, (center+width/2)/1e9, col)


def _averageTrans(prefix, angle, percentile, center, width, col=1, newFilePathFormat=False):
    filePath = f"{prefix}{percentile}/ACT_annual_{percentile}."
    if newFilePathFormat:
        filePath = f"{prefix}ACT_annual_{percentile}."
    if fullData:
        if angle >= 15 and angle <= 75 and int(angle) == angle:
            return _averageTransHelper(filePath + str(angle), center, width, col)
        elif int(angle) != angle:
            floor = _averageTransHelper(
                filePath + str(int(np.floor(angle))), center, width, col)
            ceil = _averageTransHelper(
                filePath + str(int(np.ceil(angle))), center, width, col)
            prop = angle - np.floor(angle)
            return floor * (1 - prop) + ceil * prop
        else:
            print("Angle out of range")
    else:
        if angle >= 30 and angle <= 65 and int(angle/5) == angle/5:
            return _averageTransHelper(filePath + str(angle), center, width, col)
        elif int(angle/5) != angle/5:
            floor = _averageTransHelper(
                filePath + str(int(np.floor(angle/5)*5)), center, width, col)
            ceil = _averageTransHelper(
                filePath + str(int(np.ceil(angle/5)*5)), center, width, col)
            prop = angle/5 - np.floor(angle/5)
            return floor * (1 - prop) + ceil * prop
        else:
            print("Angle out of range")


def _getEQTrans(angle, center, width):
    return np.array([[_averageTrans("CerroConfig/", angle, percentile, centeri, widthi) for percentile in [25, 50, 75]] for centeri, widthi in zip(center, width)])


def _getEQTransV2(angle, center, width):
    #print(center, width)
    lower = np.array([[_averageTrans(os.path.join(absolute_path, "data/VariablePWV/"), angle, pwvTick, centeri, widthi, newFilePathFormat=True)
                     for pwvTick in [_pwvToTick(ccatPWVQ1), _pwvToTick(ccatPWVQ2), _pwvToTick(ccatPWVQ3)]] for centeri, widthi in zip(center, width)])
    higher = np.array([[_averageTrans(os.path.join(absolute_path, "data/VariablePWV/"), angle, pwvTick, centeri, widthi, newFilePathFormat=True)
                      for pwvTick in [_pwvToTick(ccatPWVQ1)+1, _pwvToTick(ccatPWVQ3)+1]] for centeri, widthi in zip(center, width)])
    propQ1 = ccatPWVQ1/(ccatPWVQ2/20)-_pwvToTick(ccatPWVQ1)
    propQ3 = ccatPWVQ3/(ccatPWVQ2/20)-_pwvToTick(ccatPWVQ3)
    trueQ1 = lower[:, 0]*(1-propQ1)+higher[:, 0]*propQ1
    trueQ3 = lower[:, 2]*(1-propQ3)+higher[:, 1]*propQ3

    linearApprox = np.array([trueQ1, lower[:, 1], trueQ3]).T
    return linearApprox


def _trun(array, decimalPlaces):
    if decimalPlaces != None:
        for k1 in array:
            if type(array[k1]) == dict:
                for k2 in array[k1]:
                    array[k1][k2] = round(array[k1][k2], decimalPlaces)
            else:
                array[k1] = round(array[k1], decimalPlaces)
        return array
    else:
        return array


def _valueDisplayHelper(array, w, dict, outputFreq, decimalPlaces):
    if len(w) > 0:
        unit = None
        if outputFreq:
            unit = " GHz"
        else:
            unit = " um"
        dict.update({str(int(w[0])) + unit: array[0]})
        return _valueDisplayHelper(array[1:], w[1:], dict, outputFreq, decimalPlaces)
    else:
        return _trun(dict, decimalPlaces)


# Creates a dictionary to be put into a yaml file for broadband data
def _valueDisplay(array, outputFreq, centerFrequency, wavelength, decimalPlaces):
    if type(array) == np.ndarray:
        array = _arrayify(array)
    outputLabel = None
    if outputFreq:
        outputLabel = centerFrequency / 1e9
    else:
        outputLabel = wavelength
    return _valueDisplayHelper(array, outputLabel, {}, outputFreq, decimalPlaces)


def valDisplayPartial(outputFreq, centerFrequency, wavelength, decimalPlaces):
    """Returns a function that takes a 1d array as input and returns a dictionary marking the data with frequency/wavelength information. [outputFreq] is a boolean flag for frequency/wavelength output."""
    return partial(_valueDisplay, outputFreq=outputFreq, centerFrequency=centerFrequency, wavelength=wavelength, decimalPlaces=decimalPlaces)


def _quartileDisplayHelper(array, w, dict, outputFreq, decimalPlaces):
    if len(w) > 0:
        unit = None
        if outputFreq:
            unit = " GHz"
        else:
            unit = " um"
        dict.update({str(int(w[0])) + unit: {"Quartile 1": array[0][0], "Quartile 2": array[0]
                    [1], "Quartile 3": array[0][2]}})
        return _quartileDisplayHelper(array[1:], w[1:], dict, outputFreq, decimalPlaces)
    else:
        return _trun(dict, decimalPlaces)


# Creates a dictionary to be put into a yaml file for data involving quartiles
def _quartileDisplay(array, outputFreq, centerFrequency, wavelength, decimalPlaces):
    if type(array) == np.ndarray:
        array = _arrayify(array)
    outputLabel = None
    if outputFreq:
        outputLabel = centerFrequency / 1e9
    else:
        outputLabel = wavelength
    return _quartileDisplayHelper(array, outputLabel, {}, outputFreq, decimalPlaces)


def quartDisplayPartial(outputFreq, centerFrequency, wavelength, decimalPlaces):
    """Returns a function that takes a 2d array as input (representing frequency on one axis and quartiles on the other) and returns a dictionary marking frequency and quartile information. [outputFreq] is a boolean flag for frequency/wavelength output."""
    return partial(_quartileDisplay, outputFreq=outputFreq, centerFrequency=centerFrequency, wavelength=wavelength, decimalPlaces=decimalPlaces)


def calcByAngle(diameter, t, wfe, eta, doe, t_int, pixelYield, szCamNumPoln, eorSpecNumPoln, t_filter_cold, t_lens_cold, t_uhdpe_window, coldSpillOverEfficiency, singleModedAOmegaLambda2, spatialPixels, eqbw, centerFrequency, detectorNEP, backgroundSubtractionDegradationFactor, r):
    """Returns a function that takes observation zenith angle as an input and returns the following outputs in a dictionary: NET Weighted Average as "netW8Avg", NET Weighted RJ as "netW8RJ", NEI Weighted Average as "neiW8", EoR Spec NEFD as "eorNEFD", EoR Spec NEI as "eorNEI", Power per Pixel as "powerPerPixel", EoR Spec Power per Pixel as "eorPowerPerPixel", Center Wavelengths as "wavelength", Beam as "beam"."""

    partTrans = partial(_getEQTransV2, center=centerFrequency, width=eqbw)
    partCalc = partial(_calculate, diameter, t, wfe, eta, doe, t_int, pixelYield, szCamNumPoln, eorSpecNumPoln, t_filter_cold, t_lens_cold, t_uhdpe_window, coldSpillOverEfficiency,
                       singleModedAOmegaLambda2, spatialPixels, eqbw, centerFrequency, detectorNEP, backgroundSubtractionDegradationFactor, r)
    return lambda x: partCalc(partTrans(x))


def outputSensitivityFile(outputs, valueDisplay, quartileDisplay):
    """Outputs output.yaml"""
    dict_file = {"NET w8 avg": valueDisplay(outputs["netW8Avg"]),
                 "NET w8 RJ": valueDisplay(outputs["netW8RJ"]), "NEI w8 Jy/sr": valueDisplay(
        outputs["neiW8"]), "EoR Spec NEFD": quartileDisplay(outputs["eorNEFD"]), "EoR Spec NEI": quartileDisplay(outputs["eorNEI"])}
    return yaml.dump(dict_file, open("output.yaml", 'w'), sort_keys=False)


def outputPowerFile(i, outputs, calculate, quartileDisplay):
    """Outputs power.yaml"""
    outputs60 = calculate(30)
    outputs45 = calculate(45)

    def displayPower(o):
        return {"Broadband": quartileDisplay(o["powerPerPixel"]*10**12), "EoR Spec": quartileDisplay(o["eorPowerPerPixel"]*10**12)}

    def powerDiff(output1, output2):
        return {"powerPerPixel": np.array([[q1-q2 for q1, q2 in zip(f1, f2)] for f1, f2 in zip(output1["powerPerPixel"], output2["powerPerPixel"])]), "eorPowerPerPixel": np.array([[q1-q2 for q1, q2 in zip(f1, f2)] for f1, f2 in zip(output1["eorPowerPerPixel"], output2["eorPowerPerPixel"])])}

    dict_file = {str(i["observationElevationAngle"]) + " degrees": displayPower(outputs),
                 "45 degrees": displayPower(outputs45), "60 degrees": displayPower(outputs60), "45-60 Delta": displayPower(powerDiff(outputs45, outputs60))}
    return yaml.dump(dict_file, open("power.yaml", 'w'), sort_keys=False)


def _calcSpillFromData(half_angle, contractFactor, degrees, values, showPlots=False, f=None):
    def power2(db):
        return 10.0**(db/10.0)

    def invPower(x):
        return np.log10(x)*10

    myDegrees = degrees * contractFactor
    th = np.radians(myDegrees)

    # Normalize the data
    data = invPower(power2(values)/power2(np.max(values)))

    tot_cutoff = np.where(np.abs(th) <= np.radians(180))[0]
    tot_cutoff = tot_cutoff[tot_cutoff < len(data)]
    tot = np.trapz(power2(data[tot_cutoff]) *
                   np.sin(th[tot_cutoff]), th[tot_cutoff])

    th_cutoff = np.where(np.abs(th) < np.radians(half_angle))[0]

    beam = np.trapz(power2(data[th_cutoff]) *
                    np.sin(th[th_cutoff]), th[th_cutoff])

    spill_eff = beam/tot

    if showPlots:
        la = ""
        if f != None:
            la = " at %d GHz" % f
        plt.plot(np.degrees(th)[:721], power2(
            values)/np.max(power2(values)), label="Doug phi = 0" + la, linewidth=2)
        plt.axvline(x=half_angle, color='k',
                    linewidth=2, label='Lyot stop angle')
        plt.legend(loc=0)
        plt.xlim(0, 180)
        plt.yscale('log')
        plt.xlabel('angle [deg]')
        plt.ylabel('normalized beam')
        plt.show()
        plt.clf()
    return spill_eff


def _getColdSpillOverEfficiency(i, beamFreq, beamPixelSpacing, degrees, values, showPlots=False):
    return np.array([_calcSpillFromData(i["lyotStopAngle"], beamFreq / (f * (s / beamPixelSpacing)), degrees, values, showPlots=showPlots, f=f/1e9) for f, s in zip(i["centerFrequency"], i["detectorSpacing"])])


def outputSpillEfficiencyFile(i, calculate, spillEfficiency):
    """Prints a chart comparing calculated spill efficiencies to the excel values of spill efficiencies to stdout."""
    output = calculate(45)
    t = Texttable(max_width=110)
    t.set_cols_dtype(['i', (lambda x: "%.1f" % float(x)), 'f', 'i',
                     (lambda x: "%.2f" % float(x)), (lambda x: "%.1f" % float(x)), 'f'])
    t.set_cols_align(['_c', '_c', '_c', '_c', '_c', '_c', '_c'])
    excelSpillEff = [.8, .5, .7, .5, .5][::-1]
    excelNET = [241440.5, 181.8, 56.3, 11.4, 6.8][::-1]
    t.add_rows(np.concatenate((np.reshape(['Center Frequency (GHz)', 'Excel Spill Efficiency', 'Calculated Spill Efficiency', '# Pixels', 'Pixel Spacing (mm)', 'Excel NET (uK rt(s))', 'Calculated NET (uK rt(s))'], (-1, 1)),
                               np.array([[int(f/1e9), exSpillEf, spillEf, numDetect, detectSpacing, exNet, net] for f, spillEf, net, exSpillEf, numDetect, detectSpacing, exNet in zip(i['centerFrequency'], spillEfficiency, output['netW8Avg'], excelSpillEff, i["spatialPixels"], i['detectorSpacing'], excelNET)])[::-1].T), axis=1).T)
    print(t.draw())


# Load the beam file and pass the angle and value data to other functions
def getSpillEfficiency(i, oldFile=False):
    """Returns spill efficiency, using data/tolTEC_staircase_singleHorn_280GHz.txt as a reference. Is meant to be updated when better/more curves are calculated."""
    if oldFile:
        data = np.genfromtxt(os.path.join(absolute_path,
                                          'data/tolTEC_staircase_singleHorn_280GHz.txt'), skip_header=2).reshape(-1, 721, 8)
        degr = data[0, :, 0]
        vals = data[0, :, 3]
        # data = np.genfromtxt('data/beam_280.txt')
        # degr = np.degrees(data[:, 0])
        # vals = np.log10((data[:, 1]**2))*10
        return _getColdSpillOverEfficiency(
            i, 280e9, 2.75, degr, vals, showPlots=False)
    else:
        data = np.genfromtxt(
            os.path.join(absolute_path, 'data/ccat350_2p75_pitch_250um_step_v1run10_12AUG2022_beam_350GHz.txt'))
        degr = data[:, 0]
        vals = np.log10((data[:, 1]**2))*10
        return _getColdSpillOverEfficiency(i, 350e9, 2.75, degr, vals, showPlots=False)


def _averageTemp(start, end, col, filePath=None, prefix=None, angle=None, percentile=None):
    if filePath is None:
        filePath = f"{prefix}/{percentile}/ACT_annual_{percentile}.{angle}"
    file = open(filePath + ".out", "r")
    data = file.readlines()
    x = [float(i.split(" ")[0]) for i in data]
    y = [float(i.split(" ")[col]) for i in data]
    file.close()
    d_nu = x[1]-x[0]
    ghzStart = start/1e9
    ghzEnd = end/1e9

    temperature = 0
    for (f, temp) in zip(x, y):
        if f < ghzStart:
            continue
        if f >= ghzEnd:
            break
        temperature += temp*(f*10**9)**2*d_nu

    corr = 0
    for f in x:
        if f < ghzStart:
            continue
        if f >= ghzEnd:
            break
        corr += (f*10**9)**2*d_nu

    return temperature/corr


def _tangent_line_slope(cf, eqbw, col, a, p, maunaKea=False, filePath=None):
    """Depricated. Requires higher and lower calculations, should just use variable pwv."""
    higher = ""
    lower = ""
    if not maunaKea:
        higher = "Higher"
        lower = "Lower"
    else:
        higher = "MaunaKea/Higher"
        lower = "MaunaKea/Lower"
    # V1 Derivatives
    higher = np.array([_averageTemp(c-w/2, c+w/2, col, prefix=higher, angle=a, percentile=p, filePath=filePath)
                       for c, w in zip(cf, eqbw)])
    lower = np.array([_averageTemp(c-w/2, c+w/2, col, prefix=lower, angle=a, percentile=p, filePath=filePath)
                      for c, w in zip(cf, eqbw)])
    derivative = (higher - lower) / 0.02
    return derivative


def _tangent_line_slopeV2(cf, eqbw, col, a, pwv):
    lower = int(pwv * 20 / ccatPWVQ2)
    higher = lower + 1
    lower = np.array([_averageTemp(c-w/2, c+w/2, col, filePath=("VariablePWV/ACT_annual_" + str(lower) + "." + str(a)))
                      for c, w in zip(cf, eqbw)])
    higher = np.array([_averageTemp(c-w/2, c+w/2, col, filePath=("VariablePWV/ACT_annual_" + str(higher) + "." + str(a)))
                       for c, w in zip(cf, eqbw)])

    derivative = (higher - lower) / (ccatPWVQ2 / 20)

    return derivative


def _least_squares_slope(cf, eqbw, col, a, graph=False, maunaKea=False):
    # V2 Derivatives
    filePath = "VariablePWV/ACT_annual_"
    if maunaKea:
        filePath = "MaunaKea/VariablePWV/"
    wieghtedTemp = np.array([[_averageTemp(c-w/2, c+w/2, col, filePath=(filePath + str(i) + "." + str(a)))
                              for c, w in zip(cf, eqbw)] for i in np.array(range(40))+1])
    rjTemp = np.array([[_averageTransSE("VariablePWV/ACT_annual_" + str(i) + ".45", (c-w/2) /
                      1e9, (c+w/2)/1e9, col=2) for c, w in zip(cf, eqbw)] for i in np.array(range(40))+1])
    planckTemp = np.array([[_averageTransSE("VariablePWV/ACT_annual_" + str(i) + ".45", (c-w/2) /
                                            1e9, (c+w/2)/1e9, col=3) for c, w in zip(cf, eqbw)] for i in np.array(range(40))+1])
    absorption = 1-np.array([[_averageTransSE("VariablePWV/ACT_annual_" + str(i) + ".45",
                                              (c-w/2)/1e9, (c+w/2)/1e9, col=1) for c, w in zip(cf, eqbw)] for i in np.array(range(40))+1])

    def line(x, m, b):
        return m*x + b
    derivative = []
    temp = ccatPWVQ2
    if maunaKea:
        temp = ACTPWV.configPWVHelper(
            "data/MaunaKea/Default/50/.err")
        # print(ccatMedPWV)
    pwvs = (np.array(range(40))+1) / 20*temp
    for i in range(len(cf)):
        popt = op.curve_fit(line, pwvs, wieghtedTemp[:, i])[0]
        # print(popt)
        derivative.append(popt[0])
    derivative = np.array(derivative)
    if graph:
        ylabel = ["Weighted Rayleigh-Jeans Temperature (K)", "RJ Temperature (K)",
                  "Planck Temperature (K)", "Absorption"]
        if not maunaKea:
            print("Steve's PWV range:", (0.3*.7192506910, 3*.7192506910))
            print("CCAT 50th percentile PWV:", ccatPWVQ2)
        for j in range(4):
            for i in np.array(range(len(cf)))[::-1][:-1]:
                if j == 0:
                    plt.plot(pwvs, wieghtedTemp[:, i], linewidth=1,
                             label=str(int(cf[i]/1e9))+' GHz')
                if j == 1:
                    plt.plot(pwvs, rjTemp[:, i], linewidth=1,
                             label=str(int(cf[i]/1e9))+' GHz')
                if j == 2:
                    plt.plot(pwvs, planckTemp[:, i], linewidth=1,
                             label=str(int(cf[i]/1e9))+' GHz')
                if j == 3:
                    plt.plot(pwvs, absorption[:, i], linewidth=1,
                             label=str(int(cf[i]/1e9))+' GHz')
            plt.title("Brightness Temperature vs PWV")
            plt.ylim(bottom=0)
            plt.ylabel(ylabel[j])
            plt.xlim(left=0, right=ccatPWVQ2*2)
            plt.xlabel("PWV (mm)")
            plt.legend(loc='best')
            plt.grid()
            plt.show()

        for i in [1, 5]:  # 220 and 850 GHz
            plt.plot(pwvs, wieghtedTemp[:, i], linewidth=1,
                     label=str(int(cf[i]/1e9))+' GHz Weighted RJ Temp')
            plt.plot(pwvs, rjTemp[:, i], linewidth=1,
                     label=str(int(cf[i]/1e9))+' GHz RJ Temp')
            plt.plot(pwvs, planckTemp[:, i], linewidth=1,
                     label=str(int(cf[i]/1e9))+' GHz Planck Temp')
        plt.ylim(bottom=0)
        plt.ylabel("Brightness Temperature (K)")
        plt.xlim(left=0, right=ccatPWVQ2*2)
        plt.xlabel("PWV (mm)")
        plt.legend(loc='best')
        plt.title("Brightness Temperature vs PWV Comparison")
        plt.grid()
        plt.show()

        for i in np.array(range(len(cf)))[::-1]:
            plt.plot(pwvs, wieghtedTemp[:, i] - rjTemp[:, i], linewidth=1,
                     label=str(int(cf[i]/1e9))+' GHz')
        plt.ylabel("Difference Between Weighted RJ and RJ (Arbitrary Units)")
        plt.xlim(left=0, right=ccatPWVQ2*2)
        plt.xlabel("PWV (mm)")
        plt.legend(loc='best')
        plt.grid()
        plt.show()

    return derivative


def _least_squares_slopeV2(cf, eqbw, col, a, lowPWV, highPWV):
    wieghtedTemp = np.array([[_averageTemp(c-w/2, c+w/2, col, filePath=(os.path.join(absolute_path, "data/VariablePWV/ACT_annual_") + str(i) + "." + str(a)))
                              for c, w in zip(cf, eqbw)] for i in np.array(range(40))+1])

    def line(x, m, b):
        return m*x + b
    derivative = []
    pwvs = (np.array(range(40))+1) / 20*ccatPWVQ2

    valid = np.all([np.where(pwvs > lowPWV, True, False),
                   np.where(pwvs < highPWV, True, False)], axis=0)
    for i in range(len(cf)):
        popt = op.curve_fit(line, pwvs[valid], wieghtedTemp[:, i][valid])[0]
        derivative.append(popt[0])
    derivative = np.array(derivative)

    return derivative


def _data_C_calc(i, table=False, graphSlopes=False, maunaKea=False):
    P = 50
    A = 45
    COL = 2

    cf = np.append(145e9, i["centerFrequency"])
    eqbw = np.append(145*0.276e9, i["eqbw"])

    dataCs = np.array([])
    # Choose method
    derivative = _tangent_line_slope(cf, eqbw, COL, A, P, maunaKea=maunaKea)
    for n in range(1, len(cf)):
        dataCs = np.append(dataCs,
                           (derivative[n]*_a_to_CMB(cf[n]/1e9)/(derivative[0]*_a_to_CMB(cf[0]/1e9)))**2)
    if table:
        t = Texttable(max_width=110)
        table_header = np.append("Method", np.char.add(
            (cf/1e9).astype(int).astype(str), ' GHz'))
        table = np.array([table_header, np.append("Least Squares Regression Line", _least_squares_slope(cf, eqbw, COL, A, graph=graphSlopes, maunaKea=maunaKea)), np.append("Best Case Tangent Line", _tangent_line_slope(cf, eqbw, COL, A, 75, maunaKea=maunaKea)), np.append("Expected Tangent Line", _tangent_line_slope(cf, eqbw, COL, A, 50, maunaKea=maunaKea)), np.append("Worst Case Tangent Line", _tangent_line_slope(cf, eqbw, COL, A, 25, maunaKea=maunaKea)), np.append("Corrected PWV Previous Method",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     np.array([6.2, 14.7, 25.0, 48.5, 66.4, 64.4])), np.append("Uncorrected PWV Previous Method",
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               np.array([4.5, 10.6, 17.9, 34.9, 47.8, 46.3]))])
        t.add_rows(table, header=True)
        print(t.draw())
    elif graphSlopes:
        _least_squares_slope(cf, eqbw, COL, A, graph=True, maunaKea=maunaKea)
    return dataCs*1.2e4


def _data_C_calcV2(i):
    A = 45
    COL = 2

    cf = np.append(145e9, i["centerFrequency"])
    eqbw = np.append(145*0.276e9, i["eqbw"])

    dataCs = np.array([])
    derivative = _least_squares_slopeV2(cf, eqbw, COL, A, ccatPWVQ1, ccatPWVQ3)
    for n in range(1, len(cf)):
        dataCs = np.append(dataCs,
                           (derivative[n]*_a_to_CMB(cf[n]/1e9)/(derivative[0]*_a_to_CMB(cf[0]/1e9)))**2)
    return dataCs*1.2e4


def outputNoiseCurvesAndTempVsPWV(i, outputs, calculate='all', plotCurve=None, table=False, graphSlopes=False, maunaKea=False, lowFreq=False):
    """NOTE: Deprecated from removal of old am data."""
    centerFrequency = None
    beam = None
    net = None
    data_C = None
    if calculate == 'original':
        centerFrequency = [222., 280., 348., 405., 850.]
        beam = [59/60., 47/60., 37/60., 32/60., 15/60.]
        net = [6.8, 12.7, 47.7, 181.8, 305400.7]
        data_C = None  # Automatic values
    elif calculate == 'change' or calculate == 'dataC' or calculate == 'data_C':
        centerFrequency = i['centerFrequency']/1e9
        beam = outputs["beam"]/60
        net = outputs["netW8Avg"]
        data_C = None
    elif calculate == 'y' or calculate == 'yes' or calculate == 'all':
        centerFrequency = i['centerFrequency']/1e9
        beam = outputs["beam"]/60
        net = outputs["netW8Avg"]
        data_C = _data_C_calc(
            i, table=table, graphSlopes=graphSlopes, maunaKea=maunaKea)
        # print(data_C)
    else:
        print("Select a valid calculate option")
        exit(1)
    ccat = noise_file.CCAT(centerFrequency, beam, net, survey_years=4000 /
                      24./365.24, survey_efficiency=1.0, N_tubes=(1, 1, 1, 1, 1), el=45., data_C=data_C)
    fsky = 20000./(4*_pi*(180/_pi)**2)
    lat_lmax = 10000
    ell, N_ell_T_full, N_ell_P_full = ccat.get_noise_curves(
        fsky, lat_lmax, 1, full_covar=False, deconv_beam=True)

    # Formerly under the if statement
    plotTemperature = plotCurve == 'T' or plotCurve == 't' or plotCurve == 'temp' or plotCurve == 'temperature'
    if False:  # Removes 850 GHz
        N_ell_T_full = N_ell_T_full[:-1]
        N_ell_P_full = N_ell_P_full[:-1]
        centerFrequency = centerFrequency[:-1]
        plt.ylim(10**-5, 10**3)
    elif lowFreq:  # Only 280 GHz
        N_ell_T_full = [N_ell_T_full[1]]
        N_ell_P_full = [N_ell_P_full[1]]
        centerFrequency = [centerFrequency[1]]
        plt.ylim(10**-4, 10**3)
    else:  # Only 850 GHz
        N_ell_T_full = [N_ell_T_full[4]]
        N_ell_P_full = [N_ell_P_full[4]]
        centerFrequency = [centerFrequency[4]]
        plt.ylim(10**4, 10**12)
    if plotCurve is not None:
        for curve, label in zip(N_ell_T_full if plotTemperature else N_ell_P_full, centerFrequency):
            plt.plot(ell, curve, label=str(int(label))+' GHz')
            plt.yscale('log')
            plt.xscale('log')
            plt.xlim(10**2, 10**4)
            plt.title("Temperature" if plotTemperature else "Polarization")
        plt.legend(loc='upper right')
        plt.grid()
        plt.show()
    return ell, N_ell_T_full[0], N_ell_P_full[0]


def getCustNoiseCurvesSubplot(i, outputs, temp, lowFreq, ax):
    before = outputNoiseCurvesAndTempVsPWV(i, outputs, calculate='change', plotCurve=None,
                                           table=False, graphSlopes=False, maunaKea=False, lowFreq=lowFreq)
    after = outputNoiseCurvesAndTempVsPWV(i, outputs, calculate='all', plotCurve=None,
                                          table=False, graphSlopes=False, maunaKea=False, lowFreq=lowFreq)

    ell = before[0]  # Also equal to after[0]

    ax.plot(ell, before[1 if temp else 2], label='Before')
    ax.plot(ell, after[1 if temp else 2],
            ('-' if temp else '--'), label='After')

    ax.set_yscale('log')
    if lowFreq:
        ax.set_ylim(10**-4, 10**3)
    else:
        ax.set_ylim(10**4, 10**12)
    ax.set_ylabel('$N_{\ell} (\mu K^2)$')
    ax.set_xscale('log')
    ax.set_xlim(10**2, 10**4)
    ax.set_xlabel('$\ell$')
    ax.set_title(('Temperature' if temp else 'Polarization') +
                 ' Noise at ' + ('280' if lowFreq else '850') + ' GHz')
    ax.legend(loc='upper right')
    ax.grid()


def outputNoiseCurves(i, outputs):
    """NOTE: Deprecated from removal of old am data."""
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
        2, 2, sharex='col', sharey='row')
    getCustNoiseCurvesSubplot(i, outputs, True, True, ax1)
    getCustNoiseCurvesSubplot(i, outputs, False, True, ax2)
    getCustNoiseCurvesSubplot(i, outputs, True, False, ax3)
    getCustNoiseCurvesSubplot(i, outputs, False, False, ax4)
    for ax in fig.get_axes():
        ax.label_outer()
    fig.tight_layout()
    plt.show()


def outputLoadings(i, calculate):
    """Output the loadings at 45 and 60 degree angles of elevation for Q1, Q2, and Q3 PWV at all frequencies."""
    out45 = calculate(45)
    out60 = calculate(30)

    t = Texttable(max_width=110)
    t.set_cols_dtype(['t', 'e', 'e', 'e', 'e', 'e'])
    table_header = np.append("Center Frequency", np.char.add(
        (i["centerFrequency"]/1e9).astype(int).astype(str), ' GHz'))
    table = np.array([table_header, np.append("Bandwidths", i["eqbw"]), np.append("Q1 45 Degree Loading", out45["powerPerPixel"][:, 0]), np.append("Q2 45 Degree Loading", out45["powerPerPixel"][:, 1]), np.append(
        "Q3 45 Degree Loading", out45["powerPerPixel"][:, 2]), np.append("Q1 60 Degree Loading", out60["powerPerPixel"][:, 0]), np.append("Q2 60 Degree Loading", out60["powerPerPixel"][:, 1]), np.append("Q3 60 Degree Loading", out60["powerPerPixel"][:, 2])])
    t.add_rows(table, header=True)
    print(t.draw())


def useLatexFont():
    rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
    rc('text', usetex=True)


def getNoiseCurves(i, outputs):
    """Returns a tuple of (ell, N_ell_T_full, N_ell_P_full) for the given inputs and outputs from the rest of the sensitivity calculator."""
    centerFrequency = i['centerFrequency']/1e9
    beam = outputs["beam"]/60
    net = outputs["netW8Avg"]
    data_C = _data_C_calcV2(i)
    ccat = noise_file.CCAT(centerFrequency, beam, net, survey_years=4000 /
                      24./365.24, survey_efficiency=1.0, N_tubes=tuple(1 for _ in centerFrequency), el=45., data_C=data_C)
    fsky = 20000./(4*_pi*(180/_pi)**2)
    lat_lmax = 10000
    ell, N_ell_T_full, N_ell_P_full = ccat.get_noise_curves(
        fsky, lat_lmax, 1, full_covar=False, deconv_beam=True)
    return ell, N_ell_T_full, N_ell_P_full


def outputDataCChanges(i):
    """NOTE: Deprecated from removal of old am data."""
    old = _data_C_calc(i)
    new = _data_C_calcV2(i)
    change = 100*(old-new)/new

    t = Texttable(max_width=110)
    t.set_cols_dtype(['t', 'e', 'e', 'e', 'e', 'e'])
    t.set_cols_width([20, 10, 10, 10, 10, 10])
    table_header = np.append("Center Frequency", np.char.add(
        (i["centerFrequency"]/1e9).astype(int).astype(str), ' GHz'))
    table = np.array([table_header, np.append("5% to 200% of Q2 PWV", old), np.append(
        "Q1 PWV to Q3 PWV", new)])
    t.add_rows(table, header=True)
    print(t.draw())

    t = Texttable(max_width=110)
    t.set_cols_dtype(['t', 'f', 'f', 'f', 'f', 'f'])
    t.set_precision(1)
    t.set_cols_width([20, 10, 10, 10, 10, 10])
    table = np.array([np.append("Percentage Change", change)])
    t.add_rows(table, header=False)
    print(t.draw())


def _eorCalculate(diameter, t, wfe, eta, doe, eorSpecNumPoln, t_filter_cold, t_lens_cold, t_uhdpe_window, coldSpillOverEfficiency, centerFrequency, detectorNEP, backgroundSubtractionDegradationFactor, eqtrans, spatialPixels, pixelYield, eorEqBw, detIndex, fnum):
    temp = eorEqBw[fnum]
    eorEqBw = np.zeros(len(spatialPixels))
    eorEqBw[detIndex] = temp
    eorEqBw = eorEqBw[:, None]

    wavelength = _c/centerFrequency*10**6

    # Telescope
    a = _pi*(diameter/2)**2

    # Detector
    eorSpecCorForPoln = abs(3-eorSpecNumPoln)

    # Instrument, throughput
    t_cold = t_filter_cold*t_lens_cold**3
    e_window_warm = 1 - t_uhdpe_window

    # Instrument, beams/area of focal plane
    beam = 1.2*wavelength/diameter/1000000*206265
    solidAngle = _pi/4/_ln(2)*(beam/206264)**2

    # EoR Spec
    eorEqTrans = eqtrans
    eorE_warm = (t_uhdpe_window)[:, None]*((1-eorEqTrans)
                                           * eta+(1-eta))+(e_window_warm)[:, None]
    eorRuze = 1/_e**((4*_pi*wfe/(wavelength))**2)
    # Always takes from 739 um in the sheet, also weird formula
    eorT_cold = (t_cold)*0.9
    eorOccupancy = 1/(_e**(_h*_c/((wavelength)*10**(-6)*_k*t))-1)
    eorAcceptedModes = eorSpecNumPoln*a*(coldSpillOverEfficiency)*(solidAngle)/(
        (wavelength)*10**(-6))**2  # In sheet 1071 um is calculated incorrectly
    eorPhotonNoiseNEP = _h*_c/((wavelength)[:, None]*10**(-6))*(eorAcceptedModes[:, None]*eorEqBw*eorE_warm *
                                                                eorT_cold[:, None]*doe*eorOccupancy[:, None]*(1+eorE_warm*eorT_cold[:, None]*doe*eorOccupancy[:, None]))**0.5
    eorTotalNEP = _sqrt(eorPhotonNoiseNEP**2 + detectorNEP**2)
    eorColdTerminatedSpillover = (coldSpillOverEfficiency)
    eorNEF = eorSpecCorForPoln*eorTotalNEP/(a*eorT_cold[:, None]*doe*eta*eorRuze[:, None]
                                            * eorEqTrans*eorColdTerminatedSpillover[:, None]*backgroundSubtractionDegradationFactor)
    eorNEFD = eorNEF/eorEqBw*10**26*1000
    eorNEI = eorNEFD/(solidAngle)[:, None]/1000
    eorPowerPerPixel = _h*_c/(wavelength[:, None]*10**(-6))*eorE_warm*eorT_cold[:, None]*eorAcceptedModes[:,
                                                                                                          None]*eorOccupancy[:, None]*eorEqBw*doe  # Look at t_cold in excel sheet

    # Inference from non-eor stuffs
    net = eorNEFD*10**(-29)/((solidAngle)[:, None])*((wavelength)
                                                     [:, None]*0.000001)**2/(2*1.38*10**(-23))*1000
    arrayNETRJ = (net/((spatialPixels)[:, None]*pixelYield)**0.5)
    aToCMB = np.array([_a_to_CMB(i/1e9) for i in centerFrequency])
    arrayNETCMB = arrayNETRJ*aToCMB[:, None]*1000
    netw8avg = _sqrt(3/(1/arrayNETCMB[:, 0]**2+1/arrayNETCMB[:, 1]
                        ** 2+1/arrayNETCMB[:, 2]**2).astype(float))
    return {"eorNEFD": eorNEFD, "eorNEI": eorNEI, "eorPowerPerPixel": eorPowerPerPixel, "wavelength": wavelength, "beam": beam, "netW8Avg": netw8avg}


def _eorCalcByFreqR(angle, diameter, t, wfe, eta, doe, eorSpecNumPoln, t_filter_cold, t_lens_cold, t_uhdpe_window, coldSpillOverEfficiency, detectorNEP, backgroundSubtractionDegradationFactor, spatialPixels, pixelYield, eqbw):
    """Returns a function that takes observation zenith angle as an input and returns the following outputs in a dictionary: EoR Spec NEFD as "eorNEFD", EoR Spec NEI as "eorNEI", Power per Pixel as "powerPerPixel", EoR Spec Power per Pixel as "eorPowerPerPixel", Center Wavelengths as "wavelength", Beam as "beam"."""
    partTrans = partial(_getEQTransV2, angle=angle)
    partCalc = partial(_eorCalculate, diameter=diameter, t=t, wfe=wfe, eta=eta, doe=doe, eorSpecNumPoln=eorSpecNumPoln,  t_filter_cold=t_filter_cold,  t_lens_cold=t_lens_cold,  t_uhdpe_window=t_uhdpe_window,
                       coldSpillOverEfficiency=coldSpillOverEfficiency, detectorNEP=detectorNEP,  backgroundSubtractionDegradationFactor=backgroundSubtractionDegradationFactor, spatialPixels=spatialPixels, pixelYield=pixelYield, eorEqBw=eqbw)
    return lambda f, detIndex, fnum: partCalc(centerFrequency=f, detIndex=detIndex, fnum=fnum, eqtrans=partTrans(center=f, width=eqbw))


def _geteqbw(f, r):
    return f / (10 * r) * _pi / 2


def eorNoiseCurves(i, rfpairs, frequencyRanges=np.array([[210, 315], [315, 420]])*10**9):
    """Returns a list of noise curves (ell, N_ell_T_full, N_ell_P_full) corresponding to rfpairs, a list of finesse-frequency pairs, and detector arrays i. frequencyRanges is an optional parameter 
    for modifying the frequency range of the detector arrays."""
    angle = 90 - i["observationElevationAngle"]
    coldSpillOverEfficiency = getSpillEfficiency(i)

    i["eqbw"] = _geteqbw(rfpairs[:, 1], rfpairs[:, 0])
    calculate = _eorCalcByFreqR(angle, i["diameter"], i["t"], i["wfe"], i["eta"], i["doe"], i["eorSpecNumPoln"],
                                i["t_filter_cold"], i["t_lens_cold"], i["t_uhdpe_window"], coldSpillOverEfficiency,
                                i["detectorNEP"], i["backgroundSubtractionDegradationFactor"], i["spatialPixels"], i["pixelYield"], i["eqbw"])

    def _outsourcedCalculation(i, output):
        centerFrequency = i["centerFrequency"] / 1e9
        beam = output["beam"] / 60
        netw8avg = output["netW8Avg"]
        data_C = _data_C_calcV2(i)
        ccat = noise_file.CCAT(centerFrequency, beam, netw8avg, survey_years=4000 /
                          24./365.24, survey_efficiency=1.0, N_tubes=tuple(1 for _ in centerFrequency), el=45., data_C=data_C)
        fsky = 20000./(4*_pi*(180/_pi)**2)
        lat_lmax = 10000
        ell, N_ell_T_full, N_ell_P_full = ccat.get_noise_curves(
            fsky, lat_lmax, 1, full_covar=False, deconv_beam=True)
        return ell, N_ell_T_full, N_ell_P_full

    results = {}
    for count, (ri, fi) in enumerate(rfpairs):
        index = 0
        for j in range(len(frequencyRanges)):
            if frequencyRanges[j][0] <= fi <= frequencyRanges[j][1]:
                index = j
                break
        else:
            print(f"Skipped {ri}, {fi} for not being in any frequency ranges")
            continue
        centerFrequency = np.zeros(len(i["spatialPixels"]))
        centerFrequency[index] = fi
        i["centerFrequency"] = centerFrequency

        r = np.zeros(len(i["spatialPixels"]))
        r[index] = ri
        i["r"] = r

        output = calculate(i["centerFrequency"], index, count)
        ell, N_ell_T_full, N_ell_P_full = _outsourcedCalculation(i, output)
        N_ell_T_full = N_ell_T_full[index]
        N_ell_P_full = N_ell_P_full[index]
        results[(ri, fi)] = (ell, N_ell_T_full, N_ell_P_full)
    return results


def spillEfficiencyComparison(lyotStopAngle=13.4, f=350e9, ds=2.75, maxangle=180):
    def spillPlot(half_angle, contractFactor, degrees, values, label):
        def power2(db):
            return 10.0**(db/10.0)

        def invPower(x):
            return np.log10(x)*10

        myDegrees = degrees * contractFactor
        th = np.radians(myDegrees)

        # Normalize the data
        data = invPower(power2(values)/power2(np.max(values)))

        tot_cutoff = np.where(np.abs(th) <= np.radians(180))[0]
        tot_cutoff = tot_cutoff[tot_cutoff < len(data)]
        tot = np.trapz(power2(data[tot_cutoff]) *
                       np.sin(th[tot_cutoff]), th[tot_cutoff])

        th_cutoff = np.where(np.abs(th) < np.radians(half_angle))[0]

        beam = np.trapz(power2(data[th_cutoff]) *
                        np.sin(th[th_cutoff]), th[th_cutoff])

        spill_eff = beam/tot
        plt.plot(np.degrees(th)[:721], power2(
            values)/np.max(power2(values)), label=label + f", spill efficiency = {spill_eff:.2f}", linewidth=2)
        plt.legend(loc=0)
        plt.xlim(0, maxangle)
        plt.yscale('log')
        plt.xlabel('angle [deg]')
        plt.ylabel('normalized beam')

    data = np.genfromtxt(
        os.path.join(absolute_path, 'data/tolTEC_staircase_singleHorn_280GHz.txt'), skip_header=2).reshape(-1, 721, 8)
    degr = data[0, :, 0]
    vals = data[0, :, 3]
    spillPlot(lyotStopAngle, (280e9 / f)
              * (ds / 2.75), degr, vals, "280 GHz beam")
    data = np.genfromtxt(
        os.path.join(absolute_path, 'data/ccat350_2p75_pitch_250um_step_v1run10_12AUG2022_beam_350GHz.txt'))
    degr = data[:, 0]
    vals = np.log10((data[:, 1]**2))*10
    spillPlot(lyotStopAngle, (350e9 / f)
              * (ds / 2.75), degr, vals, "350 GHz beam")
    plt.axvline(x=lyotStopAngle, color='k',
                linewidth=2, label='Lyot stop angle')
    plt.title(f"Beams scaled to {f/1e9:.0f} GHz and {ds} mm detector spacing")
    plt.show()
    plt.clf()

def ccat_mapsims(i, outputs, noiseCurves, channels, pysm_components, seed, sim_cmb=False, sim_noise=False): 
    instrument_path = Path("/home/amm487/cloned_repos/Sensitivity-Calculator/src/sensitivity_calculator/data/instrument_parameters/instrument_parameters.tbl")
    hitmap_path = "/home/amm487/cloned_repos/Sensitivity-Calculator/src/sensitivity_calculator/data/ccat_uniform_coverage_nside256_201021.fits"
    NSIDE = 256
    cmb = mapsims.SOPrecomputedCMB(
        num=seed,
        nside=NSIDE,
        lensed=False,
        aberrated=False,
        has_polarization=True,
        cmb_set=0,
        cmb_dir="data/mapsimscmb",
        input_units="uK_CMB",
    )
    ccat_survey = noise_file.CCAT(i["centerFrequency"], outputs["beam"], outputs["netW8Avg"])
    channels_list = mapsims.parse_channels(instrument_parameters=instrument_path)
    noise = mapsims.noise.ExternalNoiseSimulator(
        nside=NSIDE,
        return_uK_CMB=True,
        sensitivity_mode="baseline",
        apply_beam_correction=True,
        apply_kludge_correction=True,
        survey=ccat_survey,
        channels_list=channels_list
    )
    def smooth_map(m, n_it = 5, width=0.1):
        """Helper to get_window"""
        i = 0
        m_apo = m
        while i<= n_it:
            m_apo[m_apo<0.8] =0
            m_apo = hp.smoothing(m_apo, fwhm = width)
            m_apo[m_apo<0] = 0
            m_apo /= np.max(m_apo)
            i+=1
        return m_apo
    def get_window(mapk,n_it=5,width = 0.1):
        """Helper to apodize_map"""
        m_apo = np.copy(mapk)*0
        if hp.UNSEEN in np.copy(mapk):
            m_apo[np.copy(mapk)!=hp.UNSEEN] = 1.0
        else:
            m_apo[np.copy(mapk)!=0] = 1.0

        return smooth_map(m_apo, n_it=5, width = 0.1)
    def apodize_map(map0,n_it =5):
        """Apodizes a map with hp.UNSEEN pixels as masks. n_it represents the iterations of gaussian convolutions. Higher n_it generally means bigger windows and lower n_it the inverse."""

        tmp = np.copy(map0)
        tmp2 = np.copy(map0)
        tmp1 = tmp!=hp.UNSEEN
        m_apo = get_window(tmp2, n_it = n_it)
        tmp[tmp1!=True] =0
        output = tmp*m_apo
        return output

    chs = channels
    final = []

    for ch in chs:
        simulator = mapsims.MapSim(
            channels=ch,
            nside=NSIDE,
            unit="uK_CMB",
            pysm_output_reference_frame="C",
            pysm_components_string=pysm_components,
            pysm_custom_components={"cmb": cmb} if sim_cmb else None,
            other_components={"noise": noise} if sim_noise else None,
            instrument_parameters=instrument_path,
            num=seed
        )
        output_map = simulator.execute()
        for det in output_map.keys():
            for pol in np.arange(output_map[det].shape[0]):
                output_map[det][pol] = apodize_map(output_map[det][pol])
        final.append(output_map)
    pols = ["T", "Q", "U"]
    for h in final:
        for k in h.keys():
            print(k)
            for pol in np.arange(h[k].shape[0]):
                hp.mollview(h[k][pol], title = str(k) + " " + pols[pol])
                plt.show()

def so_mapsims(channels, pysm_components, seed, sim_cmb=False, sim_noise=False): 
    NSIDE = 256
    cmb = mapsims.SOPrecomputedCMB(
        num=seed,
        nside=NSIDE,
        lensed=False,
        aberrated=False,
        has_polarization=True,
        cmb_set=0,
        cmb_dir="data/mapsimscmb",
        input_units="uK_CMB",
    )
    noise = mapsims.SONoiseSimulator(
        nside=NSIDE,
        return_uK_CMB=True,
        sensitivity_mode="baseline",
        apply_beam_correction=True,
        apply_kludge_correction=True,
    )
    def smooth_map(m, n_it = 5, width=0.1):
        """Helper to get_window"""
        i = 0
        m_apo = m
        while i<= n_it:
            m_apo[m_apo<0.8] =0
            m_apo = hp.smoothing(m_apo, fwhm = width)
            m_apo[m_apo<0] = 0
            m_apo /= np.max(m_apo)
            i+=1
        return m_apo
    def get_window(mapk,n_it=5,width = 0.1):
        """Helper to apodize_map"""
        m_apo = np.copy(mapk)*0
        if hp.UNSEEN in np.copy(mapk):
            m_apo[np.copy(mapk)!=hp.UNSEEN] = 1.0
        else:
            m_apo[np.copy(mapk)!=0] = 1.0

        return smooth_map(m_apo, n_it=5, width = 0.1)
    def apodize_map(map0,n_it =5):
        """Apodizes a map with hp.UNSEEN pixels as masks. n_it represents the iterations of gaussian convolutions. Higher n_it generally means bigger windows and lower n_it the inverse."""

        tmp = np.copy(map0)
        tmp2 = np.copy(map0)
        tmp1 = tmp!=hp.UNSEEN
        m_apo = get_window(tmp2, n_it = n_it)
        tmp[tmp1!=True] =0
        output = tmp*m_apo
        return output

    chs = channels
    final = []

    for ch in chs:
        simulator = mapsims.MapSim(
            channels=ch,
            nside=NSIDE,
            unit="uK_CMB",
            pysm_output_reference_frame="C",
            pysm_components_string=pysm_components,
            pysm_custom_components={"cmb": cmb} if sim_cmb else None,
            other_components={"noise": noise} if sim_noise else None,
            num=seed,
        )
        output_map = simulator.execute()
        for det in output_map.keys():
            for pol in np.arange(output_map[det].shape[0]):
                output_map[det][pol] = apodize_map(output_map[det][pol])
        final.append(output_map)
    pols = ["T", "Q", "U"]
    for h in final:
        for k in h.keys():
            print(k)
            for pol in np.arange(h[k].shape[0]):
                hp.mollview(h[k][pol], title = str(k) + " " + pols[pol])
                plt.show()
def _main():
    i = getInputs(os.path.join(
        absolute_path, "input.yaml"))
    angle = 90 - i["observationElevationAngle"]
    coldSpillOverEfficiency = getSpillEfficiency(i)

    calculate = calcByAngle(i["diameter"], i["t"], i["wfe"], i["eta"], i["doe"], i["t_int"], i["pixelYield"], i["szCamNumPoln"], i["eorSpecNumPoln"],
                            i["t_filter_cold"], i["t_lens_cold"], i["t_uhdpe_window"], coldSpillOverEfficiency, i["singleModedAOmegaLambda2"],
                            i["spatialPixels"], i["eqbw"], i["centerFrequency"], i["detectorNEP"],
                            i["backgroundSubtractionDegradationFactor"], i["r"])

    outputs = calculate(angle)

    # valueDisplay = valDisplayPartial(
    #    i["outputFreq"], i["centerFrequency"], outputs["wavelength"], i["decimalPlaces"])
    # quartileDisplay = quartDisplayPartial(
    #    i["outputFreq"], i["centerFrequency"], outputs["wavelength"], i["decimalPlaces"])

    # useLatexFont()
    # outputSensitivityFile(outputs, valueDisplay, quartileDisplay)
    # outputPowerFile(i, outputs, calculate, quartileDisplay)
    # outputSpillEfficiencyFile(i, calculate, coldSpillOverEfficiency)
    # outputLoadings(i, calculate)

    noiseCurves = getNoiseCurves(i, outputs)
    seed = 0
    for pysm_components, sim_noise in zip(["d1", None], [False, True]):
        if pysm_components == "d1":
            continue
        ccat_mapsims(i, outputs, noiseCurves, ["tube:LC1"], pysm_components, seed, sim_cmb=False, sim_noise=sim_noise)
        so_mapsims(["tube:LT0"], pysm_components, seed, sim_cmb=False, sim_noise=sim_noise)
        ccat_mapsims(i, outputs, noiseCurves, ["tube:LC3"], pysm_components, seed, sim_cmb=False, sim_noise=sim_noise)
    """inputs = {'diameter': 5.7, 't': 273, 'wfe': 10.7, 'eta': 0.98, 'doe': 0.8, 'pixelYield': 0.8,
              'eorSpecNumPoln': 2, 't_filter_cold': np.array([1, 1]), 't_lens_cold': np.array([.98, .98]), 't_uhdpe_window': np.array([1, 1]), 'spatialPixels': np.array([3456, 3072]),
              'centerFrequency': np.array([262.5*10**9, 367.5*10**9]), 'detectorNEP': 0,
              'backgroundSubtractionDegradationFactor': 1, 'observationElevationAngle': 45, 'detectorSpacing': np.array([2.75, 2.09]), 'lyotStopAngle': 13.4}
    rfpairs = np.array([(101, 250*10**9), (102, 350*10**9),
                       (103, 275*10**9), (104, 100*10**9)])
    # print(eorNoiseCurves(inputs, rfpairs)[(101, 250*10**9)])

    # spillEfficiencyComparison(f=280e9, maxangle=180)
    """

if __name__ == "__main__":
    _main()
