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
    """Performs the calculations done in the original sensitivity calculations excel sheet, based on "CCAT-Prime_SZ-FPI_20220117.xlsx". Returns a dictionary of [netW8Avg], [netW8RJ], [neiW8], [eorNEFD], [eorNEI], [powerPerPixel], [eorPowerPerPixel], [wavelength] (corresponding to the center frequencies of the detectors), and [beam]."""
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
            powerPerPixel, "eorPowerPerPixel": eorPowerPerPixel, "wavelength": wavelength, "beam": beam,
            "eorPhotonNoiseNEP": eorPhotonNoiseNEP, "photonNoiseNEP": photonNoiseNEP, "EqTrans": eqtrans}


def _pwvToTick(pwv):
    """Internal function that converts from [pwv] to an integer representing a fraction of the median pwv. Each number represents 5% of median pwv. I.e. 20 represents 100% of median pwv."""
    return int(pwv/(ccatPWVQ2/20))


def _averageTransSE(filePath, start, end, col=1):
    """Returns the average transmission given start and end frequencies."""
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
    """Returns the average transmission given a center frequency and width."""
    return _averageTransSE(filePath, (center-width/2)/1e9, (center+width/2)/1e9, col)


def _averageTrans(prefix, angle, percentile, center, width, col=1, newFilePathFormat=False):
    """Returns the average transmission given a center frequency and width. Should be used over _averageTransHelper directly due to computing the file path for you and not throwing errors when giving impossible observation angles."""
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


def _getEQTransV2(angle, center, width):
    """Returns the transmission at a given center frequency and width, linearly interpolating between precomputed transmission data."""
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
    """Rounds the values in a 1d or 2d array to a given precision in [decimalPlaces]."""
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
    """Outputs a dictionary with nicely formatted information. valDisplayPartial should be called instead of this function."""
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


def _valueDisplay(array, outputFreq, centerFrequency, wavelength, decimalPlaces):
    """Creates a dictionary to be put into a yaml file for broadband data."""
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
    """Outputs a dictionary with nicely formatted information. quartDisplayPartial should be called instead of this function."""
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


def _quartileDisplay(array, outputFreq, centerFrequency, wavelength, decimalPlaces):
    """Creates a dictionary to be put into a yaml file for data involving quartiles"""
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
    """Calculates the spill efficiency of the telescope assuming that the data gets contracted proportional to [contractFactor], such as for various wavelengths and detector spacings."""
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
    """Returns a 1d array of spill efficiencies for corresponding wavelengths and detector spacings."""
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
def getSpillEfficiency(i, oldFile=True):
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


def _data_C_calcV2(i):
    """Returns data_C based on the input dictionary [i] specified by getInputs."""
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
    """Changes matplotlib to use formatting that looks like latex."""
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


def _eorCalculate(diameter, t, wfe, eta, doe, eorSpecNumPoln, t_filter_cold, t_lens_cold, t_uhdpe_window, coldSpillOverEfficiency, centerFrequency, detectorNEP, backgroundSubtractionDegradationFactor, eqtrans, spatialPixels, pixelYield, eorEqBw, detIndex, fnum):
    """Equivalent to _calculate but for EoR calculations. Should be used instead of the EoR outputs from the _calculate function."""
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
    return {"eorNEFD": eorNEFD, "eorNEI": eorNEI, "eorPowerPerPixel": eorPowerPerPixel, "wavelength": wavelength, "beam": beam, "netW8Avg": netw8avg, "eorPhotonNoiseNEP":eorPhotonNoiseNEP, "eorEqTrans": eorEqTrans}


def _eorCalcByFreqR(angle, diameter, t, wfe, eta, doe, eorSpecNumPoln, t_filter_cold, t_lens_cold, t_uhdpe_window, coldSpillOverEfficiency, detectorNEP, backgroundSubtractionDegradationFactor, spatialPixels, pixelYield, eqbw):
    """Returns a function that takes observation zenith angle as an input and returns the following outputs in a dictionary: EoR Spec NEFD as "eorNEFD", EoR Spec NEI as "eorNEI", Power per Pixel as "powerPerPixel", EoR Spec Power per Pixel as "eorPowerPerPixel", Center Wavelengths as "wavelength", Beam as "beam"."""
    partTrans = partial(_getEQTransV2, angle=angle)
    partCalc = partial(_eorCalculate, diameter=diameter, t=t, wfe=wfe, eta=eta, doe=doe, eorSpecNumPoln=eorSpecNumPoln,  t_filter_cold=t_filter_cold,  t_lens_cold=t_lens_cold,  t_uhdpe_window=t_uhdpe_window,
                       coldSpillOverEfficiency=coldSpillOverEfficiency, detectorNEP=detectorNEP,  backgroundSubtractionDegradationFactor=backgroundSubtractionDegradationFactor, spatialPixels=spatialPixels, pixelYield=pixelYield, eorEqBw=eqbw)
    return lambda f, detIndex, fnum: partCalc(centerFrequency=f, detIndex=detIndex, fnum=fnum, eqtrans=partTrans(center=f, width=eqbw))


def _geteqbw(f, r):
    """Returns the equivalent bandwith for a given frequency and r for EoR calculations."""
    #return f / (10 * r) * _pi / 2
    return f / r * _pi / 2


def eorNoiseCurves(i, rfpairs, frequencyRanges=np.array([[210, 315], [315, 420]])*10**9,fsky = 20000./(4*_pi*(180/_pi)**2),survey_years=4000/24./365.24,return_outputs=False):
    """Returns a list of noise curves (ell, N_ell_T_full, N_ell_P_full) corresponding to rfpairs, a list of finesse-frequency pairs, and detector arrays i. frequencyRanges is an optional parameter 
    for modifying the frequency range of the detector arrays."""
    angle = 90 - i["observationElevationAngle"]
    coldSpillOverEfficiency = getSpillEfficiency(i)

    i["eqbw"] = _geteqbw(rfpairs[:, 1], rfpairs[:, 0])
    calculate = _eorCalcByFreqR(angle, i["diameter"], i["t"], i["wfe"], i["eta"], i["doe"], i["eorSpecNumPoln"],
                                i["t_filter_cold"], i["t_lens_cold"], i["t_uhdpe_window"], coldSpillOverEfficiency,
                                i["detectorNEP"], i["backgroundSubtractionDegradationFactor"], i["spatialPixels"], i["pixelYield"], i["eqbw"])

    def _outsourcedCalculation(i, output, fsky):
        centerFrequency = i["centerFrequency"] / 1e9
        beam = output["beam"] / 60
        netw8avg = output["netW8Avg"]
        data_C = _data_C_calcV2(i)
        ccat = noise_file.CCAT(centerFrequency, beam, netw8avg, survey_years=survey_years, survey_efficiency=1.0, N_tubes=tuple(1 for _ in centerFrequency), el=45., data_C=data_C)
        lat_lmax = 10000
        ell, N_ell_T_full, N_ell_P_full = ccat.get_noise_curves(
            fsky, lat_lmax, 1, full_covar=False, deconv_beam=True)
        return ell, N_ell_T_full, N_ell_P_full

    results = {}
    outputs = {}
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
        ell, N_ell_T_full, N_ell_P_full = _outsourcedCalculation(i,output,fsky)
        N_ell_T_full = N_ell_T_full[index]
        N_ell_P_full = N_ell_P_full[index]
        results[(ri, fi)] = (ell, N_ell_T_full, N_ell_P_full)
        outputs[(ri, fi)] = output

    if return_outputs:
        return results, outputs
    else:
        return results


def spillEfficiencyComparison(lyotStopAngle=13.4, f=350e9, ds=2.75, maxangle=180):
    """Graphs the 280 and 350 GHz data for calculating the spill efficiency."""
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


def _smooth_map(m, n_it=5, width=0.1):
    """Helper to get_window"""
    i = 0
    m_apo = m
    while i <= n_it:
        m_apo[m_apo < 0.8] = 0
        m_apo = hp.smoothing(m_apo, fwhm=width)
        m_apo[m_apo < 0] = 0
        m_apo /= np.max(m_apo)
        i += 1
    return m_apo


def _get_window(mapk, n_it=5, width=0.1):
    """Helper to apodize_map"""
    m_apo = np.copy(mapk)*0
    if hp.UNSEEN in np.copy(mapk):
        m_apo[np.copy(mapk) != hp.UNSEEN] = 1.0
    else:
        m_apo[np.copy(mapk) != 0] = 1.0

    return _smooth_map(m_apo, n_it=5, width=0.1)


def _apodize_map(map0, n_it=5):
    """Apodizes a map with hp.UNSEEN pixels as masks. n_it represents the iterations of gaussian convolutions. Higher n_it generally means bigger windows and lower n_it the inverse."""

    tmp = np.copy(map0)
    tmp2 = np.copy(map0)
    tmp1 = tmp != hp.UNSEEN
    m_apo = _get_window(tmp2, n_it=n_it)
    tmp[tmp1 != True] = 0
    output = tmp*m_apo
    return output


def ccat_mapsims(i, outputs, band, tube, pysm_components, seed, data_C, sim_cmb=False, sim_noise=False, instrument_parameters_path="/Users/stevekchoi/work/build/Sensitivity-Calculator/src/sensitivity_calculator/data/instrument_parameters/instrument_parameters.tbl", hitmap_path="/Users/stevekchoi/work/build/Sensitivity-Calculator/src/sensitivity_calculator/data/ccat_uniform_coverage_nside256_201021.fits", NSIDE=256, lmax=None):
    """Graphs and returns the map corresponding to a given [band] and [tube] in CCAT, with 
    [pysm_components] and [seed] fed into mapsims to create the map. Sensitivities from the rest of 
    the calculator are passed in through the input parameters [i] and broadband outputs [outputs]. 
    [sim_cmb] and [sim_noise] are toggles for whether the CMB and noise will be simulated respectively."""
    if lmax is None:
        lmax = 3*NSIDE-1
    lmax_over_nside = int(lmax/NSIDE)
    instrument_path = Path(instrument_parameters_path)
    channels = ["tube:" + tube]
    tag = tube + "_" + band
    if sim_cmb:
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
    ccat_survey = noise_file.CCAT(
        i["centerFrequency"], outputs["beam"]/60., outputs["netW8Avg"], hitmap_path=hitmap_path, data_C=data_C)
    channels_list = mapsims.parse_channels(
        instrument_parameters=instrument_path)
    noise = mapsims.noise.ExternalNoiseSimulator(
        nside=NSIDE,
        return_uK_CMB=True,
        sensitivity_mode="baseline",
        apply_beam_correction=True,
        apply_kludge_correction=False,
        survey=ccat_survey,
        channels_list=channels_list,
        survey_efficiency=1.0,
        rolloff_ell=30,
    )

    chs = channels
    final = []

    for ch in chs:
        simulator = mapsims.MapSim(
            channels=ch,
            nside=NSIDE,
            unit="uK_CMB",
            lmax_over_nside=lmax_over_nside,
            output_reference_frame="C",
            pysm_components_string=pysm_components,
            pysm_custom_components={"cmb": cmb} if sim_cmb else None,
            other_components={"noise": noise} if sim_noise else None,
            instrument_parameters=instrument_path,
            num=seed
        )
        output_map = simulator.execute()
        for det in output_map.keys():
            for pol in np.arange(output_map[det].shape[0]):
                output_map[det][pol] = _apodize_map(output_map[det][pol])
        final.append(output_map)
    pols = ["T", "Q", "U"]
    for h in final:
        for k in h.keys():
            if k == tag:
                for pol in np.arange(h[k].shape[0]):
                    hp.mollview(h[k][pol], title=str(k) + " " + pols[pol])
                    plt.show()
            else:
                pass
    return final[0][tag]


def so_mapsims(band, tube, pysm_components, seed, sim_cmb=False, sim_noise=False, instrument_parameters_path="/Users/stevekchoi/work/build/Sensitivity-Calculator/src/sensitivity_calculator/data/instrument_parameters/instrument_parameters.tbl", hitmap="/Users/stevekchoi/work/build/Sensitivity-Calculator/src/sensitivity_calculator/data/ccat_uniform_coverage_nside256_201021.fits", NSIDE=256, lmax=None):
    """Graphs and returns the map corresponding to a given [band] and [tube] in SO's telescopes, with 
    [pysm_components] and [seed] fed into mapsims to create the map. [sim_cmb] and [sim_noise] are toggles
    for whether the CMB and noise will be simulated respectively."""
    if lmax is None:
        lmax = 3*NSIDE-1
    lmax_over_nside = int(lmax/NSIDE)
    instrument_path = Path(instrument_parameters_path)
    channels = ["tube:" + tube]
    tag = tube + "_" + band
    if sim_cmb:
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
        apply_kludge_correction=False,
        instrument_parameters=instrument_path,
    )

    chs = channels
    final = []

    for ch in chs:
        simulator = mapsims.MapSim(
            channels=ch,
            nside=NSIDE,
            unit="uK_CMB",
            lmax_over_nside=lmax_over_nside,
            output_reference_frame="C",
            pysm_components_string=pysm_components,
            pysm_custom_components={"cmb": cmb} if sim_cmb else None,
            other_components={"noise": noise} if sim_noise else None,
            instrument_parameters=instrument_path,
            num=seed,
        )
        output_map = simulator.execute()
        for det in output_map.keys():
            for pol in np.arange(output_map[det].shape[0]):
                output_map[det][pol] = _apodize_map(output_map[det][pol])
        final.append(output_map)
    pols = ["T", "Q", "U"]
    for h in final:
        for k in h.keys():
            if k == tag:
                for pol in np.arange(h[k].shape[0]):
                    hp.mollview(h[k][pol], title=str(k) + " " + pols[pol])
                    plt.show()
            else:
                pass
    return final[0][tag]


def zeroHitmapFraction(path, nside):
    """Returns the fraction of zeros in a given hitmap given a path to the hitmap [path] and the nside [nside]."""
    hitmap = hp.ud_grade(
        hp.read_map(path, dtype=np.float64),
        nside_out=nside,
    )
    zeros = 0
    total = 0
    hitmap = np.array(hitmap)
    for i in hitmap:
        if i == 0:
            zeros += 1
        total += 1
    return zeros / total


def plotPowerSpectrum(tt, ee, bb, title, zf=1,TT_theory=None,PP_theory=None):
    """Plots the power spectrum. zf is the fraction of the hitmap with 0s"""
    plt.plot([i for i in range(len(tt))], tt/zf, linewidth=1, label="tt")
    plt.plot([i for i in range(len(ee))], ee/zf, linewidth=1, label="ee")
    plt.plot([i for i in range(len(bb))], bb/zf, linewidth=1, label="bb")
    if TT_theory is not None:
        plt.plot(TT_theory)
    if PP_theory is not None:
        plt.plot(PP_theory)
    plt.yscale("log")
    plt.ylabel("$C_{\ell}$")
    plt.xscale("log")
    plt.xlabel("$\ell$")
    plt.xlim(10**2, 10**4)
    plt.title("Power Spectrum of " + title)
    plt.legend()
    plt.show()


def _main():
    """Basic testing of functions"""
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
    NSIDE = 256
    hitmap_path = "/Users/stevekchoi/work/build/Sensitivity-Calculator/src/sensitivity_calculator/ccat_uniform_coverage_nside" + \
        str(NSIDE) + "_201021.fits"
    zeroHitmapFractio = zeroHitmapFraction(hitmap_path, NSIDE)
    zf = zeroHitmapFractio
    data_C = _data_C_calcV2(i)
    hitmap = hp.ud_grade(
        hp.read_map(hitmap_path, dtype=np.float64),
        nside_out=NSIDE,
    )
    for pysm_components, sim_noise in zip(["d2", None], [False, True]):
        if pysm_components == "d2":
            continue
        ccat280 = ccat_mapsims(
            i, outputs, "HF2", "LC1", pysm_components, seed, data_C, sim_cmb=False, sim_noise=sim_noise, hitmap_path=hitmap_path, NSIDE=NSIDE)
        ccat280cls = hp.sphtfunc.anafast(ccat280)
        plotPowerSpectrum(ccat280cls[0], ccat280cls[1], ccat280cls[2],
                          "280 GHz CCAT " + ("Noise" if sim_noise else "Signal"), zf=zf)
        so280 = so_mapsims("UHF2", "LT0", pysm_components,
                           seed, sim_cmb=False, sim_noise=sim_noise, NSIDE=NSIDE, hitmap=hitmap)
        so280cls = hp.sphtfunc.anafast(so280)
        plotPowerSpectrum(so280cls[0], so280cls[1], so280cls[2],
                          "280 GHz SO " + ("Noise" if sim_noise else "Signal"), zf=zf)
        ccat850 = ccat_mapsims(
            i, outputs, "HF5", "LC3", pysm_components, seed, data_C, sim_cmb=False, sim_noise=sim_noise, hitmap_path=hitmap_path, NSIDE=NSIDE)
        ccat850cls = hp.sphtfunc.anafast(ccat850)
        plotPowerSpectrum(ccat850cls[0], ccat850cls[1], ccat850cls[2],
                          "850 GHz CCAT " + ("Noise" if sim_noise else "Signal"), zf=zf)
    # ccat_mapsims(i, outputs, "HF2", "LC1",
        # "d1", seed, sim_cmb=False, sim_noise=True)
    # ccat_mapsims(i, outputs, "HF5", "LC3",
        # "d1", seed, sim_cmb=False, sim_noise=True)
    # inputs = {'diameter': 5.7, 't': 273, 'wfe': 10.7, 'eta': 0.98, 'doe': 0.8, 'pixelYield': 0.8,
    #          'eorSpecNumPoln': 2, 't_filter_cold': np.array([1, 1]), 't_lens_cold': np.array([.98, .98]), 't_uhdpe_window': np.array([1, 1]), 'spatialPixels': np.array([3456, 3072]),
    #          'centerFrequency': np.array([262.5*10**9, 367.5*10**9]), 'detectorNEP': 0,
    #          'backgroundSubtractionDegradationFactor': 1, 'observationElevationAngle': 45, 'detectorSpacing': np.array([2.75, 2.09]), 'lyotStopAngle': 13.4}
    # rfpairs = np.array([(101, 250*10**9), (102, 350*10**9),
    #                   (103, 275*10**9), (104, 100*10**9)])
    # print(eorNoiseCurves(inputs, rfpairs)[(101, 250*10**9)])

    # spillEfficiencyComparison(f=280e9, maxangle=180)


if __name__ == "__main__":
    _main()
