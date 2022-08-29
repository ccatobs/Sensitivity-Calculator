from sklearn.metrics import plot_det_curve
import yaml
import numpy as np
from functools import partial
import matplotlib.pyplot as plt
from texttable import Texttable
import noise
import scipy.optimize as op
import pwvCalculator as pwv

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


def _calculate(diameter, t, wfe, eta, doe, t_int, pixelYield, szCamNumPoln, eorSpecNumPoln, t_filter_cold, t_lens_cold, t_uhdpe_window, coldSpillOverEfficiency, singleModedAOmegaLambda2, spatialPixels, fpi, eqbw, centerFrequency, detectorNEP, backgroundSubtractionDegradationFactor, sensitivity, hoursPerYear, sensPerBeam, r, signal, eqtrans):
    wavelength = _c/centerFrequency*10**6

    def a_to_CMB(f):
        kb = 1.3806488e-23
        T = 2.725
        v = f*1e9
        x = _h*v/(kb*T)
        return 1./(x**2*np.exp(x)/(np.exp(x)-1)**2)

    aToCMB = np.array([a_to_CMB(i/1e9) for i in centerFrequency])

    # Telescope
    a = _pi*(diameter/2)**2

    # Detector
    szCamCorForPoln = abs(3-szCamNumPoln)
    eorSpecCorForPoln = abs(3-eorSpecNumPoln)

    #Instrument, throughput
    t_cold = t_filter_cold*t_lens_cold**3
    e_window_warm = 1 - t_uhdpe_window

    # Instrument, beams/area of focal plane
    # WindowTrans's formula doesn't make sense but is the same as the sheet, and luckily is not used here nor in the sheet
    windowTrans = np.ones(len(wavelength)) * t_uhdpe_window[1]
    beam = 1.2*wavelength/diameter/1000000*206265
    solidAngle = _pi/4/_ln(2)*(beam/206264)**2
    nativeAOmegaLambda2 = solidAngle*a/(wavelength*0.000001)**2
    fov = _pi*0.45**2
    hornFov = fov/spatialPixels
    hornDiameter = _sqrt(4*hornFov/_pi)*3600
    hornSolidAngle = _pi*(hornDiameter/2/206265)**2
    beamSolidAngle = _pi/_ln(2)*(beam/2/206265)**2
    beamsPerHorn = hornSolidAngle/beamSolidAngle
    effectiveFovInSr = solidAngle*spatialPixels
    effectiveFov = (180/_pi)**2*effectiveFovInSr
    fovFillingFactor = fov/effectiveFov
    pointingsIn1000SqDeg = 1000/effectiveFov
    secondsPerPointingIn4000Hrs = 4000*3600/pointingsIn1000SqDeg
    noiseFactor1000SqDeg4000Hrs = _sqrt(secondsPerPointingIn4000Hrs)

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


def _averageTransSE(filePath, start, end, col=1):
    file = open("data/" + filePath + ".out", "r")
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


def _averageTrans(prefix, angle, percentile, center, width, col=1):
    if angle >= 15 and angle <= 75 and int(angle) == angle:
        return _averageTransHelper(prefix + str(percentile) + "/ACT_annual_" + str(percentile) + "." + str(angle), center, width, col)
    elif int(angle) != angle:
        floor = _averageTransHelper(prefix + str(percentile) + "/ACT_annual_" +
                                    str(percentile) + "." + str(int(np.floor(angle))), center, width, col)
        ceil = _averageTransHelper(prefix + str(percentile) + "/ACT_annual_" +
                                   str(percentile) + "." + str(int(np.ceil(angle))), center, width, col)
        prop = angle - np.floor(angle)
        return floor * (1 - prop) + ceil * prop
    else:
        print("Angle out of range")


def _getEQTrans(angle, center, width):
    return np.array([[_averageTrans("CerroConfig/", angle, percentile, centeri, widthi) for percentile in [25, 50, 75]] for centeri, widthi in zip(center, width)])


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


def calcByAngle(diameter, t, wfe, eta, doe, t_int, pixelYield, szCamNumPoln, eorSpecNumPoln, t_filter_cold, t_lens_cold, t_uhdpe_window, coldSpillOverEfficiency, singleModedAOmegaLambda2, spatialPixels, fpi, eqbw, centerFrequency, detectorNEP, backgroundSubtractionDegradationFactor, sensitivity, hoursPerYear, sensPerBeam, r, signal):
    """Returns a function that takes observation zenith angle as an input and returns the following outputs in a dictionary: NET Weighted Average as "netW8Avg", NET Weighted RJ as "netW8RJ", NEI Weighted Average as "neiW8", EoR Spec NEFD as "eorNEFD", EoR Spec NEI as "eorNEI", Power per Pixel as "powerPerPixel", EoR Spec Power per Pixel as "eorPowerPerPixel", Center Wavelengths as "wavelength", Beam as "beam"."""

    partTrans = partial(_getEQTrans, center=centerFrequency, width=eqbw)
    partCalc = partial(_calculate, diameter, t, wfe, eta, doe, t_int, pixelYield, szCamNumPoln, eorSpecNumPoln, t_filter_cold, t_lens_cold, t_uhdpe_window, coldSpillOverEfficiency,
                       singleModedAOmegaLambda2, spatialPixels, fpi, eqbw, centerFrequency, detectorNEP, backgroundSubtractionDegradationFactor, sensitivity, hoursPerYear, sensPerBeam, r, signal)
    return lambda x: partCalc(partTrans(x))


def sensitivityFile(outputs, valueDisplay, quartileDisplay):
    """Outputs output.yaml"""
    dict_file = {"NET w8 avg": valueDisplay(outputs["netW8Avg"]),
                 "NET w8 RJ": valueDisplay(outputs["netW8RJ"]), "NEI w8 Jy/sr": valueDisplay(
        outputs["neiW8"]), "EoR Spec NEFD": quartileDisplay(outputs["eorNEFD"]), "EoR Spec NEI": quartileDisplay(outputs["eorNEI"])}
    return yaml.dump(dict_file, open("output.yaml", 'w'), sort_keys=False)


def powerFile(outputs, calculate, quartileDisplay):
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
        plt.axvline(x=half_angle, color='_k',
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


def spillEfficiencyFile(i, calculate, spillEfficiency):
    """Prints a chart comparing calculated spill efficiencies to the excel values of spill efficiencies to stdout."""
    output = calculate(45)
    t = Texttable(max_width=110)
    t.set_cols_dtype(['i', (lambda x: "%.1f" % float(x)), 'f', 'i',
                     (lambda x: "%.2f" % float(x)), (lambda x: "%.1f" % float(x)), 'f'])
    t.set_cols_align(['_c', '_c', '_c', '_c', '_c', '_c', '_c'])
    excelSpillEff = [.8, .5, .7, .5, .5]
    excelNET = [241440.5, 181.8, 56.3, 11.4, 6.8]
    t.add_rows(np.concatenate((np.reshape(['Center Frequency (GHz)', 'Excel Spill Efficiency', 'Calculated Spill Efficiency', '# Pixels', 'Pixel Spacing (mm)', 'Excel NET (uK rt(s))', 'Calculated NET (uK rt(s))'], (-1, 1)),
                               np.array([[int(f/1e9), exSpillEf, spillEf, numDetect, detectSpacing, exNet, net] for f, spillEf, net, exSpillEf, numDetect, detectSpacing, exNet in zip(i['centerFrequency'], spillEfficiency, output['netW8Avg'], excelSpillEff, i["spatialPixels"], i['detectorSpacing'], excelNET)])[::-1].T), axis=1).T)
    print(t.draw())


# Load the beam file and pass the angle and value data to other functions
def getSpillEfficiency(i):
    """Returns spill efficiency, using data/tolTEC_staircase_singleHorn_280GHz.txt as a reference. Is meant to be updated when better/more curves are calculated."""
    data = np.genfromtxt(
        'data/tolTEC_staircase_singleHorn_280GHz.txt', skip_header=2).reshape(-1, 721, 8)
    degr = data[0, :, 0]
    vals = data[0, :, 3]
    #data = np.genfromtxt('data/beam_280.txt')
    #degr = np.degrees(data[:, 0])
    #vals = np.log10((data[:, 1]**2))*10
    return _getColdSpillOverEfficiency(
        i, 280e9, 2.75, degr, vals, showPlots=False)


def _averageTemp(start, end, col, filePath=None, prefix=None, angle=None, percentile=None):
    if filePath is None:
        filePath = prefix + "/" + \
            str(percentile) + "/ACT_annual_" + \
            str(percentile) + "." + str(angle)
    file = open("data/" + filePath + ".out", "r")
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


def _tangent_line_slope(cf, eqbw, col, a, p, maunaKea=False):
    higher = ""
    lower = ""
    if not maunaKea:
        higher = "Higher"
        lower = "Lower"
    else:
        higher = "MaunaKea/Higher"
        lower = "MaunaKea/Lower"
    # V1 Derivatives
    higher = np.array([_averageTemp(c-w/2, c+w/2, col, prefix=higher, angle=a, percentile=p)
                       for c, w in zip(cf, eqbw)])
    lower = np.array([_averageTemp(c-w/2, c+w/2, col, prefix=lower, angle=a, percentile=p)
                      for c, w in zip(cf, eqbw)])
    derivative = (higher - lower) / 0.02
    return derivative


def _least_squares_slope(cf, eqbw, col, a, graph=False, maunaKea=False):
    # V2 Derivatives
    filePath = "VariablePWV/ACT_annual_"
    if maunaKea:
        filePath = "MaunaKea/VariablePWV/"
    temps = np.array([[_averageTemp(c-w/2, c+w/2, col, filePath=(filePath + str(i) + "." + str(a)))
                       for c, w in zip(cf, eqbw)] for i in np.array(range(40))+1])

    def line(x, m, b):
        return m*x + b
    derivative = []
    ccatMedPWV = 0.67/pwv.configPWV(50)
    if maunaKea:
        ccatMedPWV = pwv.configPWVHelper("data/MaunaKea/VariablePWV/.err")
        print(ccatMedPWV)
    pwvs = (np.array(range(40))+1) / 20*ccatMedPWV
    for i in range(len(cf)):
        popt = op.curve_fit(line, pwvs, temps[:, i])[0]
        # print(popt)
        derivative.append(popt[0])
    derivative = np.array(derivative)
    if graph:
        if not maunaKea:
            print("Steve's PWV range:", (0.3*.7192506910, 3*.7192506910))
            print("CCAT 50th percentile PWV:", ccatMedPWV)
        for i in np.array(range(len(cf)))[::-1]:
            plt.plot(pwvs, temps[:, i], linewidth=1,
                     label=str(int(cf[i]/1e9))+' GHz')
        plt.ylim(bottom=0)
        plt.ylabel("Weighted Temperature (arbitrary unit)")
        plt.xlim(left=0, right=ccatMedPWV*2)
        plt.xlabel("PWV (mm)")
        plt.legend(loc='best')
        plt.grid()
        plt.show()
    return derivative


def _data_C_calc(i, table=False, graphSlopes=False, maunaKea=False):
    P = 50
    A = 45
    COL = 2

    cf = np.append(145e9, i["centerFrequency"])
    eqbw = np.append(145*0.276e9, i["eqbw"])

    def A_to_CMB(freq_in_GHz):
        h = _h
        kb = _k
        T = 2.725
        v = freq_in_GHz*1e9
        x = h*v/(kb*T)
        return 1./(x**2*np.exp(x)/(np.exp(x)-1)**2)
    dataCs = np.array([])
    # Choose method
    derivative = _tangent_line_slope(cf, eqbw, COL, A, P, maunaKea=maunaKea)
    for n in range(1, len(cf)):
        dataCs = np.append(dataCs,
                           (derivative[n]*A_to_CMB(cf[n]/1e9)/(derivative[0]*A_to_CMB(cf[0]/1e9)))**2)
    if table:
        t = Texttable(max_width=110)
        table_header = np.append("Method", np.char.add(
            (cf/1e9).astype(int).astype(str), ' GHz'))
        table = np.array([table_header, np.append("Tangent Line", _tangent_line_slope(cf, eqbw, COL, A, P, maunaKea=maunaKea)), np.append("Least Squares Regression Line", _least_squares_slope(cf, eqbw, COL, A, graph=graphSlopes, maunaKea=maunaKea)), np.append("Corrected PWV Previous Method",
                                                                                                                                                                                                                                                                    np.array([6.2, 14.7, 25.0, 48.5, 66.4, 64.4])), np.append("Uncorrected PWV Previous Method",
                                                                                                                                                                                                                                                                                                                              np.array([4.5, 10.6, 17.9, 34.9, 47.8, 46.3]))])
        t.add_rows(table, header=True)
        print(t.draw())
    elif graphSlopes:
        _least_squares_slope(cf, eqbw, COL, A, graph=True, maunaKea=maunaKea)
    return dataCs*1.2e4


def custOutput(i, outputs, calculate='all', plotCurve=None, table=False, graphSlopes=False, maunaKea=False):
    """Temporary function for playing with mapsims"""
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
    ccat = noise.CCAT(centerFrequency, beam, net, survey_years=4000 /
                      24./365.24, survey_efficiency=1.0, N_tubes=(1, 1, 1, 1, 1), el=45., data_C=data_C)
    fsky = 20000./(4*_pi*(180/_pi)**2)
    lat_lmax = 10000
    ell, N_ell_T_full, N_ell_P_full = ccat.get_noise_curves(
        fsky, lat_lmax, 1, full_covar=False, deconv_beam=True)
    if plotCurve is not None:
        plotTemperature = plotCurve == 'T' or plotCurve == 't' or plotCurve == 'temp' or plotCurve == 'temperature'
        # [:-1] removes 850 GHz for plotting
        for curve, label in zip(N_ell_T_full[:-1] if plotTemperature else N_ell_P_full[:-1], centerFrequency[:-1]):
            plt.plot(ell, curve, label=str(int(label))+' GHz')
            plt.yscale('log')
            plt.ylim(10**-5, 10**3)
            plt.xscale('log')
            plt.xlim(10**2, 10**4)
            plt.title("Temperature" if plotTemperature else "Polarization")
        plt.legend(loc='upper right')
        plt.grid()
        plt.show()


if __name__ == "__main__":
    i = getInputs("input.yaml")
    angle = 90 - i["observationElevationAngle"]
    coldSpillOverEfficiency = getSpillEfficiency(i)

    calculate = calcByAngle(i["diameter"], i["t"], i["wfe"], i["eta"], i["doe"], i["t_int"], i["pixelYield"], i["szCamNumPoln"], i["eorSpecNumPoln"],
                            i["t_filter_cold"], i["t_lens_cold"], i["t_uhdpe_window"], coldSpillOverEfficiency, i["singleModedAOmegaLambda2"],
                            i["spatialPixels"], i["fpi"], i["eqbw"], i["centerFrequency"], i["detectorNEP"],
                            i["backgroundSubtractionDegradationFactor"], i["sensitivity"], i["hoursPerYear"], i["sensPerBeam"], i["r"], i["signal"])

    outputs = calculate(angle)

    valueDisplay = valDisplayPartial(
        i["outputFreq"], i["centerFrequency"], outputs["wavelength"], i["decimalPlaces"])
    quartileDisplay = quartDisplayPartial(
        i["outputFreq"], i["centerFrequency"], outputs["wavelength"], i["decimalPlaces"])

    #sensitivityFile(outputs, valueDisplay, quartileDisplay)
    #powerFile(outputs, calculate, quartileDisplay)
    #spillEfficiencyFile(i, calculate, coldSpillOverEfficiency)
    custOutput(i, outputs, calculate='all', plotCurve=None,
               table=True, graphSlopes=True, maunaKea=True)
