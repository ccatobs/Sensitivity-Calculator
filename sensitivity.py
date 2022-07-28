import yaml
import numpy as np
from functools import partial
import matplotlib.pyplot as plt
from texttable import Texttable
import noise

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


def custOutput(i, outputs, actuallyCalculate=False):
    """Temporary function for playing with mapsims"""
    def data_C_calc():
        p = 50
        a = 45
        col = 2

        def averageTemp(prefix, angle, percentile, start, end):
            file = open("data/" + prefix + "/" + str(percentile) +
                        "/ACT_annual_" + str(percentile) + "." + str(angle) + ".out", "r")
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

        higher = np.array([averageTemp("Higher", a, p, c-w/2, c+w/2)
                          for c, w in zip(i["centerFrequency"], i["eqbw"])])
        lower = np.array([averageTemp("Lower", a, p, c-w/2, c+w/2)
                          for c, w in zip(i["centerFrequency"], i["eqbw"])])
        derivative = (higher - lower) / 0.02
        print(derivative*10.6/13.2)

        def A_to_CMB(freq_in_GHz):
            h = _h
            kb = _k
            T = 2.725
            v = freq_in_GHz*1e9
            x = h*v/(kb*T)
            return 1./(x**2*np.exp(x)/(np.exp(x)-1)**2)
        last = (_averageTrans("Higher/", a, p, 145e9, 145*0.276e9, col) -
                _averageTrans("Lower/", a, p, 145e9, 145*0.276e9, col))/0.02
        dataCs = np.array([])
        for n in range(len(i["centerFrequency"])):
            dataCs = np.append(dataCs,
                               (derivative[n]*A_to_CMB(i["centerFrequency"][n]/1e9)/(last*A_to_CMB(145)))**2)
        print(dataCs[::-1]*1.2e4)
        return dataCs[::-1]*1.2e4
    if True:
        centerFrequency = None
        beam = None
        net = None
        data_C = None
        if not actuallyCalculate:
            centerFrequency = [222., 280., 348., 405., 850.]
            beam = [59/60., 47/60., 37/60., 32/60., 15/60.]
            net = [6.8, 12.7, 47.7, 181.8, 305400.7]
            data_C = None  # Automatic values
        else:
            centerFrequency = i['centerFrequency'][::-1]/1e9
            beam = outputs["beam"][::-1]/60
            net = outputs["netW8Avg"][::-1]
            data_C = data_C_calc()
        ccat = noise.CCAT(centerFrequency, beam, net, survey_years=4000 /
                          24./365.24, survey_efficiency=1.0, N_tubes=(1, 1, 1, 1, 1), el=45., data_C=data_C)
        fsky = 20000./(4*_pi*(180/_pi)**2)
        lat_lmax = 10000
        ell, N_ell_T_full, N_ell_P_full = ccat.get_noise_curves(
            fsky, lat_lmax, 1, full_covar=False, deconv_beam=True)
        plotTemperature = True
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
    if False:
        p = 50
        a = 45
        col = 3
        higher = np.array([_averageTrans("Higher/", a, p, c, w, col)
                           for c, w in zip(i["centerFrequency"], i["eqbw"])])
        lower = np.array([_averageTrans("Lower/", a, p, c, w, col)
                          for c, w in zip(i["centerFrequency"], i["eqbw"])])
        print("Frequencies:")
        print((i["centerFrequency"]/1e9).astype(int)[::-1])
        #print("T slightly higher/lower for derivative use:")
        # print(higher)
        # print(lower)
        derivative = (higher - lower) / 0.02
        print("Slopes:")
        print(((derivative*10).astype(int)/10)[::-1])
        tB = derivative**2

        def A_to_CMB(freq_in_GHz):
            h = _h
            kb = _k
            T = 2.725
            v = freq_in_GHz*1e9
            x = h*v/(kb*T)
            return 1./(x**2*np.exp(x)/(np.exp(x)-1)**2)
        cmb = np.array([A_to_CMB(f/1e9) for f in i["centerFrequency"]])
        tTh = tB * cmb
        print("dT/dPWV:")
        print(tTh[::-1])
        print("data_C:")
        print(np.array([
            # below factors from am_output/mult_pwv/get_derivative_ccat.py
            2.31956542e+05,
            1.61527385e+06,
            4.03473727e+07,
            2.51490116e+08,
            9.10884821e+13
        ]))
        for s, f in zip(derivative[::-1], i["centerFrequency"][::-1]/1e9):
            print("freq=%.1f; slope = %.1f" % (f, s))

        last = (_averageTrans("Higher/", a, p, 145e9, 145*0.276e9, col) -
                _averageTrans("Lower/", a, p, 145e9, 145*0.276e9, col))/0.02
        dataCs = np.array([])
        for n in range(len(i["centerFrequency"])):
            dataCs = np.append(dataCs,
                               (derivative[n]*A_to_CMB(i["centerFrequency"][n]/1e9)/(last*A_to_CMB(145)))**2)
        print(dataCs[::-1]*1.2e4)


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
    custOutput(i, outputs, actuallyCalculate=True)
