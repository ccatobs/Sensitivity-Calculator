import yaml
import numpy as np
from functools import partial

# These functions work on numpy arrays (componentwise) and help with code clarity
pi = np.pi
ln = np.log
sqrt = np.sqrt
e = np.e
h = 6.626e-34
k = 1.38e-23
c = 299792458
arrayify = np.ndarray.tolist


def getInputs(filePath):
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
    eqtrans = None
    centerFrequency = None
    detectorNEP = None
    backgroundSubtractionDegradationFactor = None
    sensitivity = None
    hoursPerYear = None
    sensPerBeam = None

    r = None
    signal = None

    decimalPlaces = None
    calculateTrans = None
    angle = None
    outputFreq = None

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
        if key == "eqtrans":
            eqtrans = np.array(value)
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
        if key == "calculateTrans":
            calculateTrans = value
        if key == "angle":
            angle = value
        if key == "outputFreq":
            outputFreq = value

    stream.close()
    return {"diameter": diameter, "t": t, "wfe": wfe, "eta": eta, "doe": doe, "t_int": t_int, "pixelYield": pixelYield, "szCamNumPoln": szCamNumPoln, "eorSpecNumPoln": eorSpecNumPoln, "t_filter_cold": t_filter_cold, "t_lens_cold": t_lens_cold, "t_uhdpe_window": t_uhdpe_window, "coldSpillOverEfficiency": coldSpillOverEfficiency, "singleModedAOmegaLambda2": singleModedAOmegaLambda2, "spatialPixels": spatialPixels, "fpi": fpi, "eqbw": eqbw, "eqtrans": eqtrans, "centerFrequency": centerFrequency, "detectorNEP": detectorNEP, "backgroundSubtractionDegradationFactor": backgroundSubtractionDegradationFactor, "sensitivity": sensitivity, "hoursPerYear": hoursPerYear, "sensPerBeam": sensPerBeam, "r": r, "signal": signal, "decimalPlaces": decimalPlaces, "calculateTrans": calculateTrans, "angle": angle, "outputFreq": outputFreq}


def calculate(diameter, t, wfe, eta, doe, t_int, pixelYield, szCamNumPoln, eorSpecNumPoln, t_filter_cold, t_lens_cold, t_uhdpe_window, coldSpillOverEfficiency, singleModedAOmegaLambda2, spatialPixels, fpi, eqbw, centerFrequency, detectorNEP, backgroundSubtractionDegradationFactor, sensitivity, hoursPerYear, sensPerBeam, r, signal, eqtrans):
    wavelength = c/centerFrequency*10**6

    def a_to_CMB(f):
        kb = 1.3806488e-23
        T = 2.725
        v = f*1e9
        x = h*v/(kb*T)
        return 1./(x**2*np.exp(x)/(np.exp(x)-1)**2)

    aToCMB = np.array([a_to_CMB(i/1e9) for i in centerFrequency])

    # Telescope
    a = pi*(diameter/2)**2

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
    solidAngle = pi/4/ln(2)*(beam/206264)**2
    nativeAOmegaLambda2 = solidAngle*a/(wavelength*0.000001)**2
    fov = pi*0.45**2
    hornFov = fov/spatialPixels
    hornDiameter = sqrt(4*hornFov/pi)*3600
    hornSolidAngle = pi*(hornDiameter/2/206265)**2
    beamSolidAngle = pi/ln(2)*(beam/2/206265)**2
    beamsPerHorn = hornSolidAngle/beamSolidAngle
    effectiveFovInSr = solidAngle*spatialPixels
    effectiveFov = (180/pi)**2*effectiveFovInSr
    fovFillingFactor = fov/effectiveFov
    pointingsIn1000SqDeg = 1000/effectiveFov
    secondsPerPointingIn4000Hrs = 4000*3600/pointingsIn1000SqDeg
    noiseFactor1000SqDeg4000Hrs = sqrt(secondsPerPointingIn4000Hrs)

    # Weather quartiles/broadband
    e_warm = (t_uhdpe_window[:, None])*((1-eqtrans)
                                        * eta+(1-eta))+(e_window_warm[:, None])
    ruze = 1/(e**((4*pi*wfe/(wavelength))**2))
    weather_t_cold = (t_cold)
    occupancy = np.ones((len(wavelength), 3)) * 1 / \
        (e**(h*c/((wavelength)*10**(-6)*k*t))-1)[:, None]
    acceptedModes = szCamNumPoln * \
        (coldSpillOverEfficiency)*(singleModedAOmegaLambda2)
    powerPerPixel = h*c*e_warm*(weather_t_cold[:, None])*(acceptedModes[:, None])*occupancy*(eqbw[:, None])*doe/(
        ((wavelength)*10**(-6))[:, None])  # Differs from excel sheet becauses of correctly using t(cold)
    photonNoiseNEP = h*c/((wavelength)[:, None]*10**(-6))*(1/t_int*acceptedModes[:, None]*eqbw[:, None]
                                                           * e_warm*weather_t_cold[:, None]*doe*occupancy*(1+e_warm*weather_t_cold[:, None]*doe*occupancy))**0.5
    totalNEP = sqrt(photonNoiseNEP**2 + detectorNEP**2)
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
    netW8Avg = sqrt(3/(1/arrayNETCMB[:, 0]**2+1/arrayNETCMB[:, 1]
                    ** 2+1/arrayNETCMB[:, 2]**2).astype(float))
    netW8RJ = netW8Avg/aToCMB
    neiW8 = sqrt(3/(1/arrayNEI[:, 0]**2+1/arrayNEI[:, 1]
                    ** 2+1/arrayNEI[:, 2]**2).astype(float))
    nefd = arrayNEI*(solidAngle)[:, None] * \
        sqrt((spatialPixels)[:, None]*pixelYield)*1000

    # EoR Spec
    eorEqBw = c/(wavelength[:, None]*10**(-6)*r)*pi/2
    eorEqTrans = eqtrans
    eorE_warm = (t_uhdpe_window)[:, None]*((1-eorEqTrans)
                                           * eta+(1-eta))+(e_window_warm)[:, None]
    eorRuze = 1/e**((4*pi*wfe/(wavelength))**2)
    eorT_cold = (t_cold)*0.9  # Always takes from 739 um in the sheet
    eorOccupancy = 1/(e**(h*c/((wavelength)*10**(-6)*k*t))-1)
    eorAcceptedModes = eorSpecNumPoln*a*(coldSpillOverEfficiency)*(solidAngle)/(
        (wavelength)*10**(-6))**2  # In sheet 1071 um is calculated incorrectly
    eorPhotonNoiseNEP = h*c/((wavelength)[:, None]*10**(-6))*(eorAcceptedModes[:, None]*eorEqBw*eorE_warm *
                                                              eorT_cold[:, None]*doe*eorOccupancy[:, None]*(1+eorE_warm*eorT_cold[:, None]*doe*eorOccupancy[:, None]))**0.5
    eorTotalNEP = sqrt(eorPhotonNoiseNEP**2 + detectorNEP**2)
    eorColdTerminatedSpillover = (coldSpillOverEfficiency)
    eorNEF = eorSpecCorForPoln*eorTotalNEP/(a*eorT_cold[:, None]*doe*eta*eorRuze[:, None]
                                            * eorEqTrans*eorColdTerminatedSpillover[:, None]*backgroundSubtractionDegradationFactor)
    eorNEFD = eorNEF/eorEqBw*10**26*1000
    eorNEI = eorNEFD/(solidAngle)[:, None]/1000
    eorPowerPerPixel = h*c/(wavelength[:, None]*10**(-6))*eorE_warm*eorT_cold[:, None]*eorAcceptedModes[:,
                                                                                                        None]*eorOccupancy[:, None]*eorEqBw*doe  # Look at t_cold in excel sheet

    return {"netW8Avg": netW8Avg, "netW8RJ":
            netW8RJ, "neiW8": neiW8, "eorNEFD": eorNEFD, "eorNEI": eorNEI, "powerPerPixel":
            powerPerPixel, "eorPowerPerPixel": eorPowerPerPixel, "wavelength": wavelength}


def averageTransSE(filePath, start, end):
    file = open("data/" + filePath + ".out", "r")
    data = file.readlines()
    x = [float(i.split(" ")[0]) for i in data]
    y = [float(i.split(" ")[1]) for i in data]
    file.close()

    transmission = 0
    for (f, trans) in zip(x, y):
        if f < start:
            continue
        if f >= end:
            break
        transmission += trans

    return transmission / ((end - start) * 100)


def averageTransHelper(filePath, center, width):
    return averageTransSE(filePath, (center-width/2)/1e9, (center+width/2)/1e9)


def averageTrans(prefix, angle, percentile, center, width):
    if angle >= 15 and angle <= 75 and int(angle) == angle:
        return averageTransHelper(prefix + str(percentile) + "/ACT_annual_" + str(percentile) + "." + str(angle), center, width)
    elif int(angle) != angle:
        floor = averageTransHelper(prefix + str(percentile) + "/ACT_annual_" +
                                   str(percentile) + "." + str(int(np.floor(angle))), center, width)
        ceil = averageTransHelper(prefix + str(percentile) + "/ACT_annual_" +
                                  str(percentile) + "." + str(int(np.ceil(angle))), center, width)
        prop = angle - np.floor(angle)
        return floor * (1 - prop) + ceil * prop
    else:
        print("Angle out of range")


def getEQTrans(angle, center, width):
    return np.array([[averageTrans("CerroConfig/", angle, percentile, centeri, widthi) for percentile in [25, 50, 75]] for centeri, widthi in zip(center, width)])


def trun(array, decimalPlaces):
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


def valueDisplayHelper(array, w, dict, outputFreq, decimalPlaces):
    if len(w) > 0:
        unit = None
        if outputFreq:
            unit = " GHz"
        else:
            unit = " um"
        dict.update({str(int(w[0])) + unit: array[0]})
        return valueDisplayHelper(array[1:], w[1:], dict, outputFreq, decimalPlaces)
    else:
        return trun(dict, decimalPlaces)


# Creates a dictionary to be put into a yaml file for broadband data
def valueDisplay(array, outputFreq, centerFrequency, wavelength, decimalPlaces):
    if type(array) == np.ndarray:
        array = arrayify(array)
    outputLabel = None
    if outputFreq:
        outputLabel = centerFrequency / 1e9
    else:
        outputLabel = wavelength
    return valueDisplayHelper(array, outputLabel, {}, outputFreq, decimalPlaces)


def valDisplayPartial(outputFreq, centerFrequency, wavelength, decimalPlaces):
    return partial(valueDisplay, outputFreq=outputFreq, centerFrequency=centerFrequency, wavelength=wavelength, decimalPlaces=decimalPlaces)


def quartileDisplayHelper(array, w, dict, outputFreq, decimalPlaces):
    if len(w) > 0:
        unit = None
        if outputFreq:
            unit = " GHz"
        else:
            unit = " um"
        dict.update({str(int(w[0])) + unit: {"Quartile 1": array[0][0], "Quartile 2": array[0]
                    [1], "Quartile 3": array[0][2]}})
        return quartileDisplayHelper(array[1:], w[1:], dict, outputFreq, decimalPlaces)
    else:
        return trun(dict, decimalPlaces)


# Creates a dictionary to be put into a yaml file for data involving quartiles
def quartileDisplay(array, outputFreq, centerFrequency, wavelength, decimalPlaces):
    if type(array) == np.ndarray:
        array = arrayify(array)
    outputLabel = None
    if outputFreq:
        outputLabel = centerFrequency / 1e9
    else:
        outputLabel = wavelength
    return quartileDisplayHelper(array, outputLabel, {}, outputFreq, decimalPlaces)


def quartDisplayPartial(outputFreq, centerFrequency, wavelength, decimalPlaces):
    return partial(quartileDisplay, outputFreq=outputFreq, centerFrequency=centerFrequency, wavelength=wavelength, decimalPlaces=decimalPlaces)


def calcByAngle(diameter, t, wfe, eta, doe, t_int, pixelYield, szCamNumPoln, eorSpecNumPoln, t_filter_cold, t_lens_cold, t_uhdpe_window, coldSpillOverEfficiency, singleModedAOmegaLambda2, spatialPixels, fpi, eqbw, centerFrequency, detectorNEP, backgroundSubtractionDegradationFactor, sensitivity, hoursPerYear, sensPerBeam, r, signal):
    partTrans = partial(getEQTrans, center=centerFrequency, width=eqbw)
    partCalc = partial(calculate, diameter, t, wfe, eta, doe, t_int, pixelYield, szCamNumPoln, eorSpecNumPoln, t_filter_cold, t_lens_cold, t_uhdpe_window, coldSpillOverEfficiency,
                       singleModedAOmegaLambda2, spatialPixels, fpi, eqbw, centerFrequency, detectorNEP, backgroundSubtractionDegradationFactor, sensitivity, hoursPerYear, sensPerBeam, r, signal)
    return lambda x: partCalc(partTrans(x))


def methodsComparisonFile(i, quartileDisplay):
    def customOutput(angle):
        actSiteTrans = np.array([[averageTrans("", angle, 25, cent, wid), averageTrans("", angle, 50, cent, wid),
                                averageTrans("", angle, 75, cent, wid)] for (cent, wid) in zip(i["centerFrequency"], i["eqbw"])])
        steveTrans = np.array([[averageTrans("Steve/", angle, 25, cent, wid), averageTrans("Steve/", angle, 50, cent, wid),
                                averageTrans("Steve/", angle, 75, cent, wid)] for (cent, wid) in zip(i["centerFrequency"], i["eqbw"])])
        cerroPlateau = np.array([[averageTrans("CerroPlateau/", angle, 25, cent, wid), averageTrans("CerroPlateau/", angle, 50, cent, wid),
                                averageTrans("CerroPlateau/", angle, 75, cent, wid)] for (cent, wid) in zip(i["centerFrequency"], i["eqbw"])])
        cerroAPEX = np.array([[averageTrans("CerroAPEX/", angle, 25, cent, wid), averageTrans("CerroAPEX/", angle, 50, cent, wid),
                               averageTrans("CerroAPEX/", angle, 75, cent, wid)] for (cent, wid) in zip(i["centerFrequency"], i["eqbw"])])
        cerroConfig = np.array([[averageTrans("CerroConfig/", angle, 25, cent, wid), averageTrans("CerroConfig/", angle, 50, cent, wid),
                                 averageTrans("CerroConfig/", angle, 75, cent, wid)] for (cent, wid) in zip(i["centerFrequency"], i["eqbw"])])
        return {"ACT Site": quartileDisplay(
            actSiteTrans), "Steve Method": quartileDisplay(steveTrans), "Plateau": quartileDisplay(cerroPlateau), "APEX": quartileDisplay(cerroAPEX), "Config": quartileDisplay(cerroConfig)}

    temp = customOutput(45)
    temp.update({"Excel Sheet": quartileDisplay(i["eqtrans"])})
    temp.update({"Recreated Excel": quartileDisplay(np.array([[averageTransHelper("ACT_MAM_50_pwv0.51", cent, wid), averageTransHelper("ACT_MAM_50_pwv0.95", cent, wid),
                                                               averageTransHelper("ACT_MAM_50_pwv1.81", cent, wid)] for (cent, wid) in zip(i["centerFrequency"], i["eqbw"])]))})
    dict_file = {"45 Degree Elevation": temp,
                 "60 Degrees Elevation": customOutput(30)}
    return yaml.dump(dict_file, open("methods comparison.yaml", 'w'), sort_keys=False)


def sensitivityFile(outputs, valueDisplay, quartileDisplay):
    dict_file = {"NET w8 avg": valueDisplay(outputs["netW8Avg"]),
                 "NET w8 RJ": valueDisplay(outputs["netW8RJ"]), "NEI w8 Jy/sr": valueDisplay(
        outputs["neiW8"]), "EoR Spec NEFD": quartileDisplay(outputs["eorNEFD"]), "EoR Spec NEI": quartileDisplay(outputs["eorNEI"])}
    return yaml.dump(dict_file, open("output.yaml", 'w'), sort_keys=False)


def powerFile(outputs, quartileDisplay):
    outputs60 = calculate(30)
    outputs45 = calculate(45)

    def displayPower(o):
        return {"Broadband": quartileDisplay(o["powerPerPixel"]*10**12), "EoR Spec": quartileDisplay(o["eorPowerPerPixel"]*10**12)}

    def powerDiff(output1, output2):
        return {"powerPerPixel": np.array([[q1-q2 for q1, q2 in zip(f1, f2)] for f1, f2 in zip(output1["powerPerPixel"], output2["powerPerPixel"])]), "eorPowerPerPixel": np.array([[q1-q2 for q1, q2 in zip(f1, f2)] for f1, f2 in zip(output1["eorPowerPerPixel"], output2["eorPowerPerPixel"])])}

    dict_file = {str(i["angle"]) + " degrees": displayPower(outputs),
                 "45 degrees": displayPower(outputs45), "60 degrees": displayPower(outputs60), "45-60 Delta": displayPower(powerDiff(outputs45, outputs60))}
    return yaml.dump(dict_file, open("power.yaml", 'w'), sort_keys=False)


if __name__ == "__main__":
    i = getInputs("input.yaml")

    calculate = calcByAngle(i["diameter"], i["t"], i["wfe"], i["eta"], i["doe"], i["t_int"], i["pixelYield"], i["szCamNumPoln"], i["eorSpecNumPoln"],
                            i["t_filter_cold"], i["t_lens_cold"], i["t_uhdpe_window"], i["coldSpillOverEfficiency"], i["singleModedAOmegaLambda2"],
                            i["spatialPixels"], i["fpi"], i["eqbw"], i["centerFrequency"], i["detectorNEP"],
                            i["backgroundSubtractionDegradationFactor"], i["sensitivity"], i["hoursPerYear"], i["sensPerBeam"], i["r"], i["signal"])

    outputs = calculate(i["angle"])

    valueDisplay = valDisplayPartial(
        i["outputFreq"], i["centerFrequency"], outputs["wavelength"], i["decimalPlaces"])
    quartileDisplay = quartDisplayPartial(
        i["outputFreq"], i["centerFrequency"], outputs["wavelength"], i["decimalPlaces"])

    methodsComparisonFile(i, quartileDisplay)
    sensitivityFile(outputs, valueDisplay, quartileDisplay)
    powerFile(outputs, quartileDisplay)
