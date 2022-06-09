import yaml
import numpy as np

#Declare all variables with line breaks repesenting section breaks in the excel file this code is based on
diameter = None
t = None
wfe = None
space = None

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

aToCMB = None
wavelength = None

#These functions work on numpy arrays (componentwise) and help with code clarity
pi = np.pi
ln = np.log
sqrt = np.sqrt
e = np.e
h = 6.626e-34
k = 1.38e-23
c = 299792458

#Read in yaml file and assign all input data to their respective variables
stream = open("input.yaml", 'r')
dictionary = yaml.safe_load(stream)
for key, value in dictionary.items():
    if key == "diameter": diameter = value
    if key == "t": t = value
    if key == "wfe": wfe = value
    if key == "space": space = value

    if key == "eta": eta = value

    if key == "doe": doe = value
    if key == "t_int": t_int = value
    if key == "pixelYield": pixelYield = value
    if key == "sz-camNumPoln": szCamNumPoln = value
    if key == "eor-specNumPoln": eorSpecNumPoln = value

    if key == "t_filter_cold": t_filter_cold = np.array(value)
    if key == "t_lens_cold": t_lens_cold = np.array(value)
    if key == "t_uhdpe_window": t_uhdpe_window = np.array(value)
    if key == "coldSpillOverEfficiency": coldSpillOverEfficiency = np.array(value)
    if key == "singleModedAOmegaLambda2": singleModedAOmegaLambda2 = np.array(value)
    if key == "spatialPixels": spatialPixels = np.array(value)
    if key == "fpi": fpi = np.array(value)

    if key == "eqbw": eqbw = np.array(value)
    if key == "eqtrans": eqtrans = np.array(value)
    if key == "centerFrequency": centerFrequency = np.array(value)
    if key == "detectorNEP": detectorNEP = value
    if key == "backgroundSubtractionDegradationFactor": backgroundSubtractionDegradationFactor = value
    if key == "sensitivity": sensitivity = value
    if key == "hoursPerYear": hoursPerYear = value
    if key == "sensPerBeam": sensPerBeam = value

    if key == "r": r = np.array(value)
    if key == "signal": signal = np.array(value)

    if key == "aToCMB": aToCMB = np.array(value)
    if key == "wavelength": wavelength = np.array(value)

a = pi*(diameter/2)**2

szCamCorForPoln = abs(3-szCamNumPoln)
eorSpecCorForPoln = abs(3-eorSpecNumPoln)

t_cold = t_filter_cold*t_lens_cold**3
e_window_warm = 1 - t_uhdpe_window

windowTrans = np.ones(5) * t_uhdpe_window[1]
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




e_warm = (t_uhdpe_window)*((1-eqtrans)*eta+(1-eta))+(e_window_warm)
for i in range(5):
    e_warm[i][-1] = 1
ruze = 1/(e**((4*pi*wfe/(wavelength))**2))
weather_t_cold = (t_cold)
occupancy_normal = np.ones((5, 4)) * 1/(e**(h*c/((wavelength)*10**(-6)*k*t))-1)[:, None]
occupancy_space = np.ones((5, 1)) * 1/(e**(h*c/((wavelength)*10**(-6)*k*space))-1)[:, None]
occupancy = np.concatenate((occupancy_normal, occupancy_space), axis=1)
acceptedModes = szCamNumPoln*(coldSpillOverEfficiency)*(singleModedAOmegaLambda2)
powerPerPixel = h*c*e_warm*(weather_t_cold[:, None])*(acceptedModes[:, None])*occupancy*(eqbw[:, None])*doe/(((wavelength)*10**(-6))[:, None]) #Differs from excel sheet becauses of correctly using t(cold)
photonNoiseNEP = h*c/((wavelength)[:, None]*10**(-6))*(1/t_int*acceptedModes[:, None]*eqbw[:, None]*e_warm*weather_t_cold[:, None]*doe*occupancy*(1+e_warm*weather_t_cold[:, None]*doe*occupancy))**0.5
totalNEP = sqrt(photonNoiseNEP**2 + detectorNEP**2)
coldTerminatedSpillover = (coldSpillOverEfficiency)
nef = szCamCorForPoln*totalNEP/(a*weather_t_cold[:, None]*doe*eta*ruze[:, None]*eqtrans*coldTerminatedSpillover[:, None]*backgroundSubtractionDegradationFactor)
nefd = nef*10**26*1000/(eqbw[:, None])
net = nefd*10**(-29)/((solidAngle)[:, None])*((wavelength)[:, None]*0.000001)**2/(2*1.38*10**(-23))*1000

arrayNETRJ = (net/((spatialPixels)[:, None]*pixelYield)**0.5)[:, :3]
arrayNETCMB = arrayNETRJ*aToCMB[:, None]*1000
arrayNEI = nefd[:, :3]/1000/((solidAngle)[:, None])/((spatialPixels)[:, None]*pixelYield)**0.5
netW8Avg = sqrt(3/(1/arrayNETCMB[:, 0]**2+1/arrayNETCMB[:, 1]**2+1/arrayNETCMB[:, 2]**2).astype(float))
netW8RJ = netW8Avg/aToCMB
neiW8 = sqrt(3/(1/arrayNEI[:, 0]**2+1/arrayNEI[:, 1]**2+1/arrayNEI[:, 2]**2).astype(float))
nefd = arrayNEI*(solidAngle)[:, None]*sqrt((spatialPixels)[:, None]*pixelYield)*1000

def eorRed(array): #Eliminates 353 and 448 um wavelengths from data
    return array[1:]

r = np.concatenate((r[:1, :4], r[2:, :4]))
eorEqBw = c/(eorRed(wavelength)[:, None]*10**(-6)*r)*pi/2
eorEqTrans = eqtrans[1:, :4]
eorE_warm = eorRed(t_uhdpe_window)[:, None]*((1-eorEqTrans)*eta+(1-eta))+eorRed(e_window_warm)[:, None]
centerFrequency = 299000000000000/eorRed(wavelength)
eorRuze = 1/e**((4*pi*wfe/eorRed(wavelength))**2)
eorT_cold = eorRed(t_cold)*0.9 #Always takes from 739 um in the sheet
eorOccupancy = 1/(e**(h*c/(eorRed(wavelength)*10**(-6)*k*t))-1)
eorAcceptedModes = eorSpecNumPoln*a*eorRed(coldSpillOverEfficiency)*eorRed(solidAngle)/(eorRed(wavelength)*10**(-6))**2 #In sheet 1071 um is calculated incorrectly
eorPhotonNoiseNEP = h*c/(eorRed(wavelength)[:, None]*10**(-6))*(eorAcceptedModes[:, None]*eorEqBw*eorE_warm*eorT_cold[:, None]*doe*eorOccupancy[:, None]*(1+eorE_warm*eorT_cold[:, None]*doe*eorOccupancy[:, None]))**0.5
eorTotalNEP = sqrt(eorPhotonNoiseNEP**2 + detectorNEP**2)
eorColdTerminatedSpillover = eorRed(coldSpillOverEfficiency)
eorNEF = eorSpecCorForPoln*eorTotalNEP/(a*eorT_cold[:, None]*doe*eta*eorRuze[:, None]*eorEqTrans*eorColdTerminatedSpillover[:, None]*backgroundSubtractionDegradationFactor)
eorNEFD = eorNEF/eorEqBw*10**26*1000
eorNEI = eorNEFD/eorRed(solidAngle)[:, None]/1000

arrayify = np.ndarray.tolist

def broadbandDisplayHelper(array, w, dict):
    if len(w) > 0:
        dict.update({str(w[0]) + " um": array[0]})
        return broadbandDisplayHelper(array[1:], w[1:], dict)
    else:
        return dict


def broadbandDisplay(array):
    array = arrayify(array)
    return broadbandDisplayHelper(array, wavelength, {})


def eoRDisplayHelper(array, w, dict):
    if len(w) > 0:
        dict.update({str(w[0]) + " um": {"Quartile 1": array[0][0], "Quartile 2": array[0][1], "Quartile 3": array[0][2], "Quartile 4": array[0][3]}})
        return eoRDisplayHelper(array[1:], w[1:], dict)
    else:
        return dict

def eoRDisplay(array):
    array = arrayify(array)
    return eoRDisplayHelper(array, eorRed(wavelength), {})

dict_file = {"NET w8 avg": broadbandDisplay(netW8Avg), "NET w8 RJ": broadbandDisplay(netW8RJ), "NEI w8 Jy/sr": broadbandDisplay(neiW8), "EoR Spec NEFD": eoRDisplay(eorNEFD), "EoR Spec NEI": eoRDisplay(eorNEI)}
documents = yaml.dump(dict_file, open("output.yaml", 'w'), sort_keys=False)