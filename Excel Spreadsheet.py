import yaml
import numpy as np

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

pi = np.pi
ln = np.log
sqrt = np.sqrt
e = np.e
h = 6.626e-34
k = 1.38e-23
c = 299792458
wavelength = np.array([353, 448, 740, 861, 1071, 1350])

stream = open("input.yaml", 'r')
dictionary = yaml.safe_load(stream)
for key, value in dictionary.items():
    match key:
        case "diameter": diameter = value
        case "t": t = value
        case "wfe": wfe = value
        case "space": space = value

        case "eta": eta = value

        case "doe": doe = value
        case "t_int": t_int = value
        case "pixelYield": pixelYield = value
        case "sz-camNumPoln": szCamNumPoln = value
        case "eor-specNumPoln": eorSpecNumPoln = value

        case "t_filter_cold": t_filter_cold = np.array(value)
        case "t_lens_cold": t_lens_cold = np.array(value)
        case "t_uhdpe_window": t_uhdpe_window = np.array(value)
        case "coldSpillOverEfficiency": coldSpillOverEfficiency = np.array(value)
        case "singleModedAOmegaLambda2": singleModedAOmegaLambda2 = np.array(value)
        case "spatialPixels": spatialPixels = np.array(value[:1] + [1] + value[1:]) #Dummy value in unused position [1] of list
        case "fpi": fpi = np.array(value)

        case "eqbw": eqbw = np.array(value)
        case "eqtrans": eqtrans = np.array(value)
        case "centerFrequency": centerFrequency = np.array(value)
        case "detectorNEP": detectorNEP = value
        case "backgroundSubtractionDegradationFactor": backgroundSubtractionDegradationFactor = value
        case "sensitivity": sensitivity = value
        case "hoursPerYear": hoursPerYear = value
        case "sensPerBeam": sensPerBeam = value

        case "r": r = np.array(value)
        case "signal": signal = np.array(value)

        case "aToCMB": aToCMB = np.array(value)

a = pi*(diameter/2)**2

szCamCorForPoln = abs(3-szCamNumPoln)
eorSpecCorForPoln = abs(3-eorSpecNumPoln)

t_cold = t_filter_cold*t_lens_cold**3
e_window_warm = 1 - t_uhdpe_window

windowTrans = np.ones(6) * t_uhdpe_window[1]
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


def reduce(array): #Eliminates 448 um wavelength from data
    return np.concatenate((np.array([array[0]]), array[2:]))


e_warm = reduce(t_uhdpe_window)*((1-eqtrans)*eta+(1-eta))+reduce(e_window_warm)
for i in range(5):
    e_warm[i][-1] = 1
ruze = 1/(e**((4*pi*wfe/reduce(wavelength))**2))
weather_t_cold = reduce(t_cold)
occupancy_normal = np.ones((5, 4)) * 1/(e**(h*c/(reduce(wavelength)*10**(-6)*k*t))-1)[:, None]
occupancy_space = np.ones((5, 1)) * 1/(e**(h*c/(reduce(wavelength)*10**(-6)*k*space))-1)[:, None]
occupancy = np.concatenate((occupancy_normal, occupancy_space), axis=1)
acceptedModes = szCamNumPoln*reduce(coldSpillOverEfficiency)*reduce(singleModedAOmegaLambda2)
powerPerPixel = h*c*e_warm*(weather_t_cold[:, None])*(acceptedModes[:, None])*occupancy*(eqbw[:, None])*doe/((reduce(wavelength)*10**(-6))[:, None]) #Differs from excel sheet becauses of correctly using t(cold)
photonNoiseNEP = h*c/(reduce(wavelength)[:, None]*10**(-6))*(1/t_int*acceptedModes[:, None]*eqbw[:, None]*e_warm*weather_t_cold[:, None]*doe*occupancy*(1+e_warm*weather_t_cold[:, None]*doe*occupancy))**0.5
totalNEP = sqrt(photonNoiseNEP**2 + detectorNEP**2)
coldTerminatedSpillover = reduce(coldSpillOverEfficiency)
nef = szCamCorForPoln*totalNEP/(a*weather_t_cold[:, None]*doe*eta*ruze[:, None]*eqtrans*coldTerminatedSpillover[:, None]*backgroundSubtractionDegradationFactor)
nefd = nef*10**26*1000/(eqbw[:, None])
net = nefd*10**(-29)/(reduce(solidAngle)[:, None])*(reduce(wavelength)[:, None]*0.000001)**2/(2*1.38*10**(-23))*1000

arrayNETRJ = (net/(reduce(spatialPixels)[:, None]*pixelYield)**0.5)[:, :3]
arrayNETCMB = arrayNETRJ*aToCMB[:, None]*1000
arrayNEI = nefd[:, :3]/1000/(reduce(solidAngle)[:, None])/(reduce(spatialPixels)[:, None]*pixelYield)**0.5
netW8Avg = sqrt(3/(1/arrayNETCMB[:, 0]**2+1/arrayNETCMB[:, 1]**2+1/arrayNETCMB[:, 2]**2).astype(float))
netW8RJ = netW8Avg/aToCMB
neiW8 = sqrt(3/(1/arrayNEI[:, 0]**2+1/arrayNEI[:, 1]**2+1/arrayNEI[:, 2]**2).astype(float))
nefd = arrayNEI*reduce(solidAngle)[:, None]*sqrt(reduce(spatialPixels)[:, None]*pixelYield)*1000

def eorRed(array): #Eliminates 353 and 448 um wavelengths from data
    return array[2:]

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

def broadbandDisplay(array):
    array = arrayify(array)
    return {"353 um": array[0], "740 um": array[1], "861 um": array[2], "1052 um": array[3], "1350 um": array[4]}

def eoRDisplay(array):
    array = arrayify(array)
    return {"740 um": {"Quartile 1": array[0][0], "Quartile 2": array[0][1], "Quartile 3": array[0][2], "Quartile 4": array[0][3]}, 
    "861 um": {"Quartile 1": array[1][0], "Quartile 2": array[1][1], "Quartile 3": array[1][2], "Quartile 4": array[1][3]}, 
    "1071 um": {"Quartile 1": array[2][0], "Quartile 2": array[2][1], "Quartile 3": array[2][2], "Quartile 4": array[2][3]}, 
    "1350 um": {"Quartile 1": array[3][0], "Quartile 2": array[3][1], "Quartile 3": array[3][2], "Quartile 4": array[3][3]}}

dict_file = {"NET w8 avg": broadbandDisplay(netW8Avg), "NET w8 RJ": broadbandDisplay(netW8RJ), "NEI w8 Jy/sr": broadbandDisplay(neiW8), "EoR Spec NEFD": eoRDisplay(eorNEFD), "EoR Spec NEI": eoRDisplay(eorNEI)}
documents = yaml.dump(dict_file, open("output.yaml", 'w'))