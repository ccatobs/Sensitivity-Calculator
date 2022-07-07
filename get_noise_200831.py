# import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import ccat_noise_200831 as CCAT_noise
import scipy.optimize as op
from astropy.table import QTable
#import mapsims

f = QTable.read("data/detectorParams.tbl", format="ascii.ipac")
print(f.colnames)
print(f)
print(f[0]["fwhm"])
f.add_index("band")
print(f.loc["70"]["center_frequency"])
print()
print()

if False:
    NSIDE = 16
    cmb = mapsims.SOPrecomputedCMB(
        num=0,
        nside=NSIDE,
        lensed=False,
        aberrated=False,
        has_polarization=True,
        cmb_set=0,
        cmb_dir="mapsims/tests/data",
        input_units="uK_CMB",
    )
    noise = mapsims.SONoiseSimulator(
        nside=NSIDE,
        return_uK_CMB=True,
        sensitivity_mode="baseline",
        apply_beam_correction=True,
        apply_kludge_correction=True,
        SA_one_over_f_mode="pessimistic",
    )
    simulator = mapsims.MapSim(
        channels="tube:ST0",
        nside=NSIDE,
        unit="uK_CMB",
        pysm_output_reference_frame="G",
        pysm_components_string="a1",
        pysm_custom_components={"cmb": cmb},
        other_components={"noise": noise},
    )
    output_map = simulator.execute()

plotAnything = False


CCAT_bands = np.array([222., 280., 348., 405., 850.])

# define sensitivity mode
# CCAT just has one mode so any input mode will work
mode = 'baseline'

# define fraction of sky
# fsky = 0.6
fsky = 15000./(4*np.pi*(180/np.pi)**2)
fsky = 20000./(4*np.pi*(180/np.pi)**2)

el = 45.
# sensitivity class
# SO nominal mission here assumes 5 years with efficiency of 0.2*0.85
# CCAT nominal mission assumes 4000 hours of observation
ccat = CCAT_noise.CcatLatv2b(mode,
                             survey_years=4000/24./365.24,
                             survey_efficiency=1.0,
                             N_tubes=(1, 1, 1, 1, 1), el=el)

# lmax for getting Nl's
lat_lmax = 10000  # Originally 8000

ell, N_ell_T_full, N_ell_P_full = ccat.get_noise_curves(
    fsky, lat_lmax, 1, full_covar=False, deconv_beam=True)

if plotAnything:
    plotTemperature = False
    for curve, label in zip(N_ell_T_full[:-1] if plotTemperature else N_ell_P_full[:-1], CCAT_bands[:-1]):
        plt.plot(ell, curve, label=str(int(label))+' GHz')
        plt.yscale('log')
        plt.ylim(10**-5, 10**3)
        plt.xscale('log')
        plt.xlim(10**2, 10**4)
        plt.title("Temperature" if plotTemperature else "Polarization")
    plt.legend(loc='upper right')
    plt.show()


WN_levels = ccat.get_white_noise(fsky)**.5

beam_sig_rad = ccat.get_beams() * np.pi/180/60 / (8.*np.log(2))**0.5
beams = np.exp(-0.5 * ell*(ell+1) * beam_sig_rad[:, None]**2)

bands = ccat.get_bands()

# N_ell_LA_T  = N_ell_LA_T_full[range(N_bands),range(N_bands)]
# N_ell_LA_Tx = [N_ell_LA_T_full[i,j] for i,j in corr_pairs]
# N_ell_LA_P  = N_ell_LA_P_full[range(N_bands),range(N_bands)]
# N_ell_LA_Px = [N_ell_LA_P_full[i,j] for i,j in corr_pairs]

print("band centers: ", ccat.get_bands()[:], "[GHz]")
print("beam sizes: ", ccat.get_beams()[:]*60, "[arcsec]")

print("white noise levels: ", WN_levels[:], "[uK-arcmin]")


def T_tot_nl(tempPython3Conversion, Nwhite, Nred):
    ell, bl = tempPython3Conversion
    ell0 = 1000.
    alpha = -3.5
    return (Nred*(ell/ell0)**(alpha) + Nwhite)/bl**2


def P_tot_nl(tempPython3Conversion, Nwhite):
    ell, bl = tempPython3Conversion
    ell0 = 700.
    alpha = -1.4
    return (Nwhite*(ell/ell0)**(alpha) + Nwhite)/bl**2


for i in range(len(N_ell_T_full)):
    guess = [0.1, 1.0]
    popt = op.curve_fit(T_tot_nl, (ell, beams[i]),
                        N_ell_T_full[i], p0=guess)[0]
    print(bands[i], 'N_white=%E' % popt[0], 'N_red=%E' % popt[1])

    # guess = 0.1
    # popt = op.curve_fit(P_tot_nl,(ell,beams[i]),
    #                     N_ell_P_full[i],p0=guess)[0]
    # print(bands[i],'N_white=%E'%popt[0])


# print N_ell_T_full.shape
# print N_ell_P_full.shape
# print N_ell_T_full[:,1]
# # plotting
# plt.subplot(131)

# for i,f in enumerate(SO_bands):
# 	plt.plot(ell_LA,N_ell_LA_T_full[i],linewidth=2,label='S %i'%(int(f)))

# for i,f in enumerate(CCAT_bands):
# 	plt.plot(ell,N_ell_T_full[i],linewidth=2,label='C %i'%(int(f)),linestyle='--')

# plt.title("LAT T")
# plt.legend(loc=0)
# plt.yscale('log')
# plt.xlabel("$\ell$",fontsize=18)
# plt.ylabel("$N_\ell \, [\mu\mathrm{K}^2]$",fontsize=18)

# plt.subplot(132)

# for i,f in enumerate(SO_bands):
# 	plt.plot(ell_LA,N_ell_LA_P_full[i],linewidth=2,label='S %i'%(int(f)))

# for i,f in enumerate(CCAT_bands):
# 	plt.plot(ell,N_ell_P_full[i],linewidth=2,label='C %i'%(int(f)),linestyle='--')

# plt.title("LAT P")
# plt.legend(loc=0)
# plt.yscale('log')
# plt.ylim(1e-6,1e9)
# plt.xlabel("$\ell$",fontsize=18)

# plt.subplot(133)

# for i,f in enumerate(SO_bands):
# 	plt.plot(ell_SA,N_ell_SA_P_full[i],label='S %i'%(int(f)),linewidth=2)

# plt.title("SAT P")
# plt.legend(loc=0)
# plt.yscale('log')
# plt.xlabel("$\ell$",fontsize=18)
# fig = mpl.pyplot.gcf()
# fig.set_size_inches(12.0,6.0)
# plt.savefig('noise_200414.png')
# plt.clf()
