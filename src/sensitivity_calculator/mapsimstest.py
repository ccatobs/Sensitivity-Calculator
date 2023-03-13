import mapsims
import healpy as hp
import logging
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
#import so_pysm_models
#log = logging.getLogger("mapsims")
#logging.basicConfig(level=logging.INFO)
#log.setLevel(logging.INFO)
instrument_path = Path("/home/amm487/cloned_repos/Sensitivity-Calculator/src/sensitivity_calculator/data/instrument_parameters/instrument_parameters.tbl")
hitmap_path = Path("/home/amm487/cloned_repos/Sensitivity-Calculator/src/sensitivity_calculator/data")
NSIDE = 256
cmb = mapsims.SOPrecomputedCMB(
    num=0,
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
    SA_one_over_f_mode="pessimistic",
    instrument_parameters=instrument_path
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

chs = ["tube:LC1"]
final = []

for ch in chs:
    simulator = mapsims.MapSim(
        channels=ch,
        nside=NSIDE,
        unit="uK_CMB",
        pysm_output_reference_frame="G",
        pysm_components_string="a1",
        pysm_custom_components={"cmb": cmb},
        other_components={"noise": noise},
        instrument_parameters=instrument_path
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
