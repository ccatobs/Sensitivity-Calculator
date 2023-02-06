import mapsims
import healpy as hp
NSIDE = 16
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
