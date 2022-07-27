# TODO: data_C, ccat_pwv_cor
import pwvCalculator as pwv
import numpy as np


def get_atmosphere_C(freqs, el=None, data_C=None):
    """
    Returns atmospheric noise power at ell=1000, for an ACTPol optics
    tube.  In units of [uK^2 sec].  This only works for a few special
    frequencies.

    Basic model assumes el=50.  A simple rescaling (proportional to
    csc(el)) is applied for other values of el.

    version=0: This atmospheric model was used in SO V3 forecasts but
    contains an error.

    version=1: This atmospheric model is better than version=0, in
    that at least one error has been corrected.  The relative
    atmospheric powers have been measured in AT model, and then
    calibrated to ACT.  Low frequency results are inflated somewhat
    because ACT sees more power at 90 GHz than predicted by this
    modeling.
    """
    # This el_correction is quite naive; atmosphere may be less
    # structured (i.e. smoothed out) at lower elevation because of the
    # increased thickness.
    if el is None:
        el = 50.
    el_correction = np.sin(50.*np.pi/180) / np.sin(el*np.pi/180)
    # Ratio between CCAT site and ACT site
    ccat_pwv_cor = 0.67/pwv.configPWV(50)
    data_bands = np.array(freqs)
    if data_C is None:
        data_C = np.array([
            # below factors from am_output/mult_pwv/get_derivative_ccat.py
            2.31956542e+05,
            1.61527385e+06,
            4.03473727e+07,
            2.51490116e+08,
            9.10884821e+13
        ])
    data = {}
    for b, C in zip(data_bands, data_C):
        data[b] = C
    return np.array([data[f] * (ccat_pwv_cor*el_correction)**2 for f in freqs])


class SOLatType:
    def __init__(self, *args, **kwargs):
        raise RuntimeError('You should subclass this.')

    def get_bands(self):
        return self.bands.copy()

    def get_beams(self):
        return self.beams.copy()

    def precompute(self, N_tubes, N_tels=1, data_C=None):
        # Accumulate total white noise level and atmospheric
        # covariance matrix for this configuration.

        band_weights = np.zeros(self.n_bands)
        for (tube_name, tube_count) in N_tubes:
            # commenting out white noise el rescaling as the white noise computed already has this in there.
            # * white_noise_el_rescale
            tube_noise = self.tube_configs[tube_name]
            s = (tube_noise != 0)
            band_weights[s] += tube_count * N_tels * tube_noise[s]**-2

        self.band_sens = np.zeros(self.n_bands) + 1e9
        s = (band_weights > 0)
        self.band_sens[s] = band_weights[s]**-0.5

        # Special for atmospheric noise model.
        self.Tatmos_C = get_atmosphere_C(
            self.bands, el=self.el, data_C=data_C) * self.FOV_mod
        self.Tatmos_ell = 1000. + np.zeros(self.n_bands)
        self.Tatmos_alpha = -3.5 + np.zeros(self.n_bands)

        # Compute covariant weight matrix (atmosphere parameters).
        cov_weight = np.zeros((self.n_bands, self.n_bands))
        pcov_weight = np.zeros((self.n_bands, self.n_bands))
        atm_rho = 0.9
        for (tube_name, tube_count) in N_tubes:
            # Get the list of coupled bands; e.g. [1,2] for MF.
            nonz = self.tube_configs[tube_name].nonzero()[0]
            for i in nonz:
                for j in nonz:
                    w = {True: 1., False: atm_rho}[i == j]
                    assert(cov_weight[i, j] == 0.)  # Can't do overlapping
                    # tubes without weights.
                    cov_weight[i, j] += tube_count * N_tels / (w * (
                        self.Tatmos_C[i] * self.Tatmos_C[j])**.5)
                    pcov_weight[i, j] = w

        # Reciprocate non-zero elements.
        s = (cov_weight != 0)
        self.Tatmos_cov = np.diag([1e9]*self.n_bands)
        self.Tatmos_cov[s] = 1./cov_weight[s]

        # Polarization is simpler...
        self.Patmos_ell = 700. + np.zeros(self.n_bands)
        self.Patmos_alpha = -1.4 + np.zeros(self.n_bands)

        self.Patmos_cov = pcov_weight

    def get_survey_time(self):
        t = self.survey_years * 365.25 * 86400.  # convert years to seconds
        return t * self.survey_efficiency

    def get_survey_spread(self, f_sky, units='arcmin2'):
        # Factor that converts uK^2 sec -> uK^2 arcmin^2.
        A = f_sky * 4*np.pi
        if units == 'arcmin2':
            A *= (60*180/np.pi)**2
        elif units != 'sr':
            raise ValueError("Unknown units '%s'." % units)
        return A / self.get_survey_time()

    def get_white_noise(self, f_sky, units='arcmin2'):
        return self.band_sens**2 * self.get_survey_spread(f_sky, units=units)

    def get_noise_curves(self, f_sky, ell_max, delta_ell, deconv_beam=True,
                         full_covar=False):
        ell = np.arange(2, ell_max, delta_ell)
        W = self.band_sens**2

        # Get full covariance; indices are [band,band,ell]
        ellf = (ell/self.Tatmos_ell[:, None])**(self.Tatmos_alpha[:, None])
        T_noise = self.Tatmos_cov[:, :, None] * \
            (ellf[:, None, :] * ellf[None, :, :])**.5

        # P noise is tied directly to the white noise level.
        P_low_noise = (2*W[:, None]) * (ell /
                                        self.Patmos_ell[:, None])**self.Patmos_alpha[:, None]
        P_noise = (self.Patmos_cov[:, :, None] *
                   (P_low_noise[:, None, :] * P_low_noise[None, :, :])**.5)

        # Add in white noise on the diagonal.
        for i in range(len(W)):
            T_noise[i, i] += W[i]
            P_noise[i, i] += W[i] * 2

        # Deconvolve beams.
        if deconv_beam:
            beam_sig_rad = self.get_beams() * np.pi/180/60 / (8.*np.log(2))**0.5
            beams = np.exp(-0.5 * ell*(ell+1) * beam_sig_rad[:, None]**2)
            T_noise /= (beams[:, None, :] * beams[None, :, :])
            P_noise /= (beams[:, None, :] * beams[None, :, :])

        # Diagonal only?
        if not full_covar:
            ii = range(self.n_bands)
            T_noise = T_noise[ii, ii]
            P_noise = P_noise[ii, ii]

        sr_per_arcmin2 = (np.pi/180/60)**2
        return (ell,
                T_noise * self.get_survey_spread(f_sky, units='sr'),
                P_noise * self.get_survey_spread(f_sky, units='sr'))


def el_noise_func(P, el):
    a, b = P
    return a + b / np.sin(el*np.pi/180)


class CCAT(SOLatType):
    """This special edition S4 LAT is equipped with an additional tube
    class, XHF, containing 280 and 350 GHz detectors and blessed with
    the calm atmosphere available at the CCAT site.  Otherwise it
    should reproduce S4Lat results exactly.  This is a candidate for
    introduction to the main noise model.

    """
    atm_version = 1

    def __init__(self, centerFrequency, beam, net, data_C=None,
                 N_tubes=None, N_tels=None,
                 survey_years=4000/24./365.24,
                 survey_efficiency=1.0,
                 el=None):
        # Define the instrument.
        self.bands = np.array(centerFrequency)
        # scaled beam for 410 GHz but need to check 350 and 410.
        self.beams = np.array(beam)
        self.n_bands = len(self.bands)

        # Set defaults for survey area, time, efficiency
        self.survey_years = survey_years
        self.survey_efficiency = survey_efficiency

        # Sensitivities of each kind of optics tube, in uK rtsec, by
        # band.  0 represents 0 weight, not 0 noise...
        nar = np.array
        self.tube_configs = {
            # 'ULF': nar([ 60.,   0,   0,   0,   0,   0,   0,    0,    0, 0]),
            # 'LF':  nar([   0,28.6,16.0,   0,   0,   0,   0,    0,    0, 0]),
            # 'MF':  nar([   0,   0,   0, 6.6, 6.8,   0,   0,    0,    0, 0]),
            # from MFH:
            # 'UHF': nar([   0,   0,   0,   0,   0,12.5,30.0,    0,    0,0]),
            # from gordon's spreadsheet for 1
            'HF1': nar([net[0],    0,    0,    0, 0]),
            'HF2': nar([0, net[1],    0,    0, 0]),
            'HF3': nar([0,   0, net[2],    0, 0]),
            'HF4': nar([0,   0,    0, net[3], 0]),
            'HF5': nar([0,   0,    0,    0, net[4]]),
        }

        # Save the elevation request.
        self.el = el

        # Factor by which to attenuate atmospheric power, given FOV
        # relative to ACT?
        # print("check fov!")
        self.FOV_mod = 0.5

        # The reference tube config.
        ref_tubes = [('HF1', 1), ('HF2', 1), ('HF3', 1),
                     ('HF4', 1), ('HF5', 1)]

        if N_tels is None:
            N_tels = 1

        if N_tubes is None:
            N_tubes = ref_tubes
        else:
            N_tubes = [(b, x) for (b, n), x in zip(ref_tubes, N_tubes)]

        self.precompute(N_tubes, N_tels, data_C=data_C)
