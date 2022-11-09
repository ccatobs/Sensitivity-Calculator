import mapsims
import os
import sys
os.environ['KMP_DUPLICATE_LIB_OK']='True'
import numpy as np
import healpy as hp
from matplotlib import pyplot as plt
from pixell import enmap, enplot, reproject, utils, curvedsky
import pysm3.units as u
import pysm3 as pysm

def get_window(mapk,n_it=5,width = 0.1):
    """Helper to apodize_map"""
    m_apo = np.copy(mapk)*0
    if hp.UNSEEN in np.copy(mapk):
        m_apo[np.copy(mapk)!=hp.UNSEEN] = 1.0
    else:
        m_apo[np.copy(mapk)!=0] = 1.0
    
    return smooth_map(m_apo, n_it=5, width = 0.1)

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
        
    

def rough_change_nside(m, nside):
    """takes any map m and roughly translates it to a corresponding map in a different nside"""
    nr = m.shape[-1]
    gr = hp.nside2npix(nside)
    nested = hp.reorder(m,inp = 'RING', out = 'NESTED')
    if nr>gr:
        step = int(nr/gr)
        B = (np.arange(gr)*step).astype(int)
        return hp.reorder(nested[B],inp='NESTED', out = 'RING')
    elif nr<gr:
        assert (gr/nr)%4==0, "bad npix side"
        mul = np.ones(int(gr/nr))
        rres = np.outer(nested,mul).flatten()
        return hp.smoothing(hp.reorder(rres,inp='NESTED', out = 'RING'),fwhm = hp.nside2resol(hp.npix2nside(nr)))
        
    else:
        return m
    

def apodize_map(map0,n_it =5):
    """Apodizes a map with hp.UNSEEN pixels as masks. n_it represents the iterations of gaussian convolutions. Higher n_it generally means bigger windows and lower n_it the inverse."""
    
    tmp = np.copy(map0)
    tmp2 = np.copy(map0)
    tmp1 = tmp!=hp.UNSEEN
    m_apo = get_window(tmp2, n_it = n_it)
    tmp[tmp1!=True] =0 
    output = tmp*m_apo
    return output


def lmax_f(NSIDE):
    """Returns $\ell$ max given NSIDE"""
    return 3*NSIDE-1

def get_w(map):
    """Returns the skyfraction of map, w"""
    m_apo = get_window(map)
    sh = m_apo.shape
    mt = m_apo.reshape(-1)
    npix = len(mt)
    non0 = np.where(mt!=0)
    w2 = np.sum(mt[non0]**2)/(npix)
    w4 = np.sum(mt[non0]**4)/(npix)
    mt = mt.reshape(sh)
    w = w2**2/w4
    return w
    
# def transform_cmb_map(cmb, noise, ch):
#     "ch is of the form 'LT6','LT5',etc... "
#     chs = noise.tubes[ch]
#     sky = pysm.Sky(nside=noise.nside, output_unit = cmb.input_units)
#     sky.components.append(cmb)
#     output_map = []
#     for each in chs:
#         print(each)
#         bandpass_integrated_map = sky.get_emission(
#             *each.bandpass
#         ).value
#         beam_width_arcmin = each.beam
#         # smoothing and coordinate rotation with 1 spherical harmonics transform
#         smoothed_map = hp.ma(
#             pysm.apply_smoothing_and_coord_transform(
#                 bandpass_integrated_map,
#                 fwhm=beam_width_arcmin,
#                 lmax=3 * noise.nside - 1,
#                 rot=None,
#                 map_dist=None
#             )
#         )
#         output_map.append(smoothed_map)
#     output_map = np.array(output_map)
#     return output_map


def write_output_map(output,output_file,filename):
    """Writes a map as a .fits file. output is the map, or set of maps. output_file is the directory in which the file will be written to. filename is the name given to the file. Make sure """
    for det in output.keys():
        hp.fitsfunc.write_map(output_file + "/" + str(det)+"_"+filename, output[det], overwrite = True)
    return None

def read_output_map(output_file, filename, noise, ch):
    num = len(noise.tubes[ch[5:]])
    res = np.zeros((num,3,12*(noise.nside)**2))
    for i in np.arange(num):
        tmp = noise.tubes[ch[5:]][i]
        det = tmp.tag
        data = hp.fitsfunc.read_map(output_file + "/" + det + "_" + filename)
        res[i] = data

    return np.array(res)

def cls_to_dls(cls):
    """Coverts a set of Cls into Dls by their canonical relationship"""
    ell = np.arange(np.array(cls).shape[-1])
    Dls = np.zeros(np.array(cls).shape)
    if len(cls.shape)>1:
        for i in np.arange(len(cls)):
            Dls[i] = cls[i]*ell*(ell+1) /2/np.pi
    else:
        Dls = cls*ell*(ell+1) /2/np.pi
    return Dls,ell


# def deconvolve_beam(cls, beam_width):
#     """beam_width must be in degrees
#        spec is the power spectrum to be deconvolved"""
#     ell = np.arange(cls.shape[-1])
#     beam_sig_rad = beam_width*np.pi/180/60 / (8.*np.log(2))**0.5
#     beam = np.exp(-0.5 * ell*(ell+1) * beam_sig_rad**2)
#     result = np.array(cls)/(beam[None,:]*beam[None,:])
#     return result[0]

def count_ccat_detectors():
    """Counts the number of files beginning in mapsims/data/simonsobs_instrument_parameters_2020.06 beginning with bandpass_HF. 
    So, make sure any detector you wish to simulate has a bandpass .tbl file in this directory"""
    unfil = os.listdir('mapsims/data/simonsobs_instrument_parameters_2020.06')
    counter=0
    for file in unfil:
        if file[0:11]=="bandpass_HF":
            counter+=1
    return counter

def get_bands(chs, noise, remove =None):
    bands = []
    for ch in chs:
        bands.append(noise.tubes[ch[5:]][0])
        bands.append(noise.tubes[ch[5:]][1])
    inds = np.ones(len(bands))
    if remove!=None:
        inds[np.array(remove)]=0
    
    return np.array(bands)[(inds.astype(int)).astype(bool)]


def read_sims(sim_num, noise, NSIDE, output_file, pysm_string, chs, apodize=False, remove=None):
    """Output shape will be of the form (5, sim_num , 3, n_pix). Reads files of the form LCX_HFX_LCX_NSIDE_XXX_TAG_pysm_string_sim_index_.fits. This function is best used to load up files from several simulations with the same configurations.
    sim_num: number of simulations you would like to load
    noise: the mapsims noise object
    NSIDE: the nside of the stored simulations
    output_file: the directroy from which to retrieve the files
    pysm_string: the pysm_string used to feed into the simulations. e.g. 'd0', 's0', 's0,d0' 
    chs: the channels you would like retrieve files of. e.g. ["tube:LC1"], ["tube:LC2","tube:LC3"], ...
    apodize: whether you would like to apodize the maps as they are retrieved
    remove: an index or a list thereof of bands you wish to remove. e.g. [5] to remove HF6 when chs = ["tube:LC1","tube:LC2","tube:LC3"]"""
    
    bands = get_bands(chs,noise, remove = remove)
    output = np.zeros((len(bands),sim_num, 3, 12*(NSIDE)**2))
    for j in np.arange(sim_num):
        for k in np.arange(len(bands)):
            tmp = hp.fitsfunc.read_map(output_file +"/"+bands[k].tag+"_" + bands[k].tag[:3]+"_NSIDE_" + str(NSIDE) + "_TAG_" + pysm_string + "_" + str(j) + "_" +".fits", field=(0,1,2))
            if apodize:
                for h in np.arange(tmp.shape[0]):
                    tmp[h,:] = apodize_map(tmp[h,:])
                    
            output[k,j,:,:] = tmp
    return output

def compute_cross_spectra(sim_data, lat_lmax=None, p = None,auto_only=True):
    """Computes cross frequency spectra for a set of maps from different simulations. 
    sim_data should be of the form (#_of_freqs,sim_num, #_of_pols, n_pix )
    Returns cross spectra of the form (#_of_freqs, #_of_freqs, sim_num, lmax).
    p: polarization to return (0=TT, 1=EE, 2=BB). Leave None to return all polarizations.
    auto_only: whether to return only auto spectra or cross spectra"""
    data = sim_data
    if p ==None:
        result = np.zeros((sim_data.shape[0],sim_data.shape[0], sim_data.shape[1], 6, lat_lmax+1))
    else:
        result = np.zeros((sim_data.shape[0],sim_data.shape[0], sim_data.shape[1], lat_lmax+1))
        data = sim_data[:,:,p,:]
    if auto_only:
        for k in np.arange(sim_data.shape[1]):
            for i in np.arange(sim_data.shape[0]):
                result[i,i,k] = hp.sphtfunc.anafast(data[i,k], lmax = lat_lmax)
    else:
        for k in np.arange(sim_data.shape[1]):
            for i in np.arange(sim_data.shape[0]):
                for j in np.arange(sim_data.shape[0]):
                    if i<=j:
                        result[i,j,k] = hp.sphtfunc.anafast(data[i,k],data[j,k], lmax = lat_lmax)
                    else:
                        result[i,j,k] = result[j,i,k]
    return result

    
def bin_array(array, binl):
    """Bins array on the last axis by taking an average of binl elements. e.g. (n1,n2,n3,...,nN) -> (n1,n2,n3,...,nN/binl).
    array: multi-axis array to be binned
    binl: bin size"""
    shape = array.shape
    assert shape[-1]/binl == int(shape[-1]/binl),"The bin length does not fit properly in the length of shape[-1]"
    n = int(shape[-1]/binl)
    T = np.outer(np.ones(n),np.concatenate((np.ones(binl), np.zeros(binl*n-binl))))
    for i in np.arange(T.shape[0]):
        T[i] = np.roll(T[i],i*binl)
    T /=binl
    return np.dot(array,T.T)

def compute_PS_matrix(ps_list, diag=True):
    """Computes a naive frequency cross spectrum matrix using Cl_ij = sqrt(Cl_i * Cl*j).
    ps_list: Should be of the shape (#_of_freqs,ell_max)
    diag: whether to make the off-diagonals 0 or not
    Returns a shape (#_of_freqs,#_of_freqs,ell_max)"""
    
    pspec_M = np.zeros((ps_list.shape[0],ps_list.shape[0],ps_list.shape[-1]))
    for i in np.arange(ps_list.shape[0]):
        for j in np.arange(ps_list.shape[0]):
            if i<=j:
                pspec_M[i,j] = np.sqrt(ps_list[i]*ps_list[j])
    for i in np.arange(ps_list.shape[0]):
        for j in np.arange(ps_list.shape[0]):
            if j>i:
                if diag:
                    pspec_M[i,j] = np.zeros(len(pspec_M[i,j]))
                pspec_M[i,j] = pspec_M[i,j]
                pspec_M[j,i] = pspec_M[i,j]
                
    return pspec_M

def get_noise_PS_list(chs, noise, p=1, remove = None):
    """Returns a list of noise power spectra at polarization p for channels chs.
    chs: the channels you would like retrieve files of. e.g. ["tube:LC1"], ["tube:LC2","tube:LC3"], ...
    noise: mapsims noise object
    remove: an index or a list thereof of bands you wish to remove. e.g. [5] to remove HF6 when chs are
    p: 0=Temperature, 1=Polarization
    """
    pspec_list = []
    for ch in chs:
        ell_sim, ps_T, ps_P = noise.get_fullsky_noise_spectra(ch[5:])
        if p>0:
            to_c = ps_P
        else:
            to_c = ps_T
        pspec_list.append(to_c[0])
        pspec_list.append(to_c[1])
    pspec_list = np.array(pspec_list)
    inds = np.ones(pspec_list.shape[0])
    if remove!=None:
        inds[np.array(remove)]=0
    
    return pspec_list[(inds.astype(int)).astype(bool)]

def get_noise_PS_matrix(chs, noise, p=1, remove=None, diag = True):
    """Combines the above 2 functions into one"""
    pspec_list = get_noise_PS_list(chs, noise, p=1, remove=remove)
    
    return compute_PS_matrix(pspec_list, diag=diag)

def get_x_from_bands(bands, x = 'freq'):
    """Returns a list of center frequencies of the list of bands bands."""
    res = []
    for b in bands:
        if x=='freq':
            res.append(b.center_frequency.value)
    return np.array(res)
    
def log_likelihood(v_par, M_A, M_noise, M_std, nu, ell_c,w_SO, auto_only = True):
    """Computes the Gaussian log likelihood function for a Markov Chain Monte Carlo run for the dust and synchroton model using fgspectra.
    v_par: is of the form (A_s, A_d, alpha_s, alpha_d, beta_s, beta_d, rho_BB).
    M_A: Power spectrum matrix of cls from the simulations. Should be of shape (#_of_freqs,#_of_freqs,#_of_bins)
    M_noise: Noise power spectrum matrix of cls. Should be of shape (#_of_freqs,#_of_freqs,#_of_bins)
    M_std: Std of cls from the simulations. Should be of shape (#_of_freqs,#_of_freqs,#_of_bins)
    nu: List of center frequencies for each band in M_A, M_noise, M_std
    ell_c: The ell values at which bins were taken
    w_SO: Skyfraction of SO maps
    auto_only: Whether one wishes to take into account only the auto spectra or not
    """
    dust_params = dict(nu=nu, beta= v_par[5], temp=20., nu_0=353.)
    sync_params = dict(nu=nu, beta= v_par[4], nu_0=23.)
    frequency_params = dict(kwseq=(dust_params, sync_params))
    
    try:
        rho = v_par[6]
    except:
        rho=0.045
    
    power_params = dict(
        ell= ell_c,
        alpha=np.array([v_par[3], v_par[2]]),  # +2 to (almost) get D_ell
        ell_0=84,
        amp=np.array([v_par[1]*_rj2cmb(353.)**2 , v_par[0]*_rj2cmb(23.)**2])#*(10**0.5)
        ,rho=rho
    )

    dust_sync = fgc.CorrelatedDustSynchrotron()
    T = dust_sync(frequency_params, power_params)*w_SO+M_noise
    
    coeff_mat = (np.identity(len(nu)) + (1.-auto_only)*np.ones((len(nu),len(nu)))).reshape(len(nu),len(nu),1)
    M_sig = M_std
    res_mat_1 = ((M_A-T)**2)/(M_std**2)
    res_mat = coeff_mat*(res_mat_1)
    if auto_only==False:
        res_f = -1/4*np.sum(res_mat)
    else:
        tmp = np.sum(res_mat[np.identity(res_mat.shape[0]).astype(bool)],axis=1)
        res_f = -1/4*np.sum(tmp*coeff_vec)
    return res_f


def log_prior(v_par):
    """The log prior for v_par.
    v_par: is of the form (A_s, A_d, alpha_s, alpha_d, beta_s, beta_d, rho_BB)."""
    allowed_devs = np.array([0.005, 0.001, 1. , 0.05 , 0.5 , 0.025, 0.1])[:len(v_par)]
    v_par_o =      np.array([0.012 , 0.002,-2.8, -2.35, -3.0, 1.53 , 0.045])[:len(v_par)]
    if (np.abs(v_par_o- v_par)<10*allowed_devs).all():
        return 0.0
    else:
        return -np.inf
    
def log_probability(v_par, M_A,M_noise, M_std, nu, ell_c,w_SO):
    """Input same as log_likelihood"""
    lp = log_prior(v_par)
    ll = log_likelihood(v_par,M_A,M_noise, M_std, nu, ell_c,w_SO)
    
    if not np.isfinite(lp+ll):
        return -np.inf
    return lp + ll

def compute_ps_from_files_auto(sim_num, noise, NSIDE, output_file, pysm_string, chs, lmax, p=2, hitmap = None,remove = None):
    """Used in the case maps do not need to tempered with, and just power spectra are needed. 
    sim_num: number of simulations you would like to load
    noise: the mapsims noise object
    NSIDE: the nside of the stored simulations
    output_file: the directroy from which to retrieve the files
    pysm_string: the pysm_string used to feed into the simulations. e.g. 'd0', 's0', 's0,d0' 
    chs: the channels you would like retrieve files of. e.g. ["tube:LC1"], ["tube:LC2","tube:LC3"],...
    lmax: max ell you wish to compute
    p: polarization of the output. 0 = TT, 1 = EE, 2 = BB
    hitmap: a hitmap you wish for the data to be masked with
    apodize: whether you would like to apodize the maps as they are retrieved
    remove: an index or a list thereof of bands you wish to remove. e.g. [5] to remove HF6 when chs = ["tube:LC1","tube:LC2","tube:LC3"]"""
    
    bands = get_bands(chs,noise, remove = remove)
    output = np.zeros((len(bands),sim_num, lmax+1))
    for j in np.arange(sim_num):
        for k in np.arange(len(bands)):
            tmp = hp.fitsfunc.read_map(output_file +"/"+bands[k].tag+"_" + bands[k].tag[:3]+"_NSIDE_" + str(NSIDE) + "_TAG_" + pysm_string + "_" + str(j) + "_" +".fits", field=(0,1,2))
            if list(hitmap)!=None:
                for h in np.arange(tmp.shape[0]):
                    tmp[h,:] = tmp[h,:]*hitmap
            ps = hp.sphtfunc.anafast(tmp, lmax = lat_lmax)
            output[k,j] = ps[p]
    return output

def compute_ps_from_files_cross(sim_num, noise, NSIDE, output_file, pysm_string, chs, lmax, p=2, hitmap = None,remove = None):
    """Same as compute_ps_from_files_auto but now computes cross-spectra as well"""
    bands = get_bands(chs,noise, remove = remove)
    output = np.zeros((len(bands),len(bands),sim_num, lmax+1))
    for j in np.arange(sim_num):
        lis = []
        for k in np.arange(len(bands)):
            tmp = hp.fitsfunc.read_map(output_file +"/"+bands[k].tag+"_" + bands[k].tag[:3]+"_NSIDE_" + str(NSIDE) + "_TAG_" + pysm_string + "_" + str(j) + "_" +".fits", field=(0,1,2))
            if list(hitmap)!=None:
                    for h in np.arange(tmp.shape[0]):
                        tmp[h,:] = tmp[h,:]*hitmap
            lis.append(tmp)
            for k_p in np.arange(k):                
                ps = hp.sphtfunc.anafast(tmp, lis[k_p],lmax = lat_lmax)
                output[k, k_p, j] = ps[p]
                output[k_p, k, j] = ps[p]
    return output  