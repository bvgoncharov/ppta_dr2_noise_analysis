import numpy as np
import enterprise.constants as const
from enterprise.signals import signal_base
from enterprise_extensions import models
import enterprise.signals.parameter as parameter
import enterprise.signals.gp_signals as gp_signals
import enterprise.signals.deterministic_signals as deterministic_signals
import enterprise.signals.selections as selections
import enterprise.signals.utils as utils

from enterprise_warp.enterprise_models import StandardModels
from enterprise_warp.enterprise_models import selection_factory

class PPTADR2Models(StandardModels):
  """
  Please follow this example to add your own models for enterprise_warp.
  """
  def __init__(self,psr=None,params=None):
    super(PPTADR2Models, self).__init__(psr=psr,params=params)
    self.priors.update({
      "fd_sys_slope_range": 1e-7,
      "event_j0437_t0": [57050., 57150.],
      "event_j1643_t0": [57050., 57150.],
      "event_j1713_1_t0": [54500., 54900.],
      "event_j1713_2_t0": [57500., 57520.],
      "event_j2145_t0": [56100., 56500.],
      "event_j1603_t0": [53710., 54070.]
    })

  def dm_annual(self, option="default"):
    if option=="default": idx = 2
    return models.dm_annual_signal(idx=idx)
  
  def fd_sys_g(self,option=[]):
  
    idx = 1 # fitting linear trend 
    for ii, fd_sys_term in enumerate(option):
      name = 'fd' + str(idx) + '_sys_' + fd_sys_term
      slope = parameter.Uniform(-self.params.fd_sys_slope_range,\
                                 self.params.fd_sys_slope_range)
      wf = fd_system(slope = slope, idx_fd = idx)
  
      selection_function_name = 'fd_system_selection_' + \
                                 str(self.sys_noise_count)
      setattr(self, selection_function_name,
              selection_factory(selection_function_name))
      if fd_sys_term == "WBCORR_10CM_512":
        self.psr.sys_flags.append('beconfig')
        self.psr.sys_flagvals.append('wbb256_512_256_3p_b')
      else:
        self.psr.sys_flags.append('group')
        self.psr.sys_flagvals.append(fd_sys_term)
  
      fd_sys_term = deterministic_signals.Deterministic( wf, name=name,
                    selection=selections.Selection(\
                    self.__dict__[selection_function_name]) )
  
      if ii == 0:
        fd_sys = fd_sys_term
      elif ii > 0:
        fd_sys += fd_sys_term
  
      self.sys_noise_count += 1
  
    return fd_sys

  def j0437_event(self, option="exp_dip"):
    return dm_exponential_dip(self.params.event_j0437_t0[0],
                              self.params.event_j0437_t0[1],
                              idx="vary", tau_min_10_pow=5, tau_max_10_pow=100)

  def j1713_event_1(self, option="exp_dip"):
    return dm_exponential_dip(self.params.event_j1713_1_t0[0],
                              self.params.event_j1713_1_t0[1], idx=2,
                              tau_min_10_pow=5, tau_max_10_pow=1000,
                              name='dmexp_1')

  def j1713_event_2(self, option="exp_dip"):
    return dm_exponential_dip(self.params.event_j1713_2_t0[0],
                              self.params.event_j1713_2_t0[1], idx="vary",
                              tau_min_10_pow=5, tau_max_10_pow=100,
                              name='dmexp_2')  
  
  def j1643_event(self, option="exp_dip"):
    return dm_exponential_dip(self.params.event_j1643_t0[0],
                              self.params.event_j1643_t0[1],
                              idx="vary", tau_min_10_pow=5, tau_max_10_pow=1000)
  
  def j2145_event(self, option="exp_dip"):
    return dm_exponential_dip(self.params.event_j2145_t0[0],
                              self.params.event_j2145_t0[1],
                              idx="vary", tau_min_10_pow=5, tau_max_10_pow=1000)
  
  def j1603_event(self, option="gaussian_bump"):
    return dm_gaussian_bump(self.params.event_j1603_t0[0],
                            self.params.event_j1603_t0[1], idx=2)

  def paired_ppta_band_noise(self, option=[]):
    """
    Including band noise terms for paired values of the PPTA "-B" flag:
    joint 1020-cm and 4050-cm red processes.
    """
    for ii, paired_band_term in enumerate(option):

      log10_A = parameter.Uniform(self.params.syn_lgA[0],self.params.syn_lgA[1])
      gamma = parameter.Uniform(self.params.syn_gamma[0],\
                                self.params.syn_gamma[1])
      pl = utils.powerlaw(log10_A=log10_A, gamma=gamma, \
                          components=self.params.red_general_nfouriercomp)
  
      setattr(self, paired_band_term, globals()[paired_band_term])
      nfreqs = self.determine_nfreqs(sel_func_name=paired_band_term)
  
      pbn_term = gp_signals.FourierBasisGP(spectrum=pl, Tspan=self.params.Tspan,
                                        name='band_noise_' + paired_band_term,
                                        selection=selections.Selection( \
                                        self.__dict__[paired_band_term] ),
                                        components=nfreqs)
      if ii == 0:
        pbn = pbn_term
      elif ii > 0:
        pbn += pbn_term

    return pbn

# PPTA DR2 signal models

@signal_base.function
def fd_system(freqs, slope = 1e-7, idx_fd = 1):
    freq_median = freqs - np.median(freqs)
    return np.sign(freq_median)**(idx_fd + 1) * slope * freq_median**idx_fd

@signal_base.function
def chrom_gaussian_bump(toas, freqs, log10_Amp=-2.5, sign_param=1.0,
                    t0=53890, sigma=81, idx=2):
    """
    Chromatic time-domain Gaussian delay term in TOAs.
    Example: J1603-7202 in Lentati et al, MNRAS 458, 2016.
    """
    t0 *= const.day
    sigma *= const.day
    wf = 10**log10_Amp * np.exp(-(toas - t0)**2/2/sigma**2)
    return np.sign(sign_param) * wf * (1400 / freqs) ** idx

# PPTA DR2 signal wrappers

def dm_exponential_dip(tmin, tmax, idx=2, sign='negative', name='dmexp',
                       lgA_min=-10., lgA_max=-2., sign_min=-1., sign_max=1.,
                       tau_min_10_pow=5, tau_max_10_pow=100):
    t0_dmexp = parameter.Uniform(tmin,tmax)
    log10_Amp_dmexp = parameter.Uniform(lgA_min, lgA_max)
    log10_tau_dmexp = parameter.Uniform(np.log10(tau_min_10_pow),
                                        np.log10(tau_max_10_pow))
    if idx=='vary':
        idx = parameter.Uniform(-7, 7)
    if sign == 'vary':
        sign_param = parameter.Uniform(sign_min, sign_max)
    elif sign == 'positive':
        sign_param = 1.0
    else:
        sign_param = -1.0
    wf = models.chrom_exp_decay(log10_Amp=log10_Amp_dmexp,
                                t0=t0_dmexp, log10_tau=log10_tau_dmexp,
                                sign_param=sign_param, idx=idx)
    dmexp = deterministic_signals.Deterministic(wf, name=name)

    return dmexp

def dm_gaussian_bump(tmin, tmax, idx=2, sigma_min=20, sigma_max=140,
    log10_A_low=-6, log10_A_high=-1, name='dm_bump'):
    """
    Returns chromatic Gaussian bump (i.e. TOA advance):
    :param tmin, tmax:
        search window for exponential cusp time.
    :param idx:
        index of radio frequency dependence (i.e. DM is 2). If this is set
        to 'vary' then the index will vary from 1 - 6
    :param sigma_min, sigma_max:
        standard deviation of a Gaussian in MJD
    :param sign:
        [boolean] allow for positive or negative exponential features.
    :param name: Name of signal
    :return dm_bump:
        chromatic Gaussian bump waveform.
    """
    sign_param = 1.0
    t0_dm_bump = parameter.Uniform(tmin,tmax)
    sigma_dm_bump = parameter.Uniform(sigma_min,sigma_max)
    log10_Amp_dm_bump = parameter.Uniform(log10_A_low, log10_A_high)
    if idx == 'vary':
        idx = parameter.Uniform(0, 6)
    wf = chrom_gaussian_bump(log10_Amp=log10_Amp_dm_bump,
                         t0=t0_dm_bump, sigma=sigma_dm_bump,
                         sign_param=sign_param, idx=idx)
    dm_bump = deterministic_signals.Deterministic(wf, name=name)

    return dm_bump

# PPTA DR2 selections

def by_B_4050CM(flags):
    """Selection function to split by PPTA frequency band under -B flag"""
    sel_40b = np.char.lower(flags['B']) == '40cm'
    sel_50b = np.char.lower(flags['B']) == '50cm'
    return {'4050CM': sel_40b + sel_50b}

def by_B_1020CM(flags):
    """Selection function to split by PPTA frequency band under -B flag"""
    sel_10b = np.char.lower(flags['B']) == '10cm'
    sel_20b = np.char.lower(flags['B']) == '20cm'
    return {'1020CM': sel_10b + sel_20b}
