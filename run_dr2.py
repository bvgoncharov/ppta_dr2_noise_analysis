#!/bin/python

import numpy as np
import sys
import os
import bilby
import inspect
from enterprise_warp import enterprise_warp
from enterprise_warp import bilby_warp
from enterprise_extensions import model_utils

import ppta_dr2_models

opts = enterprise_warp.parse_commandline()

custom = ppta_dr2_models.PPTADR2Models

params = enterprise_warp.Params(opts.prfile,opts=opts,custom_models_obj=custom)
pta = enterprise_warp.init_pta(params)

if params.sampler == 'ptmcmcsampler':
    super_model = model_utils.HyperModel(pta)
    print('Super model parameters: ', super_model.params)
    sampler = super_model.setup_sampler(resume=True, outdir=params.output_dir)
    N = params.nsamp
    x0 = super_model.initial_sample()

    # Remove extra kwargs that Bilby took from PTSampler module, not ".sample"
    ptmcmc_sample_kwargs = inspect.getargspec(sampler.sample).args
    upd_sample_kwargs = {key: val for key, val in params.sampler_kwargs.items()
                                  if key in ptmcmc_sample_kwargs}
    del upd_sample_kwargs['Niter']
    del upd_sample_kwargs['p0']

    sampler.sample(x0, N, **upd_sample_kwargs)
else:
    priors = bilby_warp.get_bilby_prior_dict(pta[0])
    parameters = dict.fromkeys(priors.keys())
    likelihood = bilby_warp.PTABilbyLikelihood(pta[0],parameters)
    label = os.path.basename(os.path.normpath(params.out))
    if opts.mpi_regime != 1:
      bilby.run_sampler(likelihood=likelihood, priors=priors,
                        outdir=params.output_dir, label=params.label,
                        sampler=params.sampler, **params.sampler_kwargs)
    else:
      print('Preparations for the MPI run are complete - now set \
             opts.mpi_regime to 2 and enjoy the speed!')

print('Finished: ',opts.num)
