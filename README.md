# ppta_dr2_noise_analysis

[arXiv:2010.06109](https://arxiv.org/abs/2010.06109) | [MNRAS](https://academic.oup.com/mnras/article-abstract/502/1/478/5957533)

Tools to reproduce the PPTA DR2 noise analysis of pulsar timing residuals. The code is based on [enterprise_warp](https://github.com/bvgoncharov/enterprise_warp "enterprise_warp: Wrapper and tools for Enterprise").

To perform parameter estimation with the most recent PPTA DR2 noise models, just run:
```
python run_dr2.py --prfile params/dr2_parameter_estimation.dat
```

## Measurements of power-law noise term parameters

To obtain values from Figure 1 of the publication without running the analysis, maximum-aposteriori values along with 1-sigma credible levels are provided in the folder `/reproduce_figures/` in `.json` format.

![Figure 1 from the publication (arXiv:2010.06109)](https://github.com/bvgoncharov/ppta_dr2_noise_analysis/blob/master/reproduce_figures/figure_1.jpg "The distribution of spin noise, band noise, system noise, chromatic noise")

## Citation

Goncharov, B., Reardon, D.J., Shannon, R.M., Zhu, X.J., Thrane, E., Bailes, M., Bhat, N.D.R., Dai, S., Hobbs, G., Kerr, M. and Manchester, R.N., 2020. Identifying and mitigating noise sources in precision pulsar timing data sets. Monthly Notices of the Royal Astronomical Society.

> @ARTICLE{2020MNRAS.tmp.3250G,\
> &nbsp;&nbsp;&nbsp;&nbsp;author = {{Goncharov}, Boris and {Reardon}, D.~J. and {Shannon}, R.~M. and {Zhu}, Xing-Jiang and {Thrane}, Eric and {Bailes}, M. and {Bhat}, N.~D.~R. and {Dai}, S. and {Hobbs}, G. and {Kerr}, M. and {Manchester}, R.~N. and {Os{\l}owski}, S. and {Parthasarathy}, A. and {Russell}, C.~J. and {Spiewak}, R. and {Thyagarajan}, N. and {Wang}, J.~B.},\
> &nbsp;&nbsp;&nbsp;&nbsp;title = "{Identifying and mitigating noise sources in precision pulsar timing data sets}",\
> &nbsp;&nbsp;&nbsp;&nbsp;journal = {\mnras},\
> &nbsp;&nbsp;&nbsp;&nbsp;keywords = {stars: neutron, pulsars: general, methods: data analysis, Astrophysics - High Energy Astrophysical Phenomena, Astrophysics - Instrumentation and Methods for Astrophysics, General Relativity and Quantum Cosmology},\
> &nbsp;&nbsp;&nbsp;&nbsp;year = 2020,\
> &nbsp;&nbsp;&nbsp;&nbsp;month = nov,\
> &nbsp;&nbsp;&nbsp;&nbsp;doi = {10.1093/mnras/staa3411},\
> &nbsp;&nbsp;&nbsp;&nbsp;archivePrefix = {arXiv},\
> &nbsp;&nbsp;&nbsp;&nbsp;eprint = {2010.06109},\
> &nbsp;&nbsp;&nbsp;&nbsp;primaryClass = {astro-ph.HE},\
> }
