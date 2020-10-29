# ppta_dr2_noise_analysis
Tools to reproduce the PPTA DR2 noise analysis of pulsar timing residuals. The code is based on [enterprise_warp](https://github.com/bvgoncharov/enterprise_warp "enterprise_warp: Wrapper and tools for Enterprise").

To perform parameter estimation with the most recent PPTA DR2 noise models, just run:
```
python run_dr2.py --prfile params/dr2_parameter_estimation.dat
```

## Citation

> @article{goncharov2020pptadr2noise,\
> &nbsp;&nbsp;title={Identifying and mitigating noise sources in precision pulsar timing data sets},\
> &nbsp;&nbsp;author={Goncharov, Boris and Reardon, DJ and Shannon, RM and Zhu, Xing-Jiang and Thrane, Eric and Bailes, M and Bhat, NDR and Dai, S and Hobbs, G and Kerr, M and others},\
> &nbsp;&nbsp;journal={arXiv preprint arXiv:2010.06109},\
> &nbsp;&nbsp;year={2020}\
> }

(Accepted in MNRAS, the citation will be updated soon)
