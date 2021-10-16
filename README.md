# NLLS fitting
Nonlinear least squares fitting for myelin water imaging data acquired by mGRE sequence (3-pool complex model) 

## How to use
- Run *preprocessing.m* first and then run *ifield_fit.m*.
- Fitting output: amplitude, t2star, and frequency shift for each water pool, and one overall initial phase. 
- Initialization parameters of the fitting may need to be tuned accordingly for best performance.
- The 3-pool complex model is desribed in this paper: https://doi.org/10.1016/j.neuroimage.2015.03.081

## File structure
- ***preprocessing.m***: correct raw bipolar mGRE k-space data, produce voxel-wise complex data, fit total field map and R2* map.
- ***ifield_fit.m***: main fitting program, require preprocessed inputs: complex signal, total field map, and R2* map.
- ***objfun_complex_model.m***: objective function for a 3-pool complex model (default in *ifield_fit.m*).
- ***objfun_magnitude_model.m***: objective function for a 3-pool magnitude model.
- ***phantom_make***: function to produce a 2D computational phantom.
- ***phantom_fit***: script to test the fitting performance using the computational phantom.
