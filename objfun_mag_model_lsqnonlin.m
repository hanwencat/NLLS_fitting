function yout = objfun_mag_model_lsqnonlin(p, t, s)
% function yout = objfun_mag_model_lsqnonlin(p, t, s)
% Function for magnitude 3 pool WM model fitting
% yout: T2/T2* decay.
% p: parameter array of this model
% t: TE, in unit of ms
% s: The magnitude of the measured decay signal
% Note that p4-6 are in unit of ms as TE



temp = p(1) * exp(-t/p(4)) + p(2) * exp(-t/p(5)) + p(3) * exp(-t/p(6));

yout = temp - s;
%yout = yout(:);
