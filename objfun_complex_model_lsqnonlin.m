function yout = objfun_complex_model_lsqnonlin(p, t, s)

% p: parameter array of this model
% t: TE, in unit of ms
% s: The measured decay signal

temp = (p(1) * exp(-(1/p(4) + 1j*2*pi*p(7)).*t)...
    + p(2) * exp(-(1/p(5) + 1j*2*pi*p(8)).*t)...
    + p(3) * exp(-(1/p(6) + 1j*2*pi*p(9)).*t)) * exp(-1j*p(10));

yout = abs(temp - s);