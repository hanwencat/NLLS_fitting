function yout = objfun_complex_model_lsqnonlin(p, t, s)

% p: parameter array of this model
% t: TE, in unit of s
% s: The measured decay signal

% without T1 compensation
temp = (p(1) * exp(-(1/p(4) + 1j*2*pi*p(7)).*t)...
    + p(2) * exp(-(1/p(5) + 1j*2*pi*p(8)).*t)...
    + p(3) * exp(-(1/p(6) + 1j*2*pi*p(9)).*t)) * exp(-1j*p(10));

% with T1 compensation
% alpha = (22/180)*pi;
% TR = 0.1;
% T1_my = 0.5;
% T1_ie = 1;
% T1_ax = 1;
% 
% temp = (p(1) * exp(-(1/p(4) + 1j*2*pi*p(7)).*t) * sin(alpha)*(1-exp(-TR/T1_my))/(1-cos(alpha)*exp(-TR/T1_my))...
%     + p(2) * exp(-(1/p(5) + 1j*2*pi*p(8)).*t) * sin(alpha)*(1-exp(-TR/T1_ie))/(1-cos(alpha)*exp(-TR/T1_ie))...
%     + p(3) * exp(-(1/p(6) + 1j*2*pi*p(9)).*t) * sin(alpha)*(1-exp(-TR/T1_ax))/(1-cos(alpha)*exp(-TR/T1_ax)))...
%     * exp(-1j*p(10)) ;


yout = abs(temp - s);