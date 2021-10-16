function yout = objfun_complex_model(p, t, s)
% objective function for NLLS fitting
% p: parameter array of the 3 pool complex model
% p(1-3): amplitudes of the 3 pools
% p(4-6): t2stars of the 3 pools
% p(7-9): frequency shifts of the 3 pools
% p(10): initial phase
% t: echo times, in unit of s
% s: The measured decay signals

% % without T1 compensation
fit = (p(1) * exp(-(1/p(4) + 1j*2*pi*p(7)).*t)...
    + p(2) * exp(-(1/p(5) + 1j*2*pi*p(8)).*t)...
    + p(3) * exp(-(1/p(6) + 1j*2*pi*p(9)).*t)) * exp(-1j*p(10));


% % with T1 compensation
% alpha = (22/180)*pi;
% TR = 0.1;
% T1_my = 0.5;
% T1_ie = 1;
% T1_ax = 1;
% 
% fit = (p(1) * exp(-(1/p(4) + 1j*2*pi*p(7)).*t) * sin(alpha)*(1-exp(-TR/T1_my))/(1-cos(alpha)*exp(-TR/T1_my))...
%     + p(2) * exp(-(1/p(5) + 1j*2*pi*p(8)).*t) * sin(alpha)*(1-exp(-TR/T1_ie))/(1-cos(alpha)*exp(-TR/T1_ie))...
%     + p(3) * exp(-(1/p(6) + 1j*2*pi*p(9)).*t) * sin(alpha)*(1-exp(-TR/T1_ax))/(1-cos(alpha)*exp(-TR/T1_ax)))...
%     * exp(-1j*p(10)) ;


yout = abs(fit - s);