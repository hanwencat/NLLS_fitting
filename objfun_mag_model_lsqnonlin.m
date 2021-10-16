function yout = objfun_mag_model_lsqnonlin(p, t, s)
% objective function for NLLS fitting
% p: parameter array of the 3 pool magnitude model
% p(1-3): amplitudes of the 3 pools
% p(4-6): t2stars of the 3 pools
% t: echo times, in unit of s
% s: The measured decay signals


% % Without T1 compensation
fit = p(1) * exp(-t/p(4)) + p(2) * exp(-t/p(5)) + p(3) * exp(-t/p(6));

% % With T1 compensation
% alpha = (22/180)*pi;
% TR = 0.1;
% T1_my = 0.3;
% T1_ie = 0.8;
% T1_ax = 1;
% temp = p(1) * exp(-t/p(4)) * sin(alpha)*(1-exp(-TR/T1_my))/(1-cos(alpha)*exp(-TR/T1_my)) ... 
% + p(2) * exp(-t/p(5)) * sin(alpha)*(1-exp(-TR/T1_ie))/(1-cos(alpha)*exp(-TR/T1_ie)) ...
% + p(3) * exp(-t/p(6)) * sin(alpha)*(1-exp(-TR/T1_ax))/(1-cos(alpha)*exp(-TR/T1_ax));

yout = fit - s;



% % Incorporate regularizations
% lambda_l1 = 0.001;
% lambda_l2 = 0;
% reg_l1 = abs(p(1)) + abs(p(2)) + abs(p(3));
% reg_l2 = p(1)^2 + p(2)^2 + p(3)^2;
% yout = [yout; lambda_l1*reg_l1 + lambda_l2*reg_l2];
% yout = [reshape(yout,[],1); lambda*reg];

