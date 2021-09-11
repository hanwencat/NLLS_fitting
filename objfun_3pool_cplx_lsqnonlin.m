function yout = objfun_3pool_cplx_lsqnonlin(p, t, s)
% function yout = objfun_3pool_cplx_lsqnonlin(p, t, s)
% Function for complex 3 pool WM model fitting with real/image splitted
% yout: T2/T2* decay.
% p: parameter array of this model
% t: TE, in unit of ms
% s: The measured decay signal
% Note that p4-6 are in unit of ms as TE, but p7-9 are in unit of Hz, so
% the t associated with cos and sin are converted to unit of second

temp = zeros(length(t),2);

temp1 = p(1) * exp(-t/p(4)) .* cos(-2*pi*p(7)*t/1000) ...
    + p(2) * exp(-t/p(5)) .* cos(-2*pi*p(8)*t/1000) ...
    + p(3) * exp(-t/p(6)) .* cos(-2*pi*p(9)*t/1000);

temp2 = p(1) * exp(-t/p(4)) .* sin(-2*pi*p(7)*t/1000) ...
    + p(2) * exp(-t/p(5)) .* sin(-2*pi*p(8)*t/1000) ...
    + p(3) * exp(-t/p(6)) .* sin(-2*pi*p(9)*t/1000);

temp(:,1) = temp1 * cos(-p(10)) - temp2 * sin(-p(10));
temp(:,2) = temp1 * sin(-p(10)) + temp2 * cos(-p(10));

yout = temp - [real(s(:)), imag(s(:))];
yout = yout(:);
