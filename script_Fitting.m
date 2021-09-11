fname_img = 'imgCoilOffsetCorr.mat';
fname_F = 'F_VSFFunction.mat';
fname_b = 'B0_midResults.mat';

used_TE = 16;

s = load(fname_img);
img = s.imgAllEcho_comb;

s = load(fname_F);
Fx = s.Fx;
Fy = s.Fy;
Fz = s.Fz;

load(fname_b, 'b'); % Field map in unit of T
gammabar = 42577482.5; % in unit of Hz/T
b_hz = b * gammabar; % Field map in unit of Hz

clear s;

mask = niftiread('img3_mag_comb_brain_mask.nii.gz');

if size(img(:,:,:,1)) ~= size(mask)
    error('Image and Mask should have the same size.');
end

xpt = size(mask,1);
ypt = size(mask,2);
zpt = size(mask,3);

TE = 2.5 : 1.5 : 49.0; % In unit of ms

if used_TE > length(TE)
    error('TE used for fitting cannot larger than total TE number.');
end

t_used = TE(1:used_TE);

tot_slice = 40;
% tot_slice = 1:zpt;

imgVSFCorr = abs(img)./ Fz;
imgVSFCorr(isnan(imgVSFCorr)) = 0;
imgVSFCorr(isinf(imgVSFCorr)) = 0;
imgVSFCorr = imgVSFCorr .* exp(1j*angle(img));
% clear img;

% objfcn = @(p,t) (p(1) * exp(-1/p(4)*t - 1j*2*pi*p(7)*t) ...
%     + p(2) * exp(-1/p(5)*t + 1j*2*pi*p(8)*t) ...
%     + p(3) * exp(-1/p(6)*t + 1j*2*pi*p(9)*t)) * exp(-1j*p(10));

% 10 parameters + 1 resnorm
fitted_param = zeros(xpt, ypt, length(tot_slice), 11);

% opts = optimoptions(@lsqcurvefit, 'Algorithm', 'trust-region-reflective',...
%     'Display', 'off', 'FunctionTolerance', 1e-5, 'StepTolerance', 1e-5);

opts = optimoptions(@lsqnonlin, 'Algorithm', 'trust-region-reflective',...
    'Display', 'off', 'FunctionTolerance', 1e-5, 'StepTolerance', 1e-5);

%% Start Fitting
disp('Started Fitting.');
tic;
parpool(4);
s = tot_slice;
% for s = tot_slice
    parfor x = 1 : xpt
        for y = 1 : ypt
            
            if mask(x,y,s) == 0
                continue;
            end
            
            decay = squeeze(imgVSFCorr(x, y, s, :));
            decay = decay / max(abs(decay));
            decay(isnan(decay)) = 0;
            decay(isinf(decay)) = 0;
            decay = double(decay(1:used_TE));
            decay = abs(decay).*exp(1j*unwrap(angle(decay)));
                                  
            p0 =...
                [0.1*abs(decay(1)); 0.6*abs(decay(1)); 0.4*abs(decay(1)); ...
                10; 64; 48;...
                b_hz(x,y,s); b_hz(x,y,s); b_hz(x,y,s); ...
                angle(decay(1))];
            
            p_lower = ...
                [0; 0; 0; ...
                3; 24; 24; ...
                b_hz(x,y,s) - 75; b_hz(x,y,s) - 25; b_hz(x,y,s) - 25; ...
                -pi];
            
            p_upper = ...
                [0.1*abs(decay(1)); 0.6*abs(decay(1)); 0.4*abs(decay(1)); ...
                24; 1000; 1000; ...
                b_hz(x,y,s) + 75; b_hz(x,y,s) + 25; b_hz(x,y,s) + 25; ...
                pi];
            
%             decay = [real(decay), imag(decay)]; % Split decay into real and imag parts
            
%             [p_est,resnorm] = lsqcurvefit(objfcn, p0, t_used(:), decay, p_lower, p_upper, opts);

%             [p_est,resnorm] = lsqcurvefit(@objfun_3pool_cplx, p0, t_used(:), decay, p_lower, p_upper, opts);
            
            objfcn = @(p)objfun_3pool_cplx_lsqnonlin(p,t_used,decay);
            [p_est,resnorm] = lsqnonlin(objfcn, p0, p_lower, p_upper, opts);

            fitted_param(x, y, :) = reshape([p_est; resnorm], [1,1,11]);
            
        end
    end
% end
delete(gcp('nocreate'));

toc;