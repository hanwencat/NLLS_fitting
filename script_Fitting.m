% fname_img = 'imgCoilOffsetCorr.mat';
% fname_F = 'F_VSFFunction.mat';
% fname_b = 'B0_midResults.mat';

used_te = 12;

% s = load(fname_img);
% img = s.imgAllEcho_comb;
%img = abs(iField);
img = abs(iField(:,:,:,1:2:end));


% s = load(fname_F);
% Fx = s.Fx;
% Fy = s.Fy;
% Fz = s.Fz;
Fz = 1;

% load(fname_b, 'b'); % Field map in unit of T
% gammabar = 42577482.5; % in unit of Hz/T
% b_hz = b * gammabar; % Field map in unit of Hz
% b_hz = tfs*1000/(2*pi);

% clear s;

% mask = niftiread('img3_mag_comb_brain_mask.nii.gz');
mask = Mask;

if size(img(:,:,:,1)) ~= size(mask)
    error('Image and Mask should have the same size.');
end

xpt = size(mask,1);
ypt = size(mask,2);
zpt = size(mask,3);

%te = TE*1000; % In unit of ms
te = TE(1:2:end)*1000;

if used_te > length(te)
    error('te used for fitting cannot larger than total te number.');
end

t_used = te(1:used_te);

% tot_slice = 40;
tot_slice = 1:zpt;

imgVSFCorr = img./ Fz;
imgVSFCorr(isnan(imgVSFCorr)) = 0;
imgVSFCorr(isinf(imgVSFCorr)) = 0;
%imgVSFCorr = imgVSFCorr .* exp(1j*angle(img));
%imgVSFCorr = imgVSFCorr .* exp(1j*iFreq_all_echo_corr);
%imgVSFCorr = imgVSFCorr .* exp(1j*iFreq_all_echo_corr(:,:,:,1:2:end));

% clear img;

% objfcn = @(p,t) (p(1) * exp(-1/p(4)*t - 1j*2*pi*p(7)*t) ...
%     + p(2) * exp(-1/p(5)*t + 1j*2*pi*p(8)*t) ...
%     + p(3) * exp(-1/p(6)*t + 1j*2*pi*p(9)*t)) * exp(-1j*p(10));

% 6 parameters + 1 resnorm
fitted_param = zeros(xpt, ypt, length(tot_slice), 7);

% opts = optimoptions(@lsqcurvefit, 'Algorithm', 'trust-region-reflective',...
%     'Display', 'off', 'FunctionTolerance', 1e-5, 'StepTolerance', 1e-5);

opts = optimoptions(@lsqnonlin, 'Algorithm', 'trust-region-reflective',...
    'Display', 'off', 'FunctionTolerance', 1e-5, 'StepTolerance', 1e-5, 'MaxIterations',1000, 'MaxFunctionEvaluations', 1000);

%% Start Fitting
disp('Started Fitting.');
tic;
parpool(8);
% s = tot_slice(40);

f = waitbar(0, 'Starting');
for s = 1 : zpt
    parfor x = 1 : xpt
        for y = 1 : ypt
            
            if mask(x,y,s) == 0
                continue;
            end
            
            decay = squeeze(imgVSFCorr(x, y, s, :));
            decay = decay / max(decay);
            decay(isnan(decay)) = 0;
            decay(isinf(decay)) = 0;
            decay = double(decay(1:used_te));
            %decay = abs(decay).*exp(1j*unwrap(angle(decay)));
                                  
            % 3-pool magnitude model
            p0 = [0.1*decay(1); 0.6*decay(1); 0.3*decay(1); ...
                10; 64; 48];
            
            p_lower = [0; 0; 0; ...
                3; 25; 25];
                
            p_upper = [5*decay(1); 5*decay(1); 5*decay(1); ...
                24; 2000; 2000];
                        
            
            
%             % 2-pool magnitude model
%             p0 = [0.3*abs(decay(1)); 0.7*abs(decay(1)); 0; ...
%                 10; 70; 0];
%             
%             p_lower = [0; 0; 0; ...
%                 1; 40; 0];
%                 
%             p_upper = [1*abs(decay(1)); 1*abs(decay(1)); 0; ...
%                 40; 2000; 0];
            
            
            
            objfcn = @(p)objfun_mag_model_lsqnonlin(p,t_used,decay);
            [p_est,resnorm] = lsqnonlin(objfcn, p0, p_lower, p_upper, opts);

            fitted_param(x, y, s, :) = reshape([p_est; resnorm], [1,1,1,7]);
            
        end
    end
    
    waitbar(s/zpt, f, sprintf('Progress: %d %%', floor(s/zpt*100)));
    %pause(0.1);
end
close(f);
delete(gcp('nocreate'));

toc;


mwf_map = fitted_param(:,:,:,1)./(fitted_param(:,:,:,1)+fitted_param(:,:,:,2)+fitted_param(:,:,:,3));
figure, subplot(131);imshow(squeeze(mwf_map(70,:,:)),[0,0.4],'border','tight'); colorbar;
subplot(132);imshow(squeeze(mwf_map(:,50,:)),[0,0.4],'border','tight'); colorbar;
subplot(133);imshow(squeeze(mwf_map(:,:,60)),[0,0.4],'border','tight'); colorbar;
% title_txt = ['number of echoes: ', num2str(used_te), ', fitting parameter (upper bound): ', num2str(2000)];
% suptitle(title_txt);