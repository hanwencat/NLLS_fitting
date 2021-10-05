%% Generate the phantom
% parameters for making the phantom
snr_range = [50, 500];
mwf_range = [0, 0.5];
t2s = [0.01,0.064,0.048]; % unit: seconds
x_dim = 100;
y_dim = 100; 
fs_mu = [-6,0,0];
% fs_sigma = [2,2,1]; 
%fs_mu = [0,0,0];
fs_sigma = [0,0,0]; 
echo_time = TE; 
% echo_time = 0.002:0.002:0.048;

% make the phantom
[phantom.signal, phantom.mwf, phantom.noise, phantom.fs_my, phantom.fs_ax, phantom.fs_ex, phantom.phi0] = ...
phantom_make(snr_range, mwf_range, t2s, x_dim, y_dim, fs_mu, fs_sigma, echo_time);


%% NLLS fitting
% fitting parameter for the 3-pool magnitude model
p0 = [0.2; 0.5; 0.3; 0.01; 0.064; 0.048]; % intialization value (unit: seconds)
p_lower = [0; 0; 0; 0.003; 0.025; 0.025]; % lower bound
p_upper = [1; 1; 1; 0.024; 1; 1]; % upper bound

% placeholder for the fitted parameters (6 parameters + 1 resnorm)
fitted_param = zeros(x_dim, y_dim, 7);

% denoise the phantom signal by mppca algorithm
%phantom_denoised = denoise(abs(phantom.signal),[5,5]);

% optimization parameters
opts = optimoptions(@lsqnonlin, 'Algorithm', 'trust-region-reflective',...
    'Display', 'off', 'FunctionTolerance', 1e-5, 'StepTolerance', 1e-5,... 
    'MaxIterations',1000, 'MaxFunctionEvaluations', 1000);

% start fitting
f = waitbar(0, 'Starting');
for x = 1:x_dim
    for y = 1:y_dim
        decay = squeeze(abs(phantom.signal(x,y,:))); % squeeze the dimension
        %decay = squeeze(phantom_denoised(x,y,:));
        %decay = decay / max(decay);
        
        % create object function and fit
        objfcn = @(p)objfun_mag_model_lsqnonlin(p, echo_time, decay);
        [p_est, resnorm] = lsqnonlin(objfcn, p0, p_lower, p_upper, opts);
        fitted_param(x,y,:) = reshape([p_est; resnorm], [1,1,7]);
        waitbar(x/x_dim, f, sprintf('Progress: %d %%', floor(x/x_dim*100)));
    end
end
close(f);

% produce MWF map
fitted_mwf = fitted_param(:,:,1)./(fitted_param(:,:,1)+fitted_param(:,:,2)+fitted_param(:,:,3));

% plot MWF maps
figure, subplot(131);imagesc(snr_range, mwf_range, phantom.mwf, [0,0.5]); xlabel('SNR'); ylabel('MWF');
colormap(jet); colorbar; set(colorbar, 'YDir', 'reverse'); title('Ground truth MWF vs SNR'); axis square;
subplot(132);imagesc(snr_range, mwf_range, fitted_mwf, [0,0.5]); xlabel('SNR'); ylabel('MWF');
colormap(jet); colorbar; set(colorbar, 'YDir', 'reverse'); title('Fitted MWF vs SNR');axis square;
subplot(133);imagesc(snr_range, mwf_range, (fitted_mwf - phantom.mwf), [-0.1,0.1]); xlabel('SNR'); ylabel('MWF');
colormap(jet); colorbar; set(colorbar, 'YDir', 'reverse'); title('MWF error vs SNR');axis square;


    
  