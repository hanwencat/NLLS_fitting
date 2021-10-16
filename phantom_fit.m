%% Generate the phantom
% parameters for making a 2D phantom
snr_range = [50, 500];
mwf_range = [0, 0.5];
t2s = [0.01,0.064,0.048]; % unit: seconds
x_dim = 100;
y_dim = 100; 
fs_mu = [20,-7,0];
%fs_sigma = [2,2,1]; 
%fs_mu = [0,0,0];
fs_sigma = [0,0,0]; 
%echo_time = TE; 
echo_time = (0.002:0.002:0.048)';

% make the phantom
[phantom.signal, phantom.mwf, phantom.noise, phantom.fs_my, phantom.fs_ax, phantom.fs_ex, phantom.phi0] = ...
phantom_make(snr_range, mwf_range, t2s, x_dim, y_dim, fs_mu, fs_sigma, echo_time);


% %% NLLS magnitude fitting for the produced phantom
% % fitting parameters for the 3-pool magnitude model (initial value, lower and upper bounds)
% % parameter order: 3 amplitudes, 3 t2stars
% p0 = [0.2; 0.5; 0.3; 0.01; 0.064; 0.048]; % intialization value (unit: seconds)
% p_lower = [0; 0; 0; 0.003; 0.025; 0.025]; % lower bound
% p_upper = [1; 1; 1; 0.024; 0.15; 0.15]; % upper bound
% 
% % placeholder for the fitted parameters (6 parameters + 1 resnorm)
% fitted_param = zeros(x_dim, y_dim, 7);
% 
% % denoise the phantom signal by mppca algorithm
% %phantom_denoised = denoise(abs(phantom.signal),[5,5]);
% 
% % optimization parameters
% opts = optimoptions(@lsqnonlin, 'Algorithm', 'trust-region-reflective',...
%     'Display', 'off', 'FunctionTolerance', 1e-5, 'StepTolerance', 1e-5,... 
%     'MaxIterations',1000, 'MaxFunctionEvaluations', 1000);
% 
% % start fitting
% %f = waitbar(0, 'Starting');
% tic
% signal = phantom.signal;
% parfor x = 1:x_dim
%     for y = 1:y_dim
%         decay = squeeze(abs(signal(x,y,:))); % squeeze the dimension
%         %decay = squeeze(phantom_denoised(x,y,:));
%         decay = decay / max(decay);
%         
%         % create object function and fit
%         objfcn = @(p)objfun_mag_model_lsqnonlin(p, echo_time, decay);
%         [p_est, resnorm] = lsqnonlin(objfcn, p0, p_lower, p_upper, opts);
%         fitted_param(x,y,:) = reshape([p_est; resnorm], [1,1,7]);   
%     end
%     %waitbar(x/x_dim, f, sprintf('Progress: %d %%', floor(x/x_dim*100)));
% end
% toc
% %close(f);


%% NLLS complex fitting for the produced phantom
% fitting parameters for the 3-pool complex model (initial value, lower and upper bounds)
% parameter order: 3 amplitudes, 3 t2stars, 3 frequency shifts, 1 initial phase factor
p0 = [0.2; 0.5; 0.3; 0.01; 0.064; 0.048; 20; -7; 0; 0]; % intialization value (unit: seconds)
p_lower = [0; 0; 0; 0.003; 0.025; 0.025; -10; -20; -10; -pi]; % lower bound
p_upper = [1; 1; 1; 0.024; 0.15; 0.15; 50; 20; 10; pi]; % upper bound

% placeholder for the fitted parameters (10 parameters + 1 resnorm)
fitted_param = zeros(x_dim, y_dim, 11);

% denoise the phantom signal by mppca algorithm
%phantom_denoised = denoise(abs(phantom.signal),[5,5]);

% optimization parameters
opts = optimoptions(@lsqnonlin, 'Algorithm', 'trust-region-reflective',...
    'Display', 'off', 'FunctionTolerance', 1e-5, 'StepTolerance', 1e-5,... 
    'MaxIterations', 1000, 'MaxFunctionEvaluations', 1000);

% start fitting for phantom
%f = waitbar(0, 'Starting'); % wait bar does not work with parfor loop
tic
signal = phantom.signal;
%signal = phantom_denoised .* exp(1j*angle(phantom.signal));
parfor x = 1:x_dim
    for y = 1:y_dim
        decay = squeeze(signal(x,y,:)); % squeeze the dimension
        %decay = squeeze(phantom_denoised(x,y,:));
        %decay = decay / abs(decay(1));
        
        % create object function and fit
        objfcn = @(p)objfun_complex_model_lsqnonlin(p, echo_time, decay);
        [p_est, resnorm] = lsqnonlin(objfcn, p0, p_lower, p_upper, opts);
        fitted_param(x,y,:) = reshape([p_est; resnorm], [1,1,11]); 
    end
    %waitbar(x/x_dim, f, sprintf('Progress: %d %%', floor(x/x_dim*100)));
end
toc
%close(f);


extract MWF values
fitted_mwf = fitted_param(:,:,1)./(fitted_param(:,:,1)+fitted_param(:,:,2)+fitted_param(:,:,3));

% plot MWF maps for the phantom
figure, subplot(431);imagesc(snr_range, mwf_range, phantom.mwf, [0,0.5]); xlabel('SNR'); ylabel('MWF');
colormap(jet); colorbar; set(colorbar, 'YDir', 'reverse'); title('Truth MWF'); axis square;
subplot(432);imagesc(snr_range, mwf_range, fitted_mwf, [0,0.5]); xlabel('SNR'); ylabel('MWF');
colormap(jet); colorbar; set(colorbar, 'YDir', 'reverse'); title('Fitted MWF');axis square;
subplot(433);imagesc(snr_range, mwf_range, (fitted_mwf - phantom.mwf), [-0.1,0.1]); xlabel('SNR'); ylabel('MWF');
colormap(jet); colorbar; set(colorbar, 'YDir', 'reverse'); title('fitting error');axis square;

subplot(434);imagesc(snr_range, mwf_range, phantom.phi0, [-pi,pi]); xlabel('SNR'); ylabel('MWF');
colormap(jet); colorbar; set(colorbar, 'YDir', 'reverse'); title('Truth phi0 (random)'); axis square;
subplot(435);imagesc(snr_range, mwf_range, fitted_param(:,:,10), [-pi,pi]); xlabel('SNR'); ylabel('MWF');
colormap(jet); colorbar; set(colorbar, 'YDir', 'reverse'); title('Fitted phi0');axis square;
subplot(436);imagesc(snr_range, mwf_range, (fitted_param(:,:,10) - phantom.phi0), [-pi,pi]); xlabel('SNR'); ylabel('MWF');
colormap(jet); colorbar; set(colorbar, 'YDir', 'reverse'); title('fitting error');axis square;
    
subplot(437);imagesc(snr_range, mwf_range, phantom.fs_my, [-30,30]); xlabel('SNR'); ylabel('MWF');
colormap(jet); colorbar; set(colorbar, 'YDir', 'reverse'); title('Truth myelin water fs (20 Hz)'); axis square;
subplot(438);imagesc(snr_range, mwf_range, fitted_param(:,:,7), [-30,30]); xlabel('SNR'); ylabel('MWF');
colormap(jet); colorbar; set(colorbar, 'YDir', 'reverse'); title('Fitted myelin water fs');axis square;
subplot(439);imagesc(snr_range, mwf_range, (fitted_param(:,:,7) - phantom.fs_my), [-30,30]); xlabel('SNR'); ylabel('MWF');
colormap(jet); colorbar; set(colorbar, 'YDir', 'reverse'); title('fitting error');axis square;

subplot(4,3,10);imagesc(snr_range, mwf_range, phantom.fs_ax, [-20,20]); xlabel('SNR'); ylabel('MWF');
colormap(jet); colorbar; set(colorbar, 'YDir', 'reverse'); title('Truth axonal water fs (-7 Hz)'); axis square;
subplot(4,3,11);imagesc(snr_range, mwf_range, fitted_param(:,:,8), [-20,20]); xlabel('SNR'); ylabel('MWF');
colormap(jet); colorbar; set(colorbar, 'YDir', 'reverse'); title('Fitted axonal water fs');axis square;
subplot(4,3,12);imagesc(snr_range, mwf_range, (fitted_param(:,:,8) - phantom.fs_ax), [-20,20]); xlabel('SNR'); ylabel('MWF');
colormap(jet); colorbar; set(colorbar, 'YDir', 'reverse'); title('fitting error');axis square;
