%% NLLS complex fitting for mGRE data
% Created by Hanwen liu on 2021-10-15

%%%%%%%%%%%%%%%%%%%%%% NAMING CONVENTION %%%%%%%%%%%%%%%%%%%%%%%%%
% iField - multiecho complex MRI dataset (the 4th dimension is the echo)
% iMag - magnitude image, square root of squares of all echoes
% Mask - binary mask denoting the region of interest
% tfs - total field map (rad/ms), generated from preprocessing
% fs - frequency shift map (Hz)
% R2s - R2* map, generated from preprocessing 

% TE - echo time, unit in sec
% B0_dir - unit vector representing direction of B0 field 
% matrix_size - sizes ([x y z]) of the imaging volume
% voxel_size - size of the voxel unit in mm

% delta_TE - echo spacing, unit in sec
% CF - center frequency, unit in Hz
% B0_strength - magnetic field strength, unit in Tesla
%%%%%%%%%%%%%%%%%%%%%% NAMING CONVENTION %%%%%%%%%%%%%%%%%%%%%%%%%

% load data, echo time, and frequency shift map
data = iField.*Mask;
echo_time = TE;
% data = iField(:,:,:,1:2:end).*Mask; % fit with odd echoes only
% echo_time = TE(1:2:end); % fit with odd echoes only
fs = -tfs(:,:,:)*1000/(2*pi); % unit in Hz

% get dimemsion
data_dim = size(data);
x_dim = data_dim(1);
y_dim = data_dim(2);
z_dim = data_dim(3);

% fitting parameters for the 3-pool complex model (initial value, lower and upper bounds)
% parameter order: 3 amplitudes, 3 t2stars, 3 frequency shifts, 1 initial phase factor

% placeholder for the fitted parameters (10 parameters + 1 resnorm)
fitted_param_cplx = zeros(x_dim, y_dim, z_dim, 11);
fitted_param_mag = zeros(x_dim, y_dim, z_dim, 7);

% optimization parameters
opts = optimoptions(@lsqnonlin, 'Algorithm', 'trust-region-reflective',...
    'Display', 'off', 'FunctionTolerance', 1e-5, 'StepTolerance', 1e-5,... 
    'MaxIterations', 1000, 'MaxFunctionEvaluations', 1000);

% start fitting for the experimental data
f = waitbar(0, 'Fitting in progress');
tic
%signal = phantom_denoised .* exp(1j*angle(phantom.signal));
for x = 1:x_dim
    parfor y = 1:y_dim
        for z = 1:z_dim
            if Mask(x,y,z) == 0 || R2s(x,y,z) == 0 ...
                    || isnan(R2s(x,y,z)) || isnan(fs(x,y,z))
                continue;
            end
            decay = squeeze(data(x,y,z,:)); % squeeze the dimension
            %decay = squeeze(phantom_denoised(x,y,:));
            decay = decay / abs(decay(1));
            
            if any(isnan(decay)) || any(isinf(decay))
                continue;
            end
            
            % disp(num2str(x)+ " " + num2str(y)+ " " + num2str(z));
            
%             % create complex fitting parameters adaptively for each voxel
%             p0_cplx = [0.2; 0.5; 0.3; 0.01; 0.064; 0.048; 
%                 fs(x,y,z); fs(x,y,z); fs(x,y,z); 0]; % intialization value (unit: seconds)
%             p_lower_cplx = [0; 0; 0; 0.003; 0.025; 0.025; 
%                 fs(x,y,z)-75; fs(x,y,z)-25; fs(x,y,z)-25; -pi]; % lower bound
%             p_upper_cplx = [2; 2; 2; 0.025; 0.08; 0.08; 
%                 fs(x,y,z)+75; fs(x,y,z)+25; fs(x,y,z)+25; pi]; % upper bound
% 
%             % create objective function (complex model) and fit
%             objfcn = @(p)objfun_complex_model(p, echo_time, decay);
%             [p_est_cplx, resnorm_cplx] = lsqnonlin(objfcn, p0_cplx, p_lower_cplx, p_upper_cplx, opts);
%             fitted_param_cplx(x,y,z,:) = reshape([p_est_cplx; resnorm_cplx], [1,1,11]);
            
            
            % create magnitude fitting parameters adaptively for each voxel
            p0_mag = [0.2; 0.5; 0.3; 0.01; 0.064; 0.048] % intialization value (unit: seconds)
            p_lower_mag = [0; 0; 0; 0.003; 0.025; 0.025]; % lower bound
            p_upper_mag = [0.5; 1; 1; 0.025; 2; 2]; % upper bound
          
            % create objective function (magnitude model) and fit
            objfcn = @(p)objfun_magnitude_model(p, echo_time, abs(decay));
            [p_est_mag, resnorm_mag] = lsqnonlin(objfcn, p0_mag, p_lower_mag, p_upper_mag, opts);
            fitted_param_mag(x,y,z,:) = reshape([p_est_mag; resnorm_mag], [1,1,7]);
        end
    end
    waitbar(x/x_dim, f, sprintf('Progress: %d %%', floor(x/x_dim*100)));
end
toc
close(f);


% extract MWF values
fitted_mwf_cplx = fitted_param_cplx(:,:,:,1)./(fitted_param_cplx(:,:,:,1)+fitted_param_cplx(:,:,:,2)+fitted_param_cplx(:,:,:,3));
fitted_mwf_mag = fitted_param_mag(:,:,:,1)./(fitted_param_mag(:,:,:,1)+fitted_param_mag(:,:,:,2)+fitted_param_mag(:,:,:,3));

% plot MWF maps for experimental data
figure; suptitle('MWF maps fitted by the complex model (top) and the magnitude model(bottom)');
subplot(231);imshow(squeeze(fitted_mwf_cplx(70,:,:)),[0,0.4]); colorbar;
subplot(232);imshow(squeeze(fitted_mwf_cplx(:,50,:)),[0,0.4]); colorbar;
subplot(233);imshow(squeeze(fitted_mwf_cplx(:,:,60)),[0,0.4]); colorbar;
subplot(234);imshow(squeeze(fitted_mwf_mag(70,:,:)),[0,0.4]); colorbar;
subplot(235);imshow(squeeze(fitted_mwf_mag(:,50,:)),[0,0.4]); colorbar;
subplot(236);imshow(squeeze(fitted_mwf_mag(:,:,60)),[0,0.4]); colorbar;

