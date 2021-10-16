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
fitted_param = zeros(x_dim, y_dim, z_dim, 11);

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
            if mask(x,y,z) == 0
                continue;
            end
            decay = squeeze(data(x,y,z,:)); % squeeze the dimension
            %decay = squeeze(phantom_denoised(x,y,:));
            decay = decay / abs(decay(1));
            
            % create fitting parameters adaptively for each voxel
            p0 = [0.2; 0.5; 0.3; 0.01; R2s(x,y,z)+10; R2s(x,y,z)-10; 
                fs(x,y,z); fs(x,y,z); fs(x,y,z); 0]; % intialization value (unit: seconds)
            p_lower = [0; 0; 0; 0.003; 0.025; 0.025; 
                fs(x,y,z)-75; fs(x,y,z)-25; fs(x,y,z)-25; -pi]; % lower bound
            p_upper = [2; 2; 2; 0.025; 0.08; 0.08; 
                fs(x,y,z)+75; fs(x,y,z)+25; fs(x,y,z)+25; pi]; % upper bound

            % create object function and fit
            objfcn = @(p)objfun_complex_model(p, echo_time, decay);
            [p_est, resnorm] = lsqnonlin(objfcn, p0, p_lower, p_upper, opts);
            fitted_param(x,y,z,:) = reshape([p_est; resnorm], [1,1,11]);
        end
    end
    waitbar(x/x_dim, f, sprintf('Progress: %d %%', floor(x/x_dim*100)));
end
toc
close(f);


% extract MWF values
fitted_mwf = fitted_param(:,:,:,1)./(fitted_param(:,:,:,1)+fitted_param(:,:,:,2)+fitted_param(:,:,:,3));

% plot MWF maps for experimental data
figure, subplot(131);imshow(squeeze(fitted_mwf(70,:,:)),[0,0.4],'border','tight'); colorbar;
subplot(132);imshow(squeeze(fitted_mwf(:,50,:)),[0,0.4],'border','tight'); colorbar;
subplot(133);imshow(squeeze(fitted_mwf(:,:,60)),[0,0.4],'border','tight'); colorbar;
