%% NLLS complex fitting for experimental data

% load data, echo time, and frequency shift map
data = iField.*Mask;
% mask = Mask(69:70,:,:);
% data = data(69:70,:,:,:);
echo_time = TE;
% data = iField(:,:,:,1:2:end).*Mask;
% echo_time = TE(1:2:end);
fs = -tfs*1000/(2*pi); % unit in Hz

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
f = waitbar(0, 'Starting');
tic
%signal = phantom_denoised .* exp(1j*angle(phantom.signal));
for x = 1:x_dim
    parfor y = 1:y_dim
        for z = 1:z_dim
            if Mask(x,y,z) == 0
                continue;
            end
            decay = squeeze(data(x,y,z,:)); % squeeze the dimension
            %decay = squeeze(phantom_denoised(x,y,:));
            decay = decay / abs(decay(1));
            
            % create fitting parameters
            p0 = [0.2; 0.5; 0.3; 0.01; 0.064; 0.048; 
                fs(x,y,z); fs(x,y,z); fs(x,y,z); 0]; % intialization value (unit: seconds)
            p_lower = [0; 0; 0; 0.003; 0.025; 0.025; 
                fs(x,y,z)-75; fs(x,y,z)-25; fs(x,y,z)-25; -pi]; % lower bound
            p_upper = [2; 2; 2; 0.025; 0.08; 0.08; 
                fs(x,y,z)+75; fs(x,y,z)+25; fs(x,y,z)+25; pi]; % upper bound

            % create object function and fit
            objfcn = @(p)objfun_complex_model_lsqnonlin(p, echo_time, decay);
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
