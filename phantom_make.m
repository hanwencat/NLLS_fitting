function [signal, mwf, noise, fs_my, fs_ax, fs_ex, phi0] = phantom_make(snr_range, mwf_range, t2s, x_dim, y_dim, fs_mu, fs_sigma, echo_time)
    


    % create placeholders
    snr = zeros(x_dim, y_dim);
    mwf = zeros(x_dim, y_dim);
    
    num_echo = size(echo_time,1);
    signal = complex(zeros(x_dim, y_dim, num_echo));
    noise = complex(zeros(x_dim, y_dim, num_echo));
    
    % produce the ground truth snr map
    for y = 1:y_dim
        snr(:,y) = snr_range(1) + y*(snr_range(2)-snr_range(1))/y_dim; % not include snr lower boundary
    end
    
    % produce the ground truth mwf map
    for x = 1:x_dim
        mwf(x,:) = mwf_range(1) + x*(mwf_range(2)-mwf_range(1))/x_dim; % not include mwf lower boundary
    end
    
    % produce the random fs maps for my-, ax-, and ex- water pools    
    fs_my_mu = fs_mu(1);
    fs_ax_mu = fs_mu(2);
    fs_ex_mu = fs_mu(3);
    fs_my_sigma = fs_sigma(1);
    fs_ax_sigma = fs_sigma(2);
    fs_ex_sigma = fs_sigma(3);
    
    fs_my = normrnd(fs_my_mu, fs_my_sigma, [x_dim, y_dim]);
    fs_ax = normrnd(fs_ax_mu, fs_ax_sigma, [x_dim, y_dim]);
    fs_ex = normrnd(fs_ex_mu, fs_ex_sigma, [x_dim, y_dim]);
    
    % produce the random initial phase map
    phi0 = -pi + 2*pi*rand(x_dim, y_dim);
    
    % produce the decay signals (no T1 compensation)
    t2_my = t2s(1);
    t2_ax = t2s(2);
    t2_ex = t2s(3);
    
    for x = 1:x_dim
        for y = 1:y_dim
            % generate a random ratio for ax/(ax+ex) 
            ax_ratio = rand(1);
            %ax_ratio = 1; % 2-pool model
            
            % produce complex noise for each echo according to SNR (independent Gaussian noise on real and imaginary axes) 
            noise_mu = 0; % noise mean is 0
            noise_sigma = 1/(snr(x,y)*((pi/2)^0.5)); % noise variance is calculated according to snr
            noise(x,y,:) = (normrnd(noise_mu,noise_sigma,[num_echo,1])+ 1j*normrnd(noise_mu,noise_sigma,[num_echo,1]));
            
            % generate signal with noise added
            signal(x,y,:) = (mwf(x,y) * exp(-(1/t2_my + 1j*2*pi*fs_my(x,y)).*echo_time)...
                + (1 - mwf(x,y)) * ax_ratio * exp(-(1/t2_ax + 1j*2*pi*fs_ax(x,y)).*echo_time)...
                + (1 - mwf(x,y)) * (1-ax_ratio) * exp(-(1/t2_ex + 1j*2*pi*fs_ex(x,y)).*echo_time)) * exp(-1j*phi0(x,y))...
                + reshape(noise(x,y,:), [num_echo,1]);
        end
    end
     
    
end
        