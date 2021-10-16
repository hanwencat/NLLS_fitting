% Preprocessing for NLLS complex fitting
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


% Set paths
addpath(genpath('/export01/data/Hanwen/matlab_tools'));
path_to_repo = '/export01/data/Hanwen/QSM/QSM_scripts'; % full path to QSM analysis pipeline repo
path_to_raw = '/export01/data/Hanwen/data/mgre_data/WT_F_21/14'; % full path to raw data directory

% Add the root directory and all sub-directories of this repo to the matlab search path, and change directory to the save folder;
addpath(genpath(path_to_repo));

% Get parameters from raw Bruker directory:
[voxel_size, matrix_size, TE, delta_TE, CF, Affine3D, B0_dir, TR, NumEcho, Read_offset, Phase1_offset, Phase2_offset, Slice_offset] = Read_Bruker_raw_hanwen(path_to_raw);

% Save parameters
% save parameters.mat voxel_size matrix_size TE delta_TE CF Affine3D B0_dir TR NumEcho Read_offset Phase1_offset Phase2_offset Slice_offset;

% Reconstruct the image with customized bipolar phase correction 
[iField] = recon_mgre_raw(append(path_to_raw,'/fid'), 256, matrix_size, voxel_size, NumEcho, Read_offset, Phase1_offset, Slice_offset);

figure;
imshow(squeeze(angle(iField(:,50,:))),[-pi,pi]);
title('Raw phase after bipolar correction');
colorbar;
drawnow;

% % denoise the data via mppca
% iField_mag_denoise = denoise(abs(iField),[5,5,5]);
% iField = iField_mag_denoise.*exp(1j*angle(iField));

% Gibbs ring removal and denoise via mppca
iField_mag = permute(abs(iField),[3,1,2,4]); % permute for better unring performance
iField_mag_unring = iField_mag;
parfor echo = 1:NumEcho
    iField_mag_unring(:,:,:,echo) = unring(iField_mag(:,:,:,echo));
end
iField_mag_unring = permute(iField_mag_unring,[2,3,1,4]); % permute back to match with original dimension

iField_mag_unring_denoise = denoise(iField_mag_unring,[5,5,5]); % denoise via mppca
iField = iField_mag_unring_denoise.*exp(1j*angle(iField)); % inplace corrected iField
clearvars iField_mag iField_mag_unring iField_mag_unring_denoise


% Compute magnitude image combining all echoes
iMag = sqrt(sum(abs(iField).^2, 4));

% Use FSL BET to extract brain mask
Mask = BET(iMag, matrix_size, voxel_size, 0.5); % fractional intensity threshold is chosen to be the deault 0.5

% Unwrap individual echoes 
iFreq_all_echo = unwrap_all_echo(iField, iMag, Mask);

figure;
imshow(squeeze(iFreq_all_echo(:,50,:)),[-pi,pi]);
title('Unwrapped individual echoes');
colorbar;
drawnow;

% Unwrap between echoes
iFreq_all_echo_corr = unwrap_between_echoes(iFreq_all_echo, Mask);

figure;
imshow(squeeze(iFreq_all_echo_corr(:,50,:)),[-5*pi,5*pi]);
title('Unwrapped between echoes');
colorbar;
drawnow;

% fit echoes to get total field map
% [tfs, res] = echofit(iFreq_all_echo_corr, abs(iField), TE*1000);
[tfs, res] = echofit(iFreq_all_echo_corr(:,:,:,1:2:end), abs(iField(:,:,:,1:2:end)), TE(1:2:end)*1000); % use half echoes

% Monoexponential fitting for R2* map
R2s = arlo(TE, abs(iField));
