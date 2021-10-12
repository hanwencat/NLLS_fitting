% QSM analysis pipeline (main program)
% Created by Hanwen liu on 2021-09-01

%%%%%%%%%%%%%%%%%%%%%% NAMING CONVENTION %%%%%%%%%%%%%%%%%%%%%%%%%
% high dimensional variable
%
% iField - 4 dimensional complex MRI dataset. 
%          the 4th dimension is echo 
%         
% 3D variables
%
% Mask - binary mask denoting the region of interest
% iMag - magnitude image, square root of squares of all echoes
% iFreq_raw - the raw field map, which may contain wrapping, 
%             unit in rad/echo
% N_std - estimated noise standard deviation on iFreq_raw
% iFreq - the unwrapped field map, aka total field (tfs)
%         unit in rad/echo
% RDF - Relative Difference Field, aka local field (lfs)
%       unit in rad/echo
% R2s - R2* map
% QSM - Quantitative Susceptibility Map, 
%       unit in parts per million, aka ppm
%
% QSM_res - Quantitative Susceptibility Map residual, 
%       unit in parts per million, aka ppm

% vectors
%
% B0_dir - unit vector representing direction of B0 field 
% matrix_size - sizes ([x y z]) of the imaging volume
% voxel_size - size of the voxel
%              unit in mm
% TE - echo time, unit in sec
%
% scalars
%
% delta_TE - echo spacing, unit in sec
% CF - center frequency, unit in Hz
% B0_strength - magnetic field strength, unit in Tesla
%%%%%%%%%%%%%%%%%%%%%% NAMING CONVENTION %%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% NECESSARY FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%

% Set paths
addpath(genpath('/export01/data/Hanwen/matlab_tools'));
path_to_repo = '/export01/data/Hanwen/QSM/QSM_scripts'; % full path to QSM analysis pipeline repo
path_to_raw = '/export01/data/Hanwen/data/mgre_data/WT_F_21/14/'; % full path to raw data directory
% path_to_save = '/export01/data/Hanwen/QSM/QSM_analysis/WT_F_21'; % full path to the folder for saving the analysis results

% Add the root directory and all sub-directories of this repo to the matlab search path, and change directory to the save folder;
addpath(genpath(path_to_repo));
% cd(path_to_save);

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

% denoise the data via mppca
iField_mag_denoise = denoise(abs(iField),[5,5,5]);
iField = iField_mag_denoise.*exp(1j*angle(iField));

% % Gibbs ring removal and denoise via mppca
% iField_mag = permute(abs(iField),[1,2,3,4]); % default: no permute applied [1,2,3,4] 
% iField_mag_unring = iField_mag;
% parfor echo = 1:NumEcho
%     iField_mag_unring(:,:,:,echo) = unring(iField_mag(:,:,:,echo));
% end
% iField_mag_unring = permute(iField_mag_unring,[1,2,3,4]); % permute back to match with iField
% iField_mag_unring_denoise = denoise(iField_mag_unring,[5,5,5]); % denoise via mppca
% iField = iField_mag_unring_denoise.*exp(1j*angle(iField)); % inplace corrected iField
% clearvars iField_mag iField_mag_unring iField_mag_unring_denoise


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

% Fit the frequency offset in each of the voxel using a complex fitting
[iFreq_raw, N_std] = Fit_ppm_complex_TE(iField,TE);
% [iFreq_raw, N_std] = Fit_ppm_complex_TE_HL(iField,iFreq_all_echo_corr,TE);


% Spatial phase unwrapping (region-growing)
% iFreq = unwrapPhase(iMag, iFreq_raw, matrix_size);

 
% Background field removal using Projection onto Dipole Fields (Default)
% RDF = PDF(iFreq, N_std, Mask, matrix_size, voxel_size, B0_dir);
lfs_PDF = PDF(tfs, N_std, Mask, matrix_size, voxel_size, B0_dir);

% Background field removal using other methods: LBV, vSHARP, reSHARP
lfs_LBV = LBV(tfs, Mask, matrix_size, voxel_size, 0.01, -1, 1);
lfs_vSHARP = vSHARP(tfs, Mask, [3,3,3]);
lfs_reSHARP = resharp(tfs,Mask,voxel_size, 0.3);
delete mask_p1.bin; % LBV method generates this file

figure;
subplot(1,4,1);
imshow(squeeze(lfs_PDF(70,:,:)),[-0.1,0.1]);
title('PDF'); colorbar('southoutside');
subplot(1,4,2);
imshow(squeeze(lfs_LBV(70,:,:)),[-0.1,0.1]);
title('LBV'); colorbar('southoutside');
subplot(1,4,3);
imshow(squeeze(lfs_vSHARP(70,:,:)),[-0.1,0.1]);
title('vSHARP'); colorbar('southoutside');
subplot(1,4,4);
imshow(squeeze(lfs_reSHARP(70,:,:)),[-0.1,0.1]);
title('reSHARP'); colorbar('southoutside');
drawnow;

RDF = lfs_PDF; % choose the best background removal, default is PDF.

% Monoexponential fitting for R2* map
R2s = arlo(TE, abs(iField));

% Dipole inversion for QSM map 
%[QSM, QSM_res] = tvdi(RDF, double(Mask), voxel_size, 5e-4, iMag, B0_dir);

% plot the results
%figure, plot_QSM_results
% 
% % Save analysis results
% save QSM_results.mat QSM R2s RDF tfs iFreq_raw iMag N_std Mask QSM_res 



