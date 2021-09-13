%%%%%%%%%%%%%%%%%%%%%%%
% Compute runoff from CryoSat-2 plane fit data in the Greenland Ice Sheet
%%%%%%%%%%%%%%%%%%%%%%%

%% initialise
close all; clear;
% addpaths
addpath(genpath('~/Documents/github/ncomms_cs2_runoff/'));
addpath(genpath('~/Documents/MATLAB/mathworks/'));

%% 1. load data
% load plane fit data
disp('loading plane fit data...')
disp('loading baseline d data...')
load('/Volumes/eartsl/gris_smb/surfacefit_cs2_greenland_5km_lrmsin_cycle1_cycle122_cumul_dz_stack_60_day.mat','cs_cumdz_stack','gridx','gridy','mean_t','ts_midpt_sampling_vec','ts_vec_out') % baseline d
% clear up time vector
t_cs = ts_vec_out; clear ts_vec_out
% mask to 2011 onwards
tr = t_cs<2011;
t_cs(tr) = [];
cs_cumdz_stack(:,:,tr) = [];
tn_cs = t_cs + (30/365.25); % shift time vector to mid point
% define seasons
summer_start = tn_cs(2:6:end);
summer_end = tn_cs(4:6:end);

% apply advection correction - use long term mean SMB as correction for long term dynamics
% load long term mean SM
load('/Volumes/eartsl/gris_smb/racmo_smb_mean_1960_1979.mat')
% regrid to cs2 grid
smb_mean_ref_cs2 = griddata(rxx,ryy,smb_mean_ref,gridx,gridy,'nearest');
dh_advection_summer_frac = 0.3285.*smb_mean_ref_cs2./917; % convert to height and get summer fraction for correction, assume density of ice for ablation zone

% zwally mask
load('/Volumes/eartsl/gris_smb/zmask_gris.mat','zmask_gris');

% dynamics mask
disp('loading velocity derived dynamics mask...')
load('/Volumes/eartsl/gris_smb/dynamics_masks_from_velocities.mat','vdmask_final')
dmask = vdmask_final; clear vdmask_final

% define ice mask
disp('defining ice mask...')
x = gridx(1,:); y = gridy(:,1);
load('/Volumes/eartsl/gris_smb/continent_outlines.mat','gx','gy');
icx=1:length(x); icy=1:length(y);
icxgl=interp1(x,icx,gx);
icygl=interp1(y,icy,gy);
cmask=poly2mask(icxgl,icygl,length(y),length(x));
% clean up
clearvars icx icy icxgl icygl

% load gimp dem for interpolation scheme
load('/Volumes/eartsl/gris_smb/z_gimp_cs2.mat')
z_gimp_cs2 = flipud(z_gimp_cs2); z_gimp_cs2(~cmask) = NaN; % mask to ice sheet

% load runoff zone definition from racmo
load('/Volumes/eartsl/gris_smb/runoff_pc_area.mat')
runoff_pc_area_cs2 = griddata(rxx,ryy,runoff_pc_area,gridx,gridy,'nearest');
clearvars rxx ryy runoff_pc_area

%% 2 compute runoff using final solution parameters
disp('computing smb anomaly from cs2 data...')
% load dynamics mask
disp('loading velocity derived dynamics mask...')
load('/Volumes/eartsl/gris_smb/dynamics_masks_from_velocities.mat','vdmask_final')
dmask = vdmask_final; clear vdmask_final

%% 2.1 ablation zone
disp('computing runoff in ablation zone...')
% initialise ablation zone mask
load('/Volumes/eartsl/gris_smb/melt_zone_mask1116.mat','melt_zone_mask1116')
% remove dynamic areas
disp('removing dynamic areas...')
ablation_runoff_mask = melt_zone_mask1116==1; ablation_runoff_mask(dmask) = 0;

% compute runoff
[cs_ablation_zone_h_adv,cs_ablation_zone_h_adv_sigma,~,~,frac_observed,~] = dz_to_runoff(cs_cumdz_stack,tn_cs,5000,ablation_runoff_mask==1,[],[],[],cmask,melt_zone_mask1116==1,'elevation_bands',z_gimp_cs2,...
    'y',dh_advection_summer_frac,917,summer_start,summer_end,0);
% sum to obtain runoff per year
cs_ablation_zone_yrs_h_adv = nan(size(summer_start));
cs_ablation_zone_yrs_h_adv_sigma = nan(size(summer_start));

for j = 1:length(summer_start) % match time period to fdm to compare
    id = tn_cs >= summer_start(j) & tn_cs <= summer_end(j);
    yr_runoff = cs_ablation_zone_h_adv(id);
    cs_ablation_zone_yrs_h_adv(j) = yr_runoff(end)-yr_runoff(1);
    cs_ablation_zone_yrs_h_adv_sigma(j) = rssq(cs_ablation_zone_h_adv_sigma(id));
    clearvars id yr_runoff
end

%% 2.2 inland runoff zone
disp('computing runoff in inland zone...')
% initalise inland runoff zone mask
inland_runoff_mask = runoff_pc_area_cs2 <= 100; % define area
inland_runoff_mask = bwmorph(inland_runoff_mask,'clean',inf); % clean up with bwmorph operations

inland_runoff_mask(melt_zone_mask1116==1) = 0; inland_runoff_mask(~cmask) = 0; % remove areas outside ablation zone and ice sheet

% remove dynamic areas
disp('removing dynamic areas...')
inland_runoff_mask_nodyn = inland_runoff_mask;
inland_runoff_mask_nodyn(dmask) = 0;

% compute runoff
[cs_inland_zone_h_adv,cs_inland_zone_h_adv_sigma,~,~,~,~] = dz_to_runoff(cs_cumdz_stack,tn_cs,5000,inland_runoff_mask_nodyn==1,[],[],[],cmask,inland_runoff_mask==1,'elevation_bands',z_gimp_cs2,...
    'y',dh_advection_summer_frac,684,summer_start,summer_end,0);

% % sum to obtain runoff per year
cs_inland_zone_yrs_h_adv = nan(size(summer_start));
cs_inland_zone_yrs_h_adv_sigma = nan(size(summer_start));

for j = 1:length(summer_start) % match time period to fdm to compare
    id = tn_cs >= summer_start(j) & tn_cs <= summer_end(j);
    yr_runoff = cs_inland_zone_h_adv(id);
    d_yr_runoff = diff(yr_runoff); % for inland area only take epochs where runoff > 0
    cs_inland_zone_yrs_h_adv(j) = sum(d_yr_runoff(d_yr_runoff>0));
    cs_inland_zone_yrs_h_adv_sigma(j) = rssq(cs_inland_zone_h_adv_sigma(id));
    clearvars yr_runoff
end

% add for total runoff
cs_runoff_combination_sum = cs_ablation_zone_yrs_h_adv + cs_inland_zone_yrs_h_adv;
cs_runoff_combination_sum_sigma = rssq([cs_ablation_zone_yrs_h_adv_sigma;cs_inland_zone_yrs_h_adv_sigma]);

%% display values
T = table([2011:2020]',round(cs_runoff_combination_sum,0)',round(cs_runoff_combination_sum_sigma,0)');
T.Properties.VariableNames = {'Year','Runoff','Error'};
disp(T)

fprintf('Average runoff = %.0f +/- %.0f Gt/yr\n',mean(cs_runoff_combination_sum),rssq(cs_runoff_combination_sum_sigma)/sqrt(length(cs_runoff_combination_sum_sigma)))
fprintf('Runoff variability = %.0f Gt/yr\n',std(cs_runoff_combination_sum))
fprintf('Runoff spread = %.0f Gt/yr\n',max(cs_runoff_combination_sum)-min(cs_runoff_combination_sum))
