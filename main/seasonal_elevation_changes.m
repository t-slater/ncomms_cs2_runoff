%%%%%%%%%%%%%%%%%%%%%%%
% Produce seasonal elevation change time series from CryoSat-2 plane fit data in the Greenland Ice Sheet Ablation Zone
% Compare to equivalent time series derived from firn modelling
%%%%%%%%%%%%%%%%%%%%%%%

%% initialise
close all; clear;
% addpaths
addpath(genpath('~/Documents/github/ncomms_cs2_runoff/'));
addpath(genpath('~/Documents/MATLAB/mathworks/'));

% define plot level
plot_level = 0;

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

% load imau fdm
load('/Volumes/eartsl/gris_smb/racmo_v2p2_zs_stack_v3_60day.mat','racmo_zs_anomaly_mthly_cs2grid', 't_mthly')
fdm_cumdz_stack = racmo_zs_anomaly_mthly_cs2grid; clear racmo_zs_anomaly_mthly_cs2grid
t_fdm = t_mthly; clear t_mthly
% mask to 2011 onwards
tr = t_fdm<2011;
t_fdm(tr) = [];
fdm_cumdz_stack(:,:,tr) = [];

clear tr

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

% load basin mask
load('/Volumes/eartsl/gris_smb/zmask_gris.mat','zmask_gris');

%  mask racmo to cryosat observed areas
fdm_cumdz_stack_m = fdm_cumdz_stack;
for i = 1:length(t_fdm);
    tmp = fdm_cumdz_stack(:,:,i);
    k = isnan(cs_cumdz_stack(:,:,i)); % mask to cs2 observed
    tmp(k) = NaN;
    fdm_cumdz_stack_m(:,:,i) = tmp;
    clearvars tmp k
end

%% 2. compute average time series
disp('computing average time series...')
%% 2.1.1 ablation zone
disp('extracting mean ablation zone time series...')

% load melt zone mask
load('/Volumes/eartsl/gris_smb/melt_zone_mask1116.mat','melt_zone_mask1116')
melt_zone_mask1116
% load dynamics mask
disp('loading velocity derived dynamics mask...')
load('/Volumes/eartsl/gris_smb/dynamics_masks_from_velocities.mat','vdmask_final')
dmask = vdmask_final; clear vdmask_final

% compute time series over the ablation zone
[cs_abl_mean_dz,cs_abl_mean_dz_sigma,cs_abl_mean_dz_sigma_cum,tn_cs] = dz_mean(cs_cumdz_stack,t_cs,melt_zone_mask1116==1,[],[],[],0); % cryosat
[fdm_abl_mean_dz,~,~,tn_fdm] = dz_mean(fdm_cumdz_stack_m,t_fdm,melt_zone_mask1116==1,[],[],[],0); % imau fdm

% clean up cryosat time series
id = isnan(cs_abl_mean_dz);
t_cs(id) = []; cs_abl_mean_dz(id) = []; cs_abl_mean_dz_sigma(id) = []; cs_abl_mean_dz_sigma_cum(id) = [];

% smooth with gaussain to preserve seasonality
cs_abl_mean_dz_sm = smoothdata(cs_abl_mean_dz,'gaussian',3);
fdm_abl_mean_dz_sm = smoothdata(fdm_abl_mean_dz,'gaussian',3);

%% 2.1.2 zwally basin means
% larger zwally basins
basin_id = [1,2,3,4,5,6,7,8]; % indices for basins

zwallybasins2_area = nan(size(basin_id));
% get basin areas
for i = 1:length(basin_id)
    tmp = floor(zmask_gris) == basin_id(i);
    tmp(melt_zone_mask1116~=1) = 0;
    zwallybasins2_area(i) = sum(tmp(:))*5*5;
    clear tmp
end

% initialise arrays to store smoothed and filtered time series
cs_zwallybasins2_mean_dz = nan(length(tn_cs),length(basin_id)); cs_zwallybasins2_sigma_dz = nan(length(tn_cs),length(basin_id)); cs_zwallybasins2_sigma_cum_dz = nan(length(tn_cs),length(basin_id));
fdm_zwallybasins2_mean_dz = nan(length(tn_fdm),length(basin_id)); fdm_zwallybasins2_sigma_dz = nan(length(tn_fdm),length(basin_id)); fdm_zwallybasins2_sigma_cum_dz = nan(length(tn_fdm),length(basin_id));

% loop through each basin
for i = 1:length(basin_id);
    [cs_basin_mean,cs_basin_sigma,cs_basin_sigma_cum,tnb_cs] = dz_mean(cs_cumdz_stack,t_cs,melt_zone_mask1116==1,floor(zmask_gris),basin_id(i),[],0);
    [fdm_basin_mean,fdm_basin_sigma,fdm_basin_sigma_cum,tnb_fdm] = dz_mean(fdm_cumdz_stack_m,t_fdm,melt_zone_mask1116==1,floor(zmask_gris),basin_id(i),[],0);

    % smooth with gaussian to preserve seasonality
    cs_dz_tmp = smoothdata(cs_basin_mean,'gaussian',3);
    fdm_dz_tmp = smoothdata(fdm_basin_mean,'gaussian',3);

    % store
    cs_zwallybasins2_mean_dz(:,i) = cs_dz_tmp - cs_dz_tmp(1); cs_zwallybasins2_sigma_dz(:,i) = cs_basin_sigma ; cs_zwallybasins2_sigma_cum_dz(:,i) = cs_basin_sigma_cum;
    fdm_zwallybasins2_mean_dz(:,i) = fdm_dz_tmp - fdm_dz_tmp(1); fdm_zwallybasins2_sigma_dz(:,i) = fdm_basin_sigma ; fdm_zwallybasins2_sigma_cum_dz(:,i) = fdm_basin_sigma_cum;
    clearvars cs_basin_mean cs_basin_sigma cs_basin_sigma_cum cs_dz_tmp fdm_basin_mean fdm_basin_sigma fdm_basin_sigma_cum fdm_dz_tmp
end

%% 3. seasonal changes in elevation
%% 3.1 derive seasonal changes
% Define summer as CryoSat-2 data points inclusive of May 1st - August 31st
summer_start = tn_cs(2:6:end);
summer_end = tn_cs(4:6:end);
% 3.1.1 ablation zone
disp('computing seasonal elevation changes across ablation zone')
% initalise
% summer
cs_abl_summer_dz = nan(1,10); cs_abl_summer_dz_err = nan(1,10);
fdm_abl_summer_dz = nan(1,6); fdm_abl_summer_dz_err = nan(1,6);
% winter
cs_abl_winter_dz = nan(1,9); cs_abl_winter_dz_err = nan(1,9);
fdm_abl_winter_dz = nan(1,5); fdm_abl_winter_dz_err = nan(1,5);

% peak to peak method
% get basin time series and error
cs_dz_tmp = cs_abl_mean_dz; cs_dz_tmp_err = cs_abl_mean_dz_sigma_cum;
fdm_dz_tmp = fdm_abl_mean_dz; fdm_dz_tmp_err = fdm_abl_mean_dz_err;
% find max and min in summer period and compute seasonal elevation changes
if plot_level == 1
    figure; hold on; plot(tn_cs,cs_dz_tmp,'-k'); plot(tn_fdm,fdm_dz_tmp,'--k'); % check
end
% summer
% cryosat
for j = 1:10
    id = tn_cs >= summer_start(j) & tn_cs <= summer_end(j);
    cs_dz_summer_max = max(cs_dz_tmp(id)); cs_dz_summer_min = min(cs_dz_tmp(id));
    cs_dz_summer_max_err = cs_dz_tmp_err(cs_dz_tmp == cs_dz_summer_max); cs_dz_summer_min_err = cs_dz_tmp_err(cs_dz_tmp == cs_dz_summer_min); % find corresponding error
    cs_abl_summer_dz(j) = cs_dz_summer_min - cs_dz_summer_max;
    cs_abl_summer_dz_err(j) = rssq([cs_dz_summer_min_err,cs_dz_summer_max_err]);

    if plot_level == 1
        plot(tn_cs(id),cs_dz_tmp(id),'rs');
    end
    clearvars id cs_dz_summer_max cs_dz_summer_min cs_dz_summer_max_err cs_dz_summer_min_err
end

% fdm
for j = 1:6
    id = tn_fdm >= summer_start(j) & tn_fdm <= summer_end(j);
    fdm_dz_summer_max = max(fdm_dz_tmp(id)); fdm_dz_summer_min = min(fdm_dz_tmp(id));
    fdm_dz_summer_max_err = fdm_dz_tmp_err(fdm_dz_tmp == fdm_dz_summer_max); fdm_dz_summer_min_err = fdm_dz_tmp_err(fdm_dz_tmp == fdm_dz_summer_min); % find corresponding error
    fdm_abl_summer_dz(j) = fdm_dz_summer_min - fdm_dz_summer_max;
    fdm_abl_summer_dz_err(j) = rssq([fdm_dz_summer_min_err,fdm_dz_summer_max_err]);

    if plot_level == 1
        plot(tn_fdm(id),fdm_dz_tmp(id),'rs');
    end
    clearvars id fdm_dz_summer_max fdm_dz_summer_min fdm_dz_summer_max_err fdm_dz_summer_min_err
end

% winter
if plot_level == 1
    figure; hold on; plot(tn_cs,cs_dz_tmp,'-k'); plot(tn_fdm,fdm_dz_tmp,'--k'); % check
end
% cryosat
for j = 1:9
    id = tn_cs >= summer_end(j) & tn_cs < (summer_start(j)-(60/365.25))+1;
    cs_dz_winter_max = max(cs_dz_tmp(id)); cs_dz_winter_min = min(cs_dz_tmp(id));
    cs_dz_winter_max_err = cs_dz_tmp_err(cs_dz_tmp == cs_dz_winter_max); cs_dz_winter_min_err = cs_dz_tmp_err(cs_dz_tmp == cs_dz_winter_min); % find corresponding error
    cs_abl_winter_dz(j) = cs_dz_winter_max - cs_dz_winter_min;
    cs_abl_winter_dz_err(j) = rssq([cs_dz_winter_max_err,cs_dz_winter_min_err]);

    if plot_level == 1
        plot(tn_cs(id),cs_dz_tmp(id),'bo');
    end
    clearvars id cs_dz_winter_max cs_dz_winter_min cs_dz_winter_max_err cs_dz_winter_min_err
end
% fdm
for j = 1:5
    id = tn_fdm >= summer_end(j) & tn_fdm < (summer_start(j)-(60/365.25))+1;
    fdm_dz_winter_max = max(fdm_dz_tmp(id)); fdm_dz_winter_min = min(fdm_dz_tmp(id));
    fdm_dz_winter_max_err = fdm_dz_tmp_err(fdm_dz_tmp == fdm_dz_winter_max); fdm_dz_winter_min_err = fdm_dz_tmp_err(fdm_dz_tmp == fdm_dz_winter_min); % find corresponding error
    fdm_abl_winter_dz(j) = fdm_dz_winter_max - fdm_dz_winter_min;
    fdm_abl_winter_dz_err(j) = rssq([fdm_dz_winter_max_err,fdm_dz_winter_min_err]);

    if plot_level == 1
        plot(tn_fdm(id),fdm_dz_tmp(id),'bo');
    end
    clearvars id fdm_dz_winter_max fdm_dz_winter_min fdm_dz_winter_max_err fdm_dz_winter_min_err
end



% 3.1.2 large zwally basins
disp('computing seasonal elevation changes for principal zwally basins...')

% initalise
% summer
cs_zwallybasins2_summer_dz = nan(size(cs_zwallybasins2_mean_dz,2),10); cs_zwallybasins2_summer_dz_err = nan(size(cs_zwallybasins2_mean_dz,2),10);
fdm_zwallybasins2_summer_dz = nan(size(fdm_zwallybasins2_mean_dz,2),6); fdm_zwallybasins2_summer_dz_err = nan(size(fdm_zwallybasins2_mean_dz,2),6);
% winter
cs_zwallybasins2_winter_dz = nan(size(cs_zwallybasins2_mean_dz,2),9); cs_zwallybasins2_winter_dz_err = nan(size(cs_zwallybasins2_mean_dz,2),9);
fdm_zwallybasins2_winter_dz = nan(size(fdm_zwallybasins2_mean_dz,2),5); fdm_zwallybasins2_winter_dz_err = nan(size(fdm_zwallybasins2_mean_dz,2),5);

% cycle through each basin
for i = 1:size(cs_zwallybasins2_mean_dz,2)
    % peak to peak method
    % get basin time series and error
    cs_dz_tmp = cs_zwallybasins2_mean_dz(:,i); cs_dz_tmp_err = cs_zwallybasins2_sigma_dz(:,i);
    fdm_dz_tmp = fdm_zwallybasins2_mean_dz(:,i); fdm_dz_tmp_err = fdm_zwallybasins2_sigma_dz(:,i);
    % find max and min in summer period and compute seasonal elevation changes
    if plot_level == 1
        figure; hold on; plot(tn_cs,cs_dz_tmp,'-k'); plot(tn_fdm,fdm_dz_tmp,'--k'); % check
    end
    % summer
    % cryosat
    for j = 1:10
        id = tn_cs >= summer_start(j) & tn_cs <= summer_end(j);
        cs_dz_summer_max = max(cs_dz_tmp(id)); cs_dz_summer_min = min(cs_dz_tmp(id));
        cs_dz_summer_max_err = cs_dz_tmp_err(cs_dz_tmp == cs_dz_summer_max); cs_dz_summer_min_err = cs_dz_tmp_err(cs_dz_tmp == cs_dz_summer_min); % find corresponding error
        cs_zwallybasins2_summer_dz(i,j) = cs_dz_summer_min - cs_dz_summer_max;
        cs_zwallybasins2_summer_dz_err(i,j) = rssq([cs_dz_summer_min_err,cs_dz_summer_max_err]);

        if plot_level == 1
            plot(tn_cs(id),cs_dz_tmp(id),'rs');
        end
        clearvars id cs_dz_summer_max cs_dz_summer_min cs_dz_summer_max_err cs_dz_summer_min_err
    end
    % fdm
    for j = 1:6
        id = tn_fdm >= summer_start(j) & tn_fdm <= summer_end(j);
        fdm_dz_summer_max = max(fdm_dz_tmp(id)); fdm_dz_summer_min = min(fdm_dz_tmp(id));
        fdm_dz_summer_max_err = fdm_dz_tmp_err(fdm_dz_tmp == fdm_dz_summer_max); fdm_dz_summer_min_err = fdm_dz_tmp_err(fdm_dz_tmp == fdm_dz_summer_min); % find corresponding error
        fdm_zwallybasins2_summer_dz(i,j) = fdm_dz_summer_min - fdm_dz_summer_max;
        fdm_zwallybasins2_summer_dz_err(i,j) = rssq([fdm_dz_summer_min_err,fdm_dz_summer_max_err]);

        if plot_level == 1
            plot(tn_fdm(id),fdm_dz_tmp(id),'rs');
        end
        clearvars id fdm_dz_summer_max fdm_dz_summer_min fdm_dz_summer_max_err fdm_dz_summer_min_err
    end

    % winter
    % cryosat
    for j = 1:9
        id = tn_cs >= summer_end(j) & tn_cs < (summer_start(j)-(60/365.25))+1;
        cs_dz_winter_max = max(cs_dz_tmp(id)); cs_dz_winter_min = min(cs_dz_tmp(id));
        cs_dz_winter_max_err = cs_dz_tmp_err(cs_dz_tmp == cs_dz_winter_max); cs_dz_winter_min_err = cs_dz_tmp_err(cs_dz_tmp == cs_dz_winter_min); % find corresponding error
        cs_zwallybasins2_winter_dz(i,j) = cs_dz_winter_max - cs_dz_winter_min;
        cs_zwallybasins2_winter_dz_err(i,j) = rssq([cs_dz_winter_max_err,cs_dz_winter_min_err]);

        if plot_level == 1
            plot(tn_cs(id),cs_dz_tmp(id),'bo');
        end
        clearvars id cs_dz_winter_max cs_dz_winter_min cs_dz_winter_max_err cs_dz_winter_min_err
    end
    % fdm
    for j = 1:5
        id = tn_fdm >= summer_end(j) & tn_fdm < (summer_start(j)-(60/365.25))+1;
        fdm_dz_winter_max = max(fdm_dz_tmp(id)); fdm_dz_winter_min = min(fdm_dz_tmp(id));
        fdm_dz_winter_max_err = fdm_dz_tmp_err(fdm_dz_tmp == fdm_dz_winter_max); fdm_dz_winter_min_err = fdm_dz_tmp_err(fdm_dz_tmp == fdm_dz_winter_min); % find corresponding error
        fdm_zwallybasins2_winter_dz(i,j) = fdm_dz_winter_max - fdm_dz_winter_min;
        fdm_zwallybasins2_winter_dz_err(i,j) = rssq([fdm_dz_winter_max_err,fdm_dz_winter_min_err]);

        if plot_level == 1
            plot(tn_fdm(id),fdm_dz_tmp(id),'bo');
        end
        clearvars id fdm_dz_winter_max fdm_dz_winter_min fdm_dz_winter_max_err fdm_dz_winter_min_err
    end

end

% 3.2 fit piecewise linear trends to examine spatial variability in seasonal elevation changes in 20 x 20 km grid cells
gridspacing = 20e3; % define coarser grid
plot_level = 0;

% regrid data

if gridspacing > 5e3
    xn = [min(gridx(1,:)):gridspacing:max(gridx(1,:))]; yn = [min(gridy(:,1)):gridspacing:max(gridy(:,1))];
    [gridxn,gridyn] = meshgrid(xn,yn); % define new grid
    clearvars xn yn
    % regrid melt zone mask
    melt_zone_mask1116n = griddata(gridx,gridy,melt_zone_mask1116,gridxn,gridyn,'nearest');
    % regrid ice sheet mask
    cmaskn = logical(griddata(gridx,gridy,double(cmask),gridxn,gridyn,'nearest'));

    % initalise new dz stack
    cs_cumdz_stackn = nan(size(gridxn,1),size(gridxn,2),size(cs_cumdz_stack,3));
    fdm_cumdz_stackn = nan(size(gridxn,1),size(gridxn,2),size(fdm_cumdz_stack,3));

    % average cum_dz_stack onto new resolution
    for i = 1:size(gridxn,1) % cycle through each cell
        fprintf('regridding %.0f percent complete...\n',(i/length(gridxn(1,:)))*100)
        for j = 1:size(gridyn,2)
            %if melt_zone_mask1116n(i,j) == 1 % only if in desired region
            if cmaskn(i,j) == 1 % only if in desired region
                % make mask for region
                mask1 = zeros(size(gridxn)); mask1(i,j) = 1;
                % regrid to parent resolution
                mask2 = logical(griddata(gridxn,gridyn,mask1,gridx,gridy,'nearest'));
                clear mask1
                % average cs and fdm time series and add to new stacks
                % cryosat
                for k = 1:size(cs_cumdz_stack,3)
                    cs_tmp = cs_cumdz_stack(:,:,k); % get stack at timestamp
                    cs_cumdz_stackn(i,j,k) = nanmean(cs_tmp(mask2));
                    clear cs_tmp
                end
                % fdm
                for k = 1:size(fdm_cumdz_stack,3)
                    fdm_tmp = fdm_cumdz_stack_m(:,:,k);
                    fdm_cumdz_stackn(i,j,k) = nanmean(fdm_tmp(mask2));
                    clear fdm_tmp
                end
                clearvars mask2
            else
                continue
            end
        end
    end
else
    gridxn = gridx; gridyn = gridy;
    melt_zone_mask1116n = melt_zone_mask1116;
    cmaskn = cmask;
    cs_cumdz_stackn = cs_cumdz_stack; fdm_cumdz_stackn = fdm_cumdz_stack;
end

% fit piecewise linear trends to examine seasonal elevation changes on new grid
% initialise
cs_grid_summer_combined_dzdt = nan(size(gridxn));
cs_grid_summer_combined_dzdt_err = nan(size(gridxn));
cs_grid_winter_combined_dzdt = nan(size(gridxn));
cs_grid_winter_combined_dzdt_err = nan(size(gridxn));

fdm_grid_summer_combined_dzdt = nan(size(gridxn));
fdm_grid_summer_combined_dzdt_err = nan(size(gridxn));
fdm_grid_winter_combined_dzdt = nan(size(gridxn));
fdm_grid_winter_combined_dzdt_err = nan(size(gridxn));

for i = 1:size(gridxn,1) % cycle through each cell
    fprintf('gridding seasonal trends %.0f percent complete...\n',(i/length(gridxn(1,:)))*100)
    for j = 1:size(gridyn,2)% form tempory mask
        if cmaskn(i,j) == 1 % only if in desired region
            % get time series in grid cell
            cs_area_mean = cs_cumdz_stackn(i,j,:); cs_area_mean = cs_area_mean(:)';
            fdm_area_mean = fdm_cumdz_stackn(i,j,:); fdm_area_mean = fdm_area_mean(:)';

            tna_cs = t_cs;
            tna_fdm = t_fdm;

            id = isnan(cs_area_mean); % clean up cryosat time series
            tna_cs(id) = []; cs_area_mean(id) = [];
            clear id

            if ~isempty(cs_area_mean) & ~isempty(fdm_area_mean) % check for empty grid cell
                % smooth time series with gaussain to reduce noise and preserve seasonality
                cs_dz_tmp = smoothdata(cs_area_mean,'gaussian',3);
                fdm_dz_tmp = smoothdata(fdm_area_mean,'gaussian',3);
                tna_cs = tna_cs + (30/365.25); % adjust so time series corresponds to midpoint of time period
                tna_fdm = tna_fdm + (30/365.25);
                cs_dz_tmp = (cs_dz_tmp - cs_dz_tmp(1))'; cs_dz_tmp = cs_dz_tmp - 10; % arbitrary offset to fix bug in centre of grav crossing zero
                fdm_dz_tmp = (fdm_dz_tmp - fdm_dz_tmp(1))'; fdm_dz_tmp = fdm_dz_tmp - 7;

                % compute average seasonal trend

                % summer
                % cryosat
                cs_dz_to_fit = [];
                tn_cs_to_fit = [];
                for k = 1:10

                    % align all data with first year using centre of gravity
                    if k == 1
                        cs_dz_summer = cs_dz_tmp(id);
                        tn_cs_summer = tna_cs(id)';

                        % get centre of gravity
                        cogx1 = sum(cs_dz_tmp(id).*tna_cs(id)')/sum(cs_dz_tmp(id));
                        cogy1 = sum(cs_dz_tmp(id).^2)/sum(cs_dz_tmp(id));

                        cs_dz_to_fit = cat(1,cs_dz_to_fit,cs_dz_summer);
                        tn_cs_to_fit = cat(1,tn_cs_to_fit,tn_cs_summer);

                        if plot_level == 1 && mod(i,25) == 1
                            figure; hold on
                            plot(tna_cs,cs_dz_tmp,'-k')
                            plot(tn_cs_summer,cs_dz_summer,'x')
                        end

                        clearvars cs_dz_summer tn_cs_summer
                    else
                        if variable_seasons
                            id = tna_cs >= (floor(summer_start(k))+melt_start(i,j)) & tna_cs <= (floor(summer_start(k))+melt_end(i,j));
                        else
                            id = tna_cs >= summer_start(k) & tna_cs <= summer_end(k);
                        end
                        cs_dz_summer = cs_dz_tmp(id);
                        tn_cs_summer = tna_cs(id)';

                        % get centre of gravity
                        cogx2 = sum(cs_dz_tmp(id).*tna_cs(id)')/sum(cs_dz_tmp(id));
                        cogy2 = sum(cs_dz_tmp(id).^2)/sum(cs_dz_tmp(id));

                        % realign and add to array
                        t_offset = cogx2 - cogx1;
                        dz_offset = cogy2 - cogy1;
                        cs_dz_to_fit = cat(1,cs_dz_to_fit,cs_dz_summer-dz_offset);
                        tn_cs_to_fit = cat(1,tn_cs_to_fit,tn_cs_summer-t_offset);

                        if plot_level == 1 && mod(i,25) == 1
                            plot(tn_cs_summer-t_offset,cs_dz_summer-dz_offset,'x')
                        end
                        clearvars cs_dz_summer tn_cs_summer cogx2 cogy2 t_offset dz_offset
                    end

                end

                if ~isempty(cs_dz_to_fit)
                    % least squares fit to get seasonal trend
                    f = fitlm(tn_cs_to_fit,cs_dz_to_fit,'linear');
                    cs_grid_summer_combined_dzdt(i,j) = f.Coefficients.Estimate(2);
                    cs_grid_summer_combined_dzdt_err(i,j) = f.Coefficients.SE(2);
                else
                    cs_grid_summer_combined_dzdt(i,j) = NaN;
                    cs_grid_summer_combined_dzdt_err(i,j) = NaN;
                end

                clearvars tn_cs_to_fit cs_dz_to_fit f

                % fdm
                fdm_dz_to_fit = [];
                tn_fdm_to_fit = [];
                for k = 1:6

                    if variable_seasons
                        id = tn_fdm >= (floor(summer_start(k))+melt_start(i,j)) & tn_fdm <= (floor(summer_start(k))+melt_end(i,j));
                    else
                        id = tn_fdm >= summer_start(k) & tn_fdm <= summer_end(k);
                    end
                    % align all data with first year using centre of gravity
                    if k == 1
                        fdm_dz_summer = fdm_dz_tmp(id);
                        tn_fdm_summer = tn_fdm(id)';

                        % get centre of gravity
                        cogx1 = sum(fdm_dz_tmp(id).*tn_fdm(id)')/sum(fdm_dz_tmp(id));
                        cogy1 = sum(fdm_dz_tmp(id).^2)/sum(fdm_dz_tmp(id));

                        fdm_dz_to_fit = cat(1,fdm_dz_to_fit,fdm_dz_summer);
                        tn_fdm_to_fit = cat(1,tn_fdm_to_fit,tn_fdm_summer);

                        if plot_level == 1 && mod(i,25) == 1
                            plot(tna_fdm,fdm_dz_tmp,'--k')
                            plot(tn_fdm_summer,fdm_dz_summer,'o')
                        end

                        clearvars fdm_dz_summer tn_fdm_summer
                    else
                        if variable_seasons
                            id = tn_fdm >= (floor(summer_start(k))+melt_start(i,j)) & tn_fdm <= (floor(summer_start(k))+melt_end(i,j));
                        else
                            id = tn_fdm >= summer_start(k) & tn_fdm <= summer_end(k);
                        end

                        fdm_dz_summer = fdm_dz_tmp(id);
                        tn_fdm_summer = tn_fdm(id)';

                        % get centre of gravity
                        cogx2 = sum(fdm_dz_tmp(id).*tn_fdm(id)')/sum(fdm_dz_tmp(id));
                        cogy2 = sum(fdm_dz_tmp(id).^2)/sum(fdm_dz_tmp(id));

                        % realign and add to array
                        t_offset = cogx2 - cogx1;
                        dz_offset = cogy2 - cogy1;
                        fdm_dz_to_fit = cat(1,fdm_dz_to_fit,fdm_dz_summer-dz_offset);
                        tn_fdm_to_fit = cat(1,tn_fdm_to_fit,tn_fdm_summer-t_offset);

                        if plot_level == 1 && mod(i,25) == 1
                            plot(tn_fdm_summer-t_offset,fdm_dz_summer-dz_offset,'o')
                        end
                        clearvars fdm_dz_summer tn_fdm_summer cogx2 cogy2 t_offset dz_offset
                    end

                end

                % least squares fit to get seasonal trend
                if ~isempty(fdm_dz_to_fit)
                    f = fitlm(tn_fdm_to_fit,fdm_dz_to_fit,'linear');
                    fdm_grid_summer_combined_dzdt(i,j) = f.Coefficients.Estimate(2);
                    fdm_grid_summer_combined_dzdt_err(i,j) = f.Coefficients.SE(2);
                else
                    fdm_grid_summer_combined_dzdt(i,j) = NaN;
                    fdm_grid_summer_combined_dzdt_err(i,j) = NaN;
                end
                clearvars tn_fdm_to_fit fdm_dz_to_fit f


                % winter
                % cryosat
                cs_dz_to_fit = [];
                tn_cs_to_fit = [];
                for k = 1:9
                    if variable_seasons
                        id = tna_cs >= (floor(summer_start(k))+melt_end(i,j)) & tna_cs < (floor(summer_start(k))+melt_start(i,j)-(60/365.25))+1;
                    else
                        id = tna_cs >= summer_end(k) & tna_cs < (summer_start(k)-(60/365.25))+1;
                    end

                    % align all data with first year using centre of gravity
                    if k == 1
                        cs_dz_winter = cs_dz_tmp(id);
                        tn_cs_winter = tna_cs(id)';

                        % get centre of gravity
                        cogx1 = sum(cs_dz_tmp(id).*tna_cs(id)')/sum(cs_dz_tmp(id));
                        cogy1 = sum(cs_dz_tmp(id).^2)/sum(cs_dz_tmp(id));

                        cs_dz_to_fit = cat(1,cs_dz_to_fit,cs_dz_winter);
                        tn_cs_to_fit = cat(1,tn_cs_to_fit,tn_cs_winter);

                        if plot_level == 1 && mod(i,25) == 1
                            figure; hold on
                            plot(tna_cs,cs_dz_tmp,'-k')
                            plot(tn_cs_winter,cs_dz_winter,'x')
                        end

                        clearvars cs_dz_winter tn_cs_winter
                    else
                        if variable_seasons
                            id = tna_cs >= (floor(summer_start(k))+melt_end(i,j)) & tna_cs < (floor(summer_start(k))+melt_start(i,j)-(60/365.25))+1;
                        else
                            id = tna_cs >= summer_end(k) & tna_cs < (summer_start(k)-(60/365.25))+1;
                        end

                        cs_dz_winter = cs_dz_tmp(id);
                        tn_cs_winter = tna_cs(id)';

                        % get centre of gravity
                        cogx2 = sum(cs_dz_tmp(id).*tna_cs(id)')/sum(cs_dz_tmp(id));
                        cogy2 = sum(cs_dz_tmp(id).^2)/sum(cs_dz_tmp(id));

                        % realign and add to array
                        t_offset = cogx2 - cogx1;
                        dz_offset = cogy2 - cogy1;
                        cs_dz_to_fit = cat(1,cs_dz_to_fit,cs_dz_winter-dz_offset);
                        tn_cs_to_fit = cat(1,tn_cs_to_fit,tn_cs_winter-t_offset);

                        if plot_level == 1 && mod(i,25) == 1
                            plot(tn_cs_winter-t_offset,cs_dz_winter-dz_offset,'x')
                        end
                        clearvars cs_dz_winter tn_cs_winter cogx2 cogy2 t_offset dz_offset
                    end

                end

                if ~isempty(cs_dz_to_fit)
                    % least squares fit to get seasonal trend
                    f = fitlm(tn_cs_to_fit,cs_dz_to_fit,'linear');
                    cs_grid_winter_combined_dzdt(i,j) = f.Coefficients.Estimate(2);
                    cs_grid_winter_combined_dzdt_err(i,j) = f.Coefficients.SE(2);
                else
                    cs_grid_winter_combined_dzdt(i,j) = NaN;
                    cs_grid_winter_combined_dzdt_err(i,j) = NaN;
                end

                clearvars tn_cs_to_fit cs_dz_to_fit f

                % fdm
                fdm_dz_to_fit = [];
                tn_fdm_to_fit = [];
                for k = 1:6
                    if variable_seasons
                        id = tn_fdm >= (floor(summer_start(k))+melt_end(i,j)) & tn_fdm < (floor(summer_start(k))+melt_start(i,j)-(60/365.25))+1;
                    else
                        id = tn_fdm >= summer_end(k) & tn_fdm < (summer_start(k)-(60/365.25))+1;
                    end
                    % align all data with first year using centre of gravity
                    if k == 1
                        fdm_dz_winter = fdm_dz_tmp(id);
                        tn_fdm_winter = tn_fdm(id)';

                        % get centre of gravity
                        cogx1 = sum(fdm_dz_tmp(id).*tn_fdm(id)')/sum(fdm_dz_tmp(id));
                        cogy1 = sum(fdm_dz_tmp(id).^2)/sum(fdm_dz_tmp(id));

                        fdm_dz_to_fit = cat(1,fdm_dz_to_fit,fdm_dz_winter);
                        tn_fdm_to_fit = cat(1,tn_fdm_to_fit,tn_fdm_winter);

                        if plot_level == 1 && mod(i,25) == 1
                            plot(tna_fdm,fdm_dz_tmp,'--k')
                            plot(tn_fdm_winter,fdm_dz_winter,'o')
                        end

                        clearvars fdm_dz_winter tn_fdm_winter
                    else
                        if variable_seasons
                            id = tn_fdm >= (floor(summer_start(k))+melt_end(i,j)) & tn_fdm < (floor(summer_start(k))+melt_start(i,j)-(60/365.25))+1;
                        else
                            id = tn_fdm >= summer_end(k) & tn_fdm < (summer_start(k)-(60/365.25))+1;
                        end
                        fdm_dz_winter = fdm_dz_tmp(id);
                        tn_fdm_winter = tn_fdm(id)';

                        % get centre of gravity
                        cogx2 = sum(fdm_dz_tmp(id).*tn_fdm(id)')/sum(fdm_dz_tmp(id));
                        cogy2 = sum(fdm_dz_tmp(id).^2)/sum(fdm_dz_tmp(id));

                        % realign and add to array
                        t_offset = cogx2 - cogx1;
                        dz_offset = cogy2 - cogy1;
                        fdm_dz_to_fit = cat(1,fdm_dz_to_fit,fdm_dz_winter-dz_offset);
                        tn_fdm_to_fit = cat(1,tn_fdm_to_fit,tn_fdm_winter-t_offset);

                        if plot_level == 1 && mod(i,25) == 1
                            plot(tn_fdm_winter-t_offset,fdm_dz_winter-dz_offset,'o')
                        end
                        clearvars fdm_dz_winter tn_fdm_winter cogx2 cogy2 t_offset dz_offset
                    end

                end

                % least squares fit to get seasonal trend
                if ~isempty(fdm_dz_to_fit)
                    f = fitlm(tn_fdm_to_fit,fdm_dz_to_fit,'linear');
                    fdm_grid_winter_combined_dzdt(i,j) = f.Coefficients.Estimate(2);
                    fdm_grid_winter_combined_dzdt_err(i,j) = f.Coefficients.SE(2);
                else
                    fdm_grid_winter_combined_dzdt(i,j) = NaN;
                    fdm_grid_winter_combined_dzdt_err(i,j) = NaN;
                end
                clearvars tn_fdm_to_fit fdm_dz_to_fit f
            end
        else
            continue
        end

        clearvars cs_area_mean tna_cs cs_dz_tmp
        clearvars fdm_area_mean tna_fdm fdm_dz_tmp

    end
end
