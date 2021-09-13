function [orig_cum_dm_ts,orig_dm_ts_sigma,orig_dm_cum_ts_sigma,ts_n,frac_observed,cum_dz_stack_ref] = dz_to_runoff(cum_dz_stack,ts,gridsize,area_mask,basin_mask,basin_id,mode_mask,gmask,err_mask,...
    scaling_method,z,advection_correction_switch,h_adv,rho_runoff,summer_start,summer_end,plot_level)
% function: takes m x n dz stack and obtains dm in gigatonnes based on assumption change occurs at a user defined density within specified region
% and scales based on unobserved area
%
% inputs:
% cum_dz_stack - residual elevation change time series
% area_mask - mask of defined area, note area to keep must == 1 (if don't want to apply, define mask = [])
% zone = 1,2,3 for ice zone
% plot_level - user specified level of plotting
% ts - time vector for dz
% gridsize - resolution (in m) of gridded dataset
% basin_mask - area mask for dz time series (if don't want to apply, define mask = [])
% basin_id - basin number to be masked if masking to basin
% mode_mask - CS2 operating mode mask, mode to keep == 1 (if don't want to apply, define mode_mask = [])
% gmask - ice sheet mask
% err_mask - region to accumulate error over
% scaling_method - switch for scaling method to apply:
% 'total_unobserved_area' = scale according to fraction of observed cells in ablation zone per epoch
% 'elevation_bands' = scale according to fraction of observed cells in elevation bands per epoch
% z - dem (at same posting as cum_dz_stack) for scaling with elevation bands
% advection_correction_switch - switch to apply advection correction - 'y' for on, 'n' for off
% h_adv - map of summer fraction of long term smb, converted to height, for dz correction
% rho_runoff - density used to convert to mass (in kg/m^3)
% summer_start, summer_end - season definitions
% plot_level - user specified level of plotting
%
% outputs:
% orig_cum_dm_ts - time series of dm for user specified region scaled according to unobserved area
% orig_dm_ts_sigma - per epoch error for time series
% orig_dm_cum_ts_sigma - cumulative error for time series
% ts_n - new time vector with data only
% frac_observed - fraction of input region observed per epoch

% initialise outputs
dz_basin_mean = nan(size(ts));
dz_basin_sigma = nan(size(ts));
frac_observed = nan(size(ts));
ts_n = ts;
rmask = ones(size(area_mask)); % master mask for computing observed area

% apply masks and compute master mask
disp('applying masks...')
mask_cum_dz_stack = nan(size(cum_dz_stack)); % initialise
for i = 1:length(ts);
    tmp = cum_dz_stack(:,:,i);
    % mask to mode
    if isempty(mode_mask) == 0
        %disp('applying mode mask...')
        tmp(~mode_mask) = NaN;
        rmask(~mode_mask) = 0;
    end
    % mask to basin
    if isempty(basin_mask) == 0
        %disp('applying basin mask...')
        tmp(basin_mask~=basin_id) = NaN;
        rmask(basin_mask~=basin_id) = 0;
    end
    % mask zone
    if isempty(area_mask) == 0;
        %disp('applying area mask...')
        tmp(~area_mask) = NaN;
        rmask(~area_mask) = 0;
    end
    mask_cum_dz_stack(:,:,i) = tmp;
    rmask = logical(rmask);

    % check if plotting level specified
    if plot_level == 1 & mod(i,10) == 0;
        figure; h = imagesc(tmp); set(h,'alphadata',~isnan(tmp)); axis xy;
    end

    clear tmp
end

% smooth first layer to suppress noise and then fill in missing data to use for referencing
disp('referencing elevation stack to zero...')
cum_dz_t0_layer_sm = nanmedfilt2(cum_dz_stack(:,:,1),25);
filled_cum_dz_t0_layer_sm = inpaint_nans(cum_dz_t0_layer_sm,1);

% mask ice sheet
filled_cum_dz_t0_layer = filled_cum_dz_t0_layer_sm.*(gmask~=0);
filled_cum_dz_t0_layer(filled_cum_dz_t0_layer==0) = NaN;

% reference each grid cell time series so it initiates from zero
cum_dz_stack_ref = nan(size(cum_dz_stack));
for i = 1:size(cum_dz_stack,1)
    for j = 1:size(cum_dz_stack,2)

        % if a time series has data in the first layer then reference to this
        if ~isnan(cum_dz_stack(i,j,1))

            cum_dz_stack_ref(i,j,:) = cum_dz_stack(i,j,:) - cum_dz_stack(i,j,1);

            % if first layer empty then use filled layer and replace empty ice sheet first layer value
        else

            cum_dz_stack_ref(i,j,:) = cum_dz_stack(i,j,:) - filled_cum_dz_t0_layer(i,j);

            if gmask(i,j) == 1

                cum_dz_stack_ref(i,j,1) = 0;

            end

        end

    end
end

% apply firn air content correction
switch advection_correction_switch
    case 'y'
        disp('applying mean smb correction...')
        years = [min(floor(ts)):max(floor(ts))];
        cum_dz_stack_ref_h_adv = nan(size(cum_dz_stack_ref));
        % remove winter firn air content from summer dz
        for i = 1:size(cum_dz_stack_ref,1)
            for j = 1:size(cum_dz_stack_ref,2)
                dz_tmp = cum_dz_stack_ref(i,j,:); dz_tmp = dz_tmp(:);
                if sum(~isnan(dz_tmp)) > 0 & gmask(i,j) == 1 % in continent and data with time series
                    % get summer time series in year
                    for k = 1:length(years)
                        id = find(ts >= summer_start(k) & ts <= summer_end(k));

                        if ~isnan(dz_tmp(3)) % check if data
                            dz_tmp(id(3):end) = dz_tmp(id(3):end)+(h_adv(i,j)); % apply stepwise to summer time series
                            dz_tmp(id(2):end) = dz_tmp(id(2):end)+(h_adv(i,j)); % apply stepwise to summer time series
                        elseif ~isnan(dz_tmp(1))
                            dz_tmp(id(1):end) = dz_tmp(id(1):end)+(h_adv(i,j)); % apply stepwise to summer time series
                        end

                        clearvars id
                    end
                    cum_dz_stack_ref_h_adv(i,j,:) = dz_tmp;
                else
                    continue
                end
                clearvars dz_tmp
            end
        end
        cum_dz_stack_ref = cum_dz_stack_ref_h_adv;
        clearvars cum_dz_stack_ref_h_adv
end

dz_stack_ref = nan(size(cum_dz_stack));

disp('computing volume change...')

% for each x,y,t point in stack compute elevation change since previous measurement and mean density of elevation over the same period
for i = 1:size(cum_dz_stack,1)
    for j = 1:size(cum_dz_stack,2)

        for k = 1:length(ts)

            % set first dz layer to zero
            if k == 1 && gmask(i,j) == 1

                dz_stack_ref(i,j,k) = 0;

            else

                if gmask(i,j)~= 1
                    continue % skip if outside ice sheet

                    % if current and preceeding timestep both have data and ice sheet pixel then compute elevation difference over the time step
                elseif ~isnan(cum_dz_stack_ref(i,j,k)) && ~isnan(cum_dz_stack_ref(i,j,k-1)) && gmask(i,j) == 1

                    % compute elevation difference over the single time step
                    dz_stack_ref(i,j,k) = cum_dz_stack_ref(i,j,k) - cum_dz_stack_ref(i,j,k-1);

                    % if current timestep has data but previous timestep does not have data
                elseif ~isnan(cum_dz_stack_ref(i,j,k)) && isnan(cum_dz_stack_ref(i,j,k-1)) && gmask(i,j) == 1

                    % identify the most recent preceeding entry with data
                    first_previous_t_ind_temp = find(~isnan(cum_dz_stack_ref(i,j,1:k-1)),1,'last');

                    % compute elevation change between the two timesteps
                    dz_stack_ref(i,j,k) = cum_dz_stack_ref(i,j,k) - cum_dz_stack_ref(i,j,first_previous_t_ind_temp);

                    clear first_previous_t_ind_temp

                    % if current timestep does not have data
                elseif isnan(cum_dz_stack_ref(i,j,k))

                    % set timestep dz as empty
                    dz_stack_ref(i,j,k) = nan;

                end
            end
        end
    end
end

% convert to mass
disp('converting to mass...')

% convert elevation change at each timestep from elevation to volume
dv_stack = (dz_stack_ref/1000)*(gridsize/1000)*(gridsize/1000);

% convert volume change at each timestep to mass
dm_stack = dv_stack.*(rho_runoff/1000);

% preallocate mass change vector
orig_dm_ts = nan(1,length(ts)-1);
orig_dm_ts_sigma = nan(1,length(ts)-1);

switch scaling_method
    case 'total_unobserved_area'
        disp('scaling according to unobserved area...')

        % extract mass change for each timestep in desired region
        for i = 1:length(ts)

            orig_layer_temp = dm_stack(:,:,i);

            % mask to melt zone
            orig_layer_temp(~rmask) = NaN;
            orig_dm_ts(i) = sum(orig_layer_temp(~isnan(orig_layer_temp)));

            % compute standard error of dm pixel values in current layer
            orig_dm_ts_sigma(i) = sum(err_mask(:))*(std(orig_layer_temp(~isnan(orig_layer_temp)))/sqrt(sum(~isnan(orig_layer_temp(:)))));

            clear orig_layer_temp

        end

        % get unobserved percentage in specidied region per epoch in order to scale mass time series
        for i = 1:length(ts)
            tmp = mask_cum_dz_stack(:,:,i);
            n_cells = sum(~isnan(tmp(:)));
            frac_observed(i) = (n_cells/sum(rmask(:))); % compute percentage observed
            clear tmp
        end

        % scale according to unobserved area
        orig_dm_ts = orig_dm_ts./(frac_observed*-1); % invert sign so positive = mass loss
        orig_dm_ts_sigma = orig_dm_ts_sigma./frac_observed;

    case 'elevation_bands'
        disp('scaling in elevation bands...')

        % define elevation bands
        elev_bands = [500:500:4000];
        frac_observed = nan(length(ts),length(elev_bands)); % initialise vector to store scaling values
        orig_dm_ts_eb = nan(length(ts),length(elev_bands)); % intialise vector to store mass change per elevation band

        % extract mass change for each timestep in desired region
        for i = 1:length(ts)

            orig_layer_temp = dm_stack(:,:,i);

            % mask to melt zone
            orig_layer_temp(~rmask) = NaN;

            % scale according to unobserved area in elevation band
            for j = 1:length(elev_bands)
                idx = z >= elev_bands(j) - 500 & z < elev_bands(j); idx(~rmask) = 0 ;
                n_cells = sum(~isnan(orig_layer_temp(idx)));
                frac_observed(i,j) = (n_cells/sum(idx(:))); % compute percentage observed
                orig_layer_temp(idx) = (orig_layer_temp(idx))./(frac_observed(i,j));

                clearvars idx n_cells
            end

            orig_dm_ts(i) = sum(orig_layer_temp(~isnan(orig_layer_temp)))*-1; % invert sign so positive = mass loss

            % compute standard error of dm pixel values in current layer
            orig_dm_ts_sigma(i) = sum(err_mask(:))*(std(orig_layer_temp(~isnan(orig_layer_temp)))/sqrt(sum(~isnan(orig_layer_temp(:)))));

            clear orig_layer_temp
        end
end

% clean up time series
id = isnan(orig_dm_ts);
ts_n(id) = []; orig_dm_ts(id) = []; orig_dm_ts_sigma(id) = [];

% accumulate over time
orig_cum_dm_ts = cumsum(orig_dm_ts);
orig_dm_cum_ts_sigma = nan(size(orig_dm_ts_sigma));
% accumulate error
for i = 1:length(ts_n)
    orig_dm_cum_ts_sigma(i) = rssq(orig_dm_ts_sigma(1:i));
end
end % function
