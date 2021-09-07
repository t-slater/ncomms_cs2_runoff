function [dz_basin_mean,dz_basin_sigma,dz_basin_sigma_cum,ts_n] = dz_mean(cum_dz_stack,ts,area_mask,basin_mask,basin_id,mode_mask,plot_level)
  % function: takes m x n dz stack and obtains mean and error dz time series for user specified mask
  %
  % inputs:
  % cum_dz_stack - residual elevation change time series
  % ts - time vector for dz
  % area_mask - mask of defined area, note area to keep must == 1 (if don't want to apply, define mask = [])
  % basin_mask - area mask for dz time series
  % basin_id - basin number to be masked if masking to basin
  % mode_mask - CS2 operating mode mask, mode to keep == 1 (if don't want to apply, define mode_mask = [])
  % plot_level - user specified level of plotting

  % outputs:
  % basin_mean_dz - mean time series of dz for user specified basin
  % dz_basin_sigma - per epoch error for time series
  % dz_basin_sigma_cum - accumulated error for time series
  % ts_n - new time vector with data only
% initialise outputs
dz_basin_mean = nan(size(ts));
dz_basin_sigma = nan(size(ts));
ts_n = ts;

% apply masks
mask_cum_dz_stack = nan(size(cum_dz_stack)); % initialise
for i = 1:length(ts);
    tmp = cum_dz_stack(:,:,i);
    % mask to mode
    if isempty(mode_mask) == 0
        %disp('applying mode mask...')
        tmp(~mode_mask) = NaN;
    end
    % mask to basin
    if isempty(basin_mask) == 0
      %disp('applying basin mask...')
      tmp(basin_mask~=basin_id) = NaN;
    end
    % mask zone
    if isempty(area_mask) == 0;
        %disp('applying area mask...')
        tmp(~area_mask) = NaN;
    end
    mask_cum_dz_stack(:,:,i) = tmp;

    % check if plotting level specified
    if plot_level == 1 & mod(i,10) == 0;
        figure; h = imagesc(tmp); set(h,'alphadata',~isnan(tmp)); axis xy;
    end

    clear tmp
end

% mean for zone and error
for i = 1:length(ts)
    tmp = mask_cum_dz_stack(:,:,i);
    dz_basin_mean(i) = nanmean(tmp(:));
    n_cells = sum(~isnan(tmp(:)));
    dz_basin_sigma(i) = nanstd(tmp(:))/sqrt(n_cells);
    clear tmp
end

% clean up time series
id = isnan(dz_basin_mean);
ts_n(id) = []; dz_basin_mean(id) = []; dz_basin_sigma(id) = [];

dz_basin_sigma_cum = nan(size(ts_n));
% accumulate error in time
% accumulate error
    for i = 1:length(ts_n)
        dz_basin_sigma_cum(i) = rssq(dz_basin_sigma(1:i));
    end

end % function
