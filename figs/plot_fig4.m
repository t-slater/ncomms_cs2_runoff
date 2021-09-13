%% initialise
close all; clear;
% addpaths
addpath('~/Documents/MATLAB/mathworks/');
addpath('~/Documents/MATLAB/mathworks/arctic_mapping_tools/');
addpath('~/Documents/MATLAB/mathworks/climate_data_toolbox/');
addpath('~/Documents/MATLAB/mathworks/bedmachine/');
addpath('~/Documents/MATLAB/mathworks/climate_data_toolbox/cdt_data/');
addpath('~/Documents/MATLAB/mathworks/cptcmap/');

% load data
% imau fdm May - August elevation change by basin from 1980
A = readmatrix('/Users/thomas/OneDrive - University of Leeds/manuscripts/gris_runoff/accept_in_principle/data_files/slater_et_al_2021_ncomms_fig4_ablation_zone_dh_may_aug_imau_fdm.csv','Range','B3:AF10');
imau_fdm_may_aug_dh_zwallybasins = A; clear A

% cryosat May - August elevation change by basin
A = readmatrix('/Users/thomas/OneDrive - University of Leeds/manuscripts/gris_runoff/accept_in_principle/data_files/slater_et_al_2021_ncomms_fig3_panelc_ablation_zone_dh_may_aug_cryosat.csv','Range','B3:K10');
cryosat_may_aug_dh_zwallybasins = A; clear A

% ice sheet runoff
A = readmatrix('/Users/thomas/OneDrive - University of Leeds/manuscripts/gris_runoff/accept_in_principle/data_files/slater_et_al_2021_ncomms_fig4_gris_runoff_cryosat.csv','Range','B2:C11');
cryosat_runoff = A(:,1);
cryosat_runoff_err = A(:,2);
clear A

A = readmatrix('/Users/thomas/OneDrive - University of Leeds/manuscripts/gris_runoff/accept_in_principle/data_files/slater_et_al_2021_ncomms_fig4_gris_runoff_racmo.csv','Range','A2:B64');
racmo_runoff = A; clear A

A = readmatrix('/Users/thomas/OneDrive - University of Leeds/manuscripts/gris_runoff/accept_in_principle/data_files/slater_et_al_2021_ncomms_fig4_gris_runoff_mar.csv','Range','A2:B62');
mar_runoff = A; clear A


% aux data
    load('/Volumes/eartsl/gris_smb/cs_grn_dhdt_fill.mat','gridx','gridy','gmask') % ice sheet mask
    load('/Volumes/eartsl/gris_smb/zmask_gris.mat','zmask_gris'); % zwally basin definitions
    load('/Volumes/eartsl/gris_smb/continent_outlines.mat','gx','gy'); % ice sheet outline
    load('/Volumes/eartsl/gris_smb/GIS_forTom_7dec16.mat','melt_zone_mask1114'); % ablation zone area

 basin_id = [1,2,3,4,5,6,7,8]; % indices for basins
% make seasonal dh grids for mapping
% fdm
% may - august
imau_fdm_may_aug_dh_zwallybasins_grid = nan(size(gridx,1),size(gridx,2),size(imau_fdm_may_aug_dh_zwallybasins,2));
for i = 1:size(imau_fdm_may_aug_dh_zwallybasins,2)
    tmp = ones(size(zmask_gris));
    for j = 1:length(basin_id)
        tmp(floor(zmask_gris)==basin_id(j)) = tmp(floor(zmask_gris)==basin_id(j))*imau_fdm_may_aug_dh_zwallybasins(j,i); % get seasonal dzdt for basin and year
    end
    tmp(melt_zone_mask1114>1) = NaN; tmp(~gmask) = NaN; % mask to ablation zone, continent and observed basins
    imau_fdm_may_aug_dh_zwallybasins_grid(:,:,i) = tmp; % assign
    clearvars tmp
end

% cryosat
% may - august
cryosat_may_aug_dh_zwallybasins_grid = nan(size(gridx,1),size(gridx,2),size(cryosat_may_aug_dh_zwallybasins,2));
for i = 1:size(cryosat_may_aug_dh_zwallybasins,2)
    tmp = ones(size(zmask_gris));
    for j = 1:length(basin_id)
        tmp(floor(zmask_gris)==basin_id(j)) = tmp(floor(zmask_gris)==basin_id(j))*cryosat_may_aug_dh_zwallybasins(j,i); % get seasonal dz for basin and year
    end
    tmp(melt_zone_mask1114>1) = NaN; tmp(~gmask) = NaN;  % mask to ablation zone, continent and observed basins
    cryosat_may_aug_dh_zwallybasins_grid(:,:,i) = tmp; % assign
    clearvars tmp
end

dz_basins_grid = cat(3,imau_fdm_may_aug_dh_zwallybasins_grid,cryosat_may_aug_dh_zwallybasins_grid); % concatenate together

%% make plot
years = [2011:2020];
fig = figure('units','centimeters','position',[0 0 180 73]/10);
delete(gca)
ha = tight_subplot(3,1,.01,[.12 .12],[.1 .1]);
delete(ha(1))

hapos2 = get(ha(2),'position'); hapos3 = get(ha(3),'position');
ax = axes('position',[hapos3(1),hapos3(2),hapos2(3),hapos3(4)*2.25]); delete(ha(2)); delete(ha(3)); hold on
axis([1980 2020 0 650])
yticks([0:100:600])
% indicate decades
vline(1990,'-','color',rgb('gray'))
vline(2000,'-','color',rgb('gray'))
vline(2010,'-','color',rgb('gray'))

% cryosat runoff shaded region
[l,p] = boundedline(years,cryosat_runoff,cryosat_runoff_err,'color',rgb('ocean blue'),'alpha','linewidth',1);
ho = outlinebounds(l,p);
set(ho,'linewidth',.5)

% racmo runoff
id = racmo_runoff(:,1) >= 1980;
plot([racmo_runoff(id,1);2020],[racmo_runoff(id,2);373.98],'.-','markersize',9,'color',rgb('orange'),'linewidth',1);

% mar runoff
id = mar_runoff(:,1) >= 1980;
plot([mar_runoff(id,1);2020],[mar_runoff(id,2);368.65],'.-','markersize',9,'color',rgb('reddish'),'linewidth',1);

% cryosat runoff (overplotted for visibility)
errorbar(years(1),cryosat_runoff(1),cryosat_runoff_err(1),'color',rgb('ocean blue'),'linewidth',.5,'capsize',0);
errorbar(years(end),cryosat_runoff(end),cryosat_runoff_err(end),'color',rgb('ocean blue'),'linewidth',.5,'capsize',0);
plot(years,cryosat_runoff,'.-','markersize',9,'color',rgb('ocean blue'),'linewidth',1);

xlabel('Year'); ylabel('Runoff (Gt/yr)')
axis([1980 2020 0 650])
ntitle({'   CryoSat-2'},'location','nw','color',rgb('ocean blue'),'fontsize',8)
ntitle({'','    RACMO2.3p2'},'location','nw','color',rgb('orange'),'fontsize',8)
ntitle({'','','    MARv3.11'},'location','nw','color',rgb('reddish'),'fontsize',8)
box on

set(gca,'clipping','off','fontsize',8)
plot([1980,1980,2010,2010],([490,465,465,490]-5)+200,'k')
plot([2011,2011,2020,2020],([490,465,465,490]-5)+200,'k')
text(1994,685,'IMAU-FDM','fontsize',6)
text(2014,685,'CryoSat-2','fontsize',6)

% seasonal elevation change maps
vert_stagger = repmat([.01:.01:.05]*2.5,1,(2020-1980)/5); vert_stagger = [vert_stagger,.01*2.5]; % stagger vertically

for i = 1:length([1980:2020])
  %axb = axes('position',[hapos3(1)-.01+((i-1)*0.02),0.7,.02 .3]); hold on;
  axb = axes('position',[hapos3(1)-.01+((i-1)*0.02),0.625-vert_stagger(i),.02 .6]); hold on;

  greenland('patch','facecolor',rgb('gray'),'edgecolor','none'); patch(gx,gy,'w','edgecolor','none');
  imagescn(gridx(1,:),gridy(:,1),dz_basins_grid(:,:,i));
  %[con,h] = contour(gridx(1,:),gridy(:,1),floor(zmask_gris),[1:8],'k'); h.LineWidth = .05;
  cmocean('-matter'); caxis([-1.5,0]);
  xlim([-.65e6 .25e6])
  ylim([-3.4e6 -.65e6])
  set(axb,'visible','off','clipping','off')

end

% colorbar
axc = axes('position',[hapos3(1)-.075,0.725,.175 .2]); hold on;
colormap(cmocean('-matter')); caxis([-1.5,0]);
c = colorbar('westoutside','fontsize',6);
title(c,'(m)','fontsize',6);
%c.Label.Position(1) = 0.3;
%c.Label.Position(2) = 1.5;
set(axc,'visible','off');
cbarrow('down')
