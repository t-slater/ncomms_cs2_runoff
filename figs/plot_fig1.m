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
filename = '/Users/thomas/OneDrive - University of Leeds/manuscripts/gris_runoff/accept_in_principle/data_files/slater_et_al_2021_ncomms_fig1.nc';
dhdt_cryosat_2011_2020 = ncread(filename,'dhdt_cryosat_2011_2020')';
dhdt_cryosat_may_aug_avg = ncread(filename,'dhdt_cryosat_may_aug_avg')';
dhdt_cryosat_sep_apr_avg = ncread(filename,'dhdt_cryosat_sep_apr_avg')';

% aux data
load('/Volumes/eartsl/gris_smb/cs_grn_dhdt_fill.mat','gridx','gridy','gmask') % ice sheet mask
load('/Volumes/eartsl/gris_smb/melt_zone_mask1116.mat','melt_zone_mask1116') % facies mask
load('/Volumes/eartsl/gris_smb/dynamics_masks_from_velocities.mat','vdmask_final') % dynamics mask
dmask = vdmask_final; clear vdmask_final

% create colormap
cmap = cptcmap('~/Documents/MATLAB/mathworks/cptcmap/cptfiles/cbcRdBu.cpt');
cmap(1:100,:) = cmap(1+10:100+10,:);
cmap(220-100:220,:) = cmap(220-100-10:220-10,:);
cmap(100:120,:) = 0.96;

% load bedmachine
[bed,x,y] = bedmachine_data('bed','xy');

%% make plot
fig = figure('units','centimeters','position',[0 0 180 100]/10);
delete(gca)
ha = tight_subplot(1,3);

% panel a long term dhdt
hapos = get(ha(1),'position');
axb = axes('position',hapos); hold on;

% plot bedmachine background
sc = .1;
imagesc(imresize(x,sc),imresize(y,sc),imresize(bed,sc))
caxis([-1 0]*10000); cmocean('ice');
set(axb,'visible','off');

ax = axes('position',hapos); hold on; delete(ha(1));
ax.Color = 'none';
set(ax,'xtick',[],'ytick',[]);
greenland('patch','facecolor',rgb('gray'),'linewidth',.2,'edgecolor','none');
imagescn(gridx(1,:),gridy(:,1),dhdt_cryosat_2011_2020);

% dynamics mask contour
contour(gridx(1,:),gridy(:,1),dmask,[1 1],'color',rgb('light purple'),'linewidth',0.75);
colormap(ax,cmap); caxis([-1 1])

linkaxes([axb,ax])
c = colorbar('eastoutside','position',[(hapos(1)+hapos(3))-0.06 hapos(2)+0.09, 0.025 0.18]); title(c,'(m/yr)','fontsize',7)
set(c,'ticks',[-1:1:1])
axis([min(x),max(x),min(y),max(y)])
daspect([1 1 1]); set(ax,'visible','off');
ntitle('2011-2020','location','se','fontsize',8);
ntitle('a','location','nw','fontsize',10,'fontweight','bold');

% panel b May - August average dhdt
hapos = get(ha(2),'position');
axb = axes('position',hapos); hold on;

% plot bedmachine background
sc = .1;
imagesc(imresize(x,sc),imresize(y,sc),imresize(bed,sc))
caxis([-1 0]*10000); cmocean('ice');
set(axb,'visible','off');

ax = axes('position',hapos); hold on; delete(ha(2));
ax.Color = 'none';
set(ax,'xtick',[],'ytick',[]);
greenland('patch','facecolor',rgb('gray'),'linewidth',.2,'edgecolor','none');
imagescn(gridx(1,:),gridy(:,1),dhdt_cryosat_may_aug_avg);

% ablation zone contour
contour(gridx(1,:),gridy(:,1),melt_zone_mask1116==1,[1 1],'color',rgb('light purple'),'linewidth',0.75);
colormap(ax,cmap); caxis([-3 3])

linkaxes([axb,ax])
c = colorbar('eastoutside','position',[(hapos(1)+hapos(3))-0.06 hapos(2)+0.09, 0.025 0.18]); title(c,'(m/yr)','fontsize',7)
set(c,'ticks',[-3:3:3])
axis([min(x),max(x),min(y),max(y)])
daspect([1 1 1]); set(ax,'visible','off');
ntitle('May-August','location','se','fontsize',8);
ntitle('b','location','nw','fontsize',10,'fontweight','bold');

% panel c September - April dhdt
hapos = get(ha(3),'position');
axb = axes('position',hapos); hold on;

% plot bedmachine background
sc = .1;
imagesc(imresize(x,sc),imresize(y,sc),imresize(bed,sc))
caxis([-1 0]*10000); cmocean('ice');
set(axb,'visible','off');

ax = axes('position',hapos); hold on; delete(ha(3));
ax.Color = 'none';
set(ax,'xtick',[],'ytick',[]);
greenland('patch','facecolor',rgb('gray'),'linewidth',.2,'edgecolor','none');
imagescn(gridx(1,:),gridy(:,1),dhdt_cryosat_sep_apr_avg);

colormap(ax,cmap); caxis([-3 3])
linkaxes([axb,ax])
c = colorbar('eastoutside','position',[(hapos(1)+hapos(3))-0.06 hapos(2)+0.09, 0.025 0.18]); title(c,'(m/yr)','fontsize',7)
set(c,'ticks',[-3:3:3])
axis([min(x),max(x),min(y),max(y)])
daspect([1 1 1]); set(ax,'visible','off');
ntitle('September-April','location','se','fontsize',8);
ntitle('c','location','nw','fontsize',10,'fontweight','bold');

set(gcf,'color','w','InvertHardCopy', 'off');

set(fig,'PaperOrientation','landscape');
set(fig,'PaperUnits','centimeters');
set(fig,'PaperSize',[180 100]/10)
print('/Users/thomas/OneDrive - University of Leeds/manuscripts/gris_runoff/accept_in_principle/vector_figs/fig1.pdf','-painters','-dpdf');

close
clearvars fig ha