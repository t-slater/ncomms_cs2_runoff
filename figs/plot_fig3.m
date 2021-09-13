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
% panel a - dh time series
% cryosat
A = readmatrix('/Users/thomas/OneDrive - University of Leeds/manuscripts/gris_runoff/accept_in_principle/data_files/slater_et_al_2021_ncomms_fig2_cryosat.csv','Range',2);
t_cryosat = A(:,1); dh_cryosat = A(:,2); dh_cryosat_err = A(:,3); clear A
% imau fdm
A = readmatrix('/Users/thomas/OneDrive - University of Leeds/manuscripts/gris_runoff/accept_in_principle/data_files/slater_et_al_2021_ncomms_fig2_imau_fdm.csv','Range',2);
t_imau_fdm = A(:,1); dh_imau_fdm = A(:,2); clear A

% panel - b dh icebridge comparison
A = readmatrix('/Users/thomas/OneDrive - University of Leeds/manuscripts/gris_runoff/accept_in_principle/data_files/slater_et_al_2021_ncomms_fig3_panelb_basin6.2.csv','Range','B2:D12');
oib_dh_b62 = A(:,1);
cryosat_dh_b62 = A(:,2);
imau_fdm_dh_b62 = A(1:7,3);
clear A

A = readmatrix('/Users/thomas/OneDrive - University of Leeds/manuscripts/gris_runoff/accept_in_principle/data_files/slater_et_al_2021_ncomms_fig3_panelb_basin7.2.csv','Range','B2:D11');
oib_dh_b72 = A(:,1);
cryosat_dh_b72 = A(:,2);
imau_fdm_dh_b72 = A(1:7,3);
clear A

% panel c - seasonal dh
% cryosat
A = readmatrix('/Users/thomas/OneDrive - University of Leeds/manuscripts/gris_runoff/accept_in_principle/data_files/slater_et_al_2021_ncomms_fig3_panelc_ablation_zone_dh_may_aug_cryosat.csv','Range','B3:K10');
cryosat_may_aug_dh_zwallybasins = A; clear A

A = readmatrix('/Users/thomas/OneDrive - University of Leeds/manuscripts/gris_runoff/accept_in_principle/data_files/slater_et_al_2021_ncomms_fig3_panelc_ablation_zone_dh_sep_apr_cryosat.csv','Range','B3:J10');
cryosat_sep_apr_dh_zwallybasins = A; clear A

% imau fdm
A = readmatrix('/Users/thomas/OneDrive - University of Leeds/manuscripts/gris_runoff/accept_in_principle/data_files/slater_et_al_2021_ncomms_fig3_panelc_ablation_zone_dh_may_aug_imau_fdm.csv','Range','B3:G10');
imau_fdm_may_aug_dh_zwallybasins = A; clear A

A = readmatrix('/Users/thomas/OneDrive - University of Leeds/manuscripts/gris_runoff/accept_in_principle/data_files/slater_et_al_2021_ncomms_fig3_panelc_ablation_zone_dh_sep_apr_imau_fdm.csv','Range','B3:F10');
imau_fdm_sep_apr_dh_zwallybasins = A; clear A

% panel d - seasonal dh icebridge comparison
A = readmatrix('/Users/thomas/OneDrive - University of Leeds/manuscripts/gris_runoff/accept_in_principle/data_files/slater_et_al_2021_ncomms_fig3_paneld_basin6.2.csv','Range','B2:D5');
oib_dhdt_b62 = A(:,1);
cryosat_dhdt_b62 = A(:,2);
imau_fdm_dhdt_b62 = A(1:3,3);
clear A

A = readmatrix('/Users/thomas/OneDrive - University of Leeds/manuscripts/gris_runoff/accept_in_principle/data_files/slater_et_al_2021_ncomms_fig3_paneld_basin7.2.csv','Range','B2:D5');
oib_dhdt_b72 = A(:,1);
cryosat_dhdt_b72 = A(:,2);
imau_fdm_dhdt_b72 = A(1:3,3);
clear A

%% make plot
fig = figure('units','centimeters','position',[0 0 88 88]/10);

% panel a
subplot(2,2,1); hold on
plot([-10 10],[-10 10],'-k');
plot(dh_cryosat(2:35),dh_imau_fdm(2:end),'k.')
daspect([1 1 1]); box on;
axis([-4 0,-4,0])
set(gca,'xtick',[-4:1:0],'ytick',[-4:1:0],'fontsize',6)
xlabel('CryoSat-2 (m)'); ylabel('IMAU-FDM (m)');
ntitle('RMS = 23 cm','location','se','fontsize',6)
ntitle('a','location','nw','fontweight','bold','fontsize',6)

% panel b
subplot(2,2,3); hold on
plot([-10 10],[-10 10],'-k');
p1 = plot(oib_dh_b62,cryosat_dh_b62,'.','color',rgb('ocean blue')); plot(cryosat_dh_b72,oib_dh_b72,'.','color',rgb('ocean blue'));
p2 = plot(oib_dh_b62(1:end-4),imau_fdm_dh_b62,'.','color',rgb('light teal')); plot(imau_fdm_dh_b72,oib_dh_b72(1:end-3),'.','color',rgb('light teal'));
axis([-7 0 -7 0]); daspect([1 1 1]);
xlabel('IceBridge (m)'); ylabel('CryoSat-2 & IMAU-FDM (m)');
set(gca,'layer','top'); box on;
set(gca,'fontsize',6)
ntitle('b','location','nw','fontweight','bold','fontsize',6)
ntitle({'RMS = 39 cm',''},'location','se','color',rgb('ocean blue'),'fontsize',6)
ntitle({'RMS = 36 cm'},'location','se','color',rgb('light teal'),'fontsize',6)

% panel c
% seasonal dh scatter
subplot(2,2,2); hold on
plot([-6 6],[-6 6],'k')
fill([-6 -6 0 0 -6],[-6 6 6 -6 -6],rgb('reddish'),'facealpha',.2,'edgecolor','none')
fill([0 0 6 6 0],[-6 6 6 -6 -6],rgb('bluish'),'facealpha',.2,'edgecolor','none')

cmapg = cmocean('gray',10);

for i = 1:8
    % plot each basin in turn
    if i == 1
        p1 = plot(cryosat_may_aug_dh_zwallybasins(i,1:6),imau_fdm_may_aug_dh_zwallybasins(i,:),'o','color',cmapg(i,:),'markersize',2.5);
        plot(cryosat_sep_apr_dh_zwallybasins(i,1:5),imau_fdm_sep_apr_dh_zwallybasins(i,:),'o','color',cmapg(i,:),'markersize',2.5);
    end

    if i == 2
        p2 = plot(cryosat_may_aug_dh_zwallybasins(i,1:6),imau_fdm_may_aug_dh_zwallybasins(i,:),'s','color',cmapg(i,:),'markersize',2.5);
        plot(cryosat_sep_apr_dh_zwallybasins(i,1:5),imau_fdm_sep_apr_dh_zwallybasins(i,:),'s','color',cmapg(i,:),'markersize',2.5);
    end

    if i == 3
        p3 = plot(cryosat_may_aug_dh_zwallybasins(i,1:6),imau_fdm_may_aug_dh_zwallybasins(i,:),'^','color',cmapg(i,:),'markersize',2.5);
        plot(cryosat_sep_apr_dh_zwallybasins(i,1:5),imau_fdm_sep_apr_dh_zwallybasins(i,:),'^','color',cmapg(i,:),'markersize',2.5);
    end

    if i == 4
        scatter(cryosat_may_aug_dh_zwallybasins(i,1:6),imau_fdm_may_aug_dh_zwallybasins(i,:),7,'d','filled','markerfacecolor',cmapg(i,:),'markerfacealpha',.3);
        scatter(cryosat_sep_apr_dh_zwallybasins(i,1:5),imau_fdm_sep_apr_dh_zwallybasins(i,:),7,'d','filled','markerfacecolor',cmapg(i,:),'markerfacealpha',.3);
        p4 = plot(NaN,NaN,'d','color',cmapg(i,:),'markersize',2.5,'markerfacecolor',cmapg(i,:))
    end

    if i == 5
        p5 = plot(cryosat_may_aug_dh_zwallybasins(i,1:6),imau_fdm_may_aug_dh_zwallybasins(i,:),'v','color',cmapg(i,:),'markersize',2.5);
        plot(cryosat_sep_apr_dh_zwallybasins(i,1:5),imau_fdm_sep_apr_dh_zwallybasins(i,:),'v','color',cmapg(i,:),'markersize',2.5);
    end

    if i == 6
        p6 = plot(cryosat_may_aug_dh_zwallybasins(i,1:6),imau_fdm_may_aug_dh_zwallybasins(i,:),'>','color',cmapg(i,:),'markersize',2.5);
        plot(cryosat_sep_apr_dh_zwallybasins(i,1:5),imau_fdm_sep_apr_dh_zwallybasins(i,:),'>','color',cmapg(i,:),'markersize',2.5);
    end

    if i == 7
        p7 = plot(cryosat_may_aug_dh_zwallybasins(i,1:6),imau_fdm_may_aug_dh_zwallybasins(i,:),'<','color',cmapg(i,:),'markersize',2.5);
        plot(cryosat_sep_apr_dh_zwallybasins(i,1:5),imau_fdm_sep_apr_dh_zwallybasins(i,:),'>','color',cmapg(i,:),'markersize',2.5);
    end

    if i == 8
        p8 = plot(cryosat_may_aug_dh_zwallybasins(i,1:6),imau_fdm_may_aug_dh_zwallybasins(i,:),'x','color',cmapg(i,:),'markersize',4);
        plot(cryosat_sep_apr_dh_zwallybasins(i,1:5),imau_fdm_sep_apr_dh_zwallybasins(i,:),'x','color',cmapg(i,:),'markersize',4);
    end
end
axis([-6 6 -6 6]); daspect([1 1 1]); box on
xlabel('CryoSat-2 (m)'); ylabel('IMAU-FDM (m)');
set(gca,'fontsize',6,'xtick',[-6:2:6])
% NB basin 4 ignored for stats
ntitle({'',' May-August',' RMS = 36 cm'},'location','nw','color',rgb('reddish'),'fontsize',5)
ntitle({'September-April','RMS = 45 cm ','',''},'location','se','color',rgb('bluish'),'fontsize',5)
ntitle('c','location','nw','fontweight','bold','fontsize',6)

l = legend([p1 p2 p3 p4 p5 p6 p7 p8],'1','2','3','4','5','6','7','8','location','eastoutside','fontsize',5);
title(l,'       Basin ID','fontweight','normal'); legend box off

% panel d
subplot(2,2,4); hold on
plot([-10 10],[-10 10],'-k');
p1 = plot(oib_dhdt_b62,cryosat_dhdt_b62,'o','color',rgb('ocean blue'),'markersize',2.5); plot(oib_dhdt_b72,cryosat_dhdt_b72,'o','color',rgb('ocean blue'),'markersize',2.5);
p2 = plot(oib_dhdt_b62(1:3),imau_fdm_dhdt_b62,'o','color',rgb('light teal'),'markersize',2.5); plot(oib_dhdt_b72(1:3),imau_fdm_dhdt_b72,'o','color',rgb('light teal'),'markersize',2.5);
axis([-2.5 2.5 -2.5 2.5]); daspect([1 1 1]);
set(gca,'xtick',[-2:1:2])
xlabel('IceBridge (m)'); ylabel('CryoSat-2 & IMAU-FDM (m)');
set(gca,'layer','top'); box on;
set(gca,'fontsize',6)
ntitle('d','location','nw','fontweight','bold','fontsize',6)
ntitle({'RMS = 29 cm',''},'location','se','color',rgb('ocean blue'),'fontsize',6)
ntitle({'RMS = 18 cm'},'location','se','color',rgb('light teal'),'fontsize',6)

return % manually adjust legend position
