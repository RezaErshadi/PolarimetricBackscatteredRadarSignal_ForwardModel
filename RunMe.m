%%
% This code can model the backscattered radio signal for isotropic and anisotropic ice column 
% according to:
% Fujita et al. (2006)
% Ershadi et al. (2022) 
% written by M. Reza Ershadi (Geophysics & Glaciology - University of Tübingen)
%%
clear
close all
clc
obj = FM_class;
%% Primary parameters
obj.d_depth = 1; % vertical resolution [m]
obj.d_alpha_mdl = 0; % the angle between measurement frame and T antenna
obj.d_alpha_synt = 0:1:179; % horizontal resolution [°] (the angle between measurement frame and T antenna)
obj.t_layer = [2001]; % depth of the horizontal layer boundary
obj.theta = [45]; % fabric orientation (abgle between v1 and TR)
obj.AzOfst = 0; % Offset azimuth for theta
obj.HA = [0.026]; % Horizintal anisotropy = (largest horizontal eigenvalue or Y - largest horizontal eigenvalue or X)
obj.rdB = [0]; % reflection ratio in dB (0 = no anisotropic scattering)
obj.f = 300e6; % radar center frequency [Hz]
%% Constant parameters
obj.epsilon0 = 8.85e-12; % dielectric permitivity in a vacuum [F/m]
obj.DelEpsP = 0.034; % Delta Epsilon Prime -- > dielectric anisotropy 
obj.c = 299792458; % speed of light [m/s]
obj.mu0 = pi * 4e-7; % magnetic permeability in a vacuum [Henry/m]
obj.epsX = 3.15; % absolute value of dielectric permitivity in X
obj.GammaX = 1e-12; % absolute value of Complex amplitudes reflection coefficients in X
obj.delX = 1e-5; % absolute value of Conductivity X
obj.delY = 1e-5; % absolute value of Conductivity Y
%% Secondary parameters
obj = FM_class.CalculateSecondaryParameters(obj);
%% Run The Model
obj = FM_class.RunTheModel(obj);
%% Synthesize
obj = FM_class.AzimuthSynthesizer(obj);
%% Extract parameters
obj = FM_class.Signal2Param(obj);
alpha = obj.d_alpha_synt;
PAHH = obj.Dta{5};
PAHV = obj.Dta{7};
cohPhs = obj.Dta{14};
Psi = obj.Dta{18};
Z = obj.Z;
%% Plot
fig = FM_class.SetFigureSize(0,0,0.7,0.65);
ax{1} = subplot(1,3,1);
ax{2} = subplot(1,3,2);
ax{3} = subplot(1,3,3);
load('SeismicColorMap100.mat');
% Plot HH Power Anomaly
ii = 1;
p1 = pcolor(ax{ii},alpha,Z,PAHH);
set(p1, 'EdgeColor', 'none');
set(ax{ii},'YDir','reverse')
colorbar(ax{ii})
colormap(ax{ii},SeismicColorMap100)
caxis(ax{ii},[-5 5])
xticks(ax{ii},[45:45:180])
title(ax{ii},'Power Anomaly HH')
xlabel(ax{ii},'\alpha [deg]')
ylabel(ax{ii},'Depth [m]')
set(ax{ii},'FontSize',18,'Layer','top');
xlim(ax{ii},[0 179])
ylim(ax{ii},[Z(1) Z(end)])
% Plot HV Power Anomaly
ii = 2;
p1 = pcolor(ax{ii},alpha,Z,PAHV);
set(p1, 'EdgeColor', 'none');
set(ax{ii},'YDir','reverse')
colorbar(ax{ii})
colormap(ax{ii},SeismicColorMap100)
caxis(ax{ii},[-5 5])
xticks(ax{ii},[45:45:180])
title(ax{ii},'Power Anomaly HV')
xlabel(ax{ii},'\alpha [deg]')
ylabel(ax{ii},'Depth [m]')
set(ax{ii},'FontSize',18,'Layer','top');
xlim(ax{ii},[0 179])
ylim(ax{ii},[Z(1) Z(end)])
% Plot HHVV Coherence Phase
ii = 3;
p1 = pcolor(ax{ii},alpha,Z,cohPhs);
set(p1, 'EdgeColor', 'none');
set(ax{ii},'YDir','reverse')
colorbar(ax{ii})
colormap(ax{ii},SeismicColorMap100)
caxis(ax{ii},[-pi pi])
xticks(ax{ii},[45:45:180])
title(ax{ii},'Coherence Phase')
xlabel(ax{ii},'\alpha [deg]')
ylabel(ax{ii},'Depth [m]')
set(ax{ii},'FontSize',18,'Layer','top');
xlim(ax{ii},[0 179])
ylim(ax{ii},[Z(1) Z(end)])
%% Save Plot
saveasname = 'FM';
f.InvertHardcopy = 'off';
% print(f,saveasname+".png",'-dpng','-r300');