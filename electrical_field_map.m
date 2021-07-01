clc;
clear all;
close all;

figure
data_50_with = load('antenna-modeling/mat_files/antenna_50_with.mat');

%
% Tutorials / helical antenna
%
% Describtion at:
% http://openems.de/index.php/Tutorial:_Helical_Antenna
%
% Tested with
%  - Matlab 2011a / Octave 4.0
%  - openEMS v0.0.33
%
% (C) 2012-2015 Thorsten Liebig <thorsten.liebig@uni-due.de>

close all
clear
clc



addpath('~/opt/openEMS/share/openEMS/matlab');
addpath('~/opt/openEMS/share/CSXCAD/matlab');
addpath('~/opt/openEMS/share/hyp2mat/matlab'); % hyp2mat package
addpath('~/opt/openEMS/share/CTB/matlab'); % circuit toolbox

post_proc_only = 0;

close all

%% setup the simulation
physical_constants;
unit = 1e-3; % all length in mm
f0 = 2.25e9; % center frequency, frequency of interest!
lambda0 = round(c0/f0/unit); % wavelength in mm
fc = 1.75e9; % 20 dB corner frequency
feed.heigth = 2;
%rf = round(7.5/2);
rf = 2;
Monocone.a = 50;
Monocone.theta0 = 34*pi/180;
Monocone.sphere_radius = (rf+Monocone.a*sin(Monocone.theta0))/cos(Monocone.theta0)
Monocone.sphere_center = Monocone.a*cos(Monocone.theta0)+Monocone.sphere_radius*sin(Monocone.theta0)+feed.heigth
Helix.mesh_res = 3;

drawnow
%% NFFF contour plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find resonance frequncy from s11
f_res = f0;

% get accepted antenna power at frequency f0

% calculate the far field at phi=0 degrees and at phi=90 degrees
thetaRange = unique([0:0.5:180]);
phiRange = (0:2:360) - 180;
disp( 'calculating the 3D far field...' );

nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_res, thetaRange*pi/180, phiRange*pi/180,'Mode',0,'Outfile','3D_Pattern.h5','Verbose',1);
figure
plotFF3D(nf2ff,'logscale',-20);

  % plot liear 3D far field
theta_HPBW = interp1(nf2ff.E_norm{1}(:,1)/max(nf2ff.E_norm{1}(:,1)),thetaRange,1/sqrt(2))*2;

% display power and directivity
disp( ['radiated power: Prad = ' num2str(nf2ff.Prad) ' Watt']);
disp( ['directivity: Dmax = ' num2str(nf2ff.Dmax) ' (' num2str(10*log10(nf2ff.Dmax)) ' dBi)'] );
disp( ['efficiency: nu_rad = ' num2str(100*nf2ff.Prad./P_in_0) ' %']);
disp( ['theta_HPBW = ' num2str(theta_HPBW) ' °']);


%%
directivity = nf2ff.P_rad{1}/nf2ff.Prad*4*pi;
directivity_CPRH = abs(nf2ff.E_cprh{1}).^2./max(nf2ff.E_norm{1}(:)).^2*nf2ff.Dmax;
directivity_CPLH = abs(nf2ff.E_cplh{1}).^2./max(nf2ff.E_norm{1}(:)).^2*nf2ff.Dmax;

%%
figure
plot(thetaRange, 10*log10(directivity(:,1)'),'k-','LineWidth',2);
hold on
grid on
xlabel('theta (deg)');
ylabel('directivity (dBi)');
plot(thetaRange, 10*log10(directivity_CPRH(:,1)'),'g--','LineWidth',2);
plot(thetaRange, 10*log10(directivity_CPLH(:,1)'),'r-.','LineWidth',2);
legend('norm','CPRH','CPLH');

%% dump to vtk
DumpFF2VTK([Sim_Path '/3D_Pattern.vtk'],directivity,thetaRange,phiRange,'scale',1e-3);
DumpFF2VTK([Sim_Path '/3D_Pattern_CPRH.vtk'],directivity_CPRH,thetaRange,phiRange,'scale',1e-3);
DumpFF2VTK([Sim_Path '/3D_Pattern_CPLH.vtk'],directivity_CPLH,thetaRange,phiRange,'scale',1e-3);