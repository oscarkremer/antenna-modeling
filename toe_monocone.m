
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

post_proc_only = 0;

close all

%% setup the simulation
physical_constants;
unit = 1e-3; % all length in mm
Monocone.a = 50;
Monocone.theta0 = 46.9795*pi/180;
f0 = 2.25e9; % center frequency, frequency of interest!
lambda0 = round(c0/f0/unit); % wavelength in mm
fc = 1.75e9; % 20 dB corner frequency

Helix.radius = 20; % --> diameter is ~ lambda/pi
Helix.turns = 10;  % --> expected gain is G ~ 4 * 10 = 40 (16dBi)
Helix.pitch = 30;  % --> pitch is ~ lambda/4
Helix.mesh_res = 4;

gnd.radius = 3*lambda0;

% feeding
feed.heigth = 1;
feed.R = 50;    %feed impedance

% size of the simulation box
SimBox = [1 1 1.5]*2.5*lambda0;

%% setup FDTD parameter & excitation function
FDTD = InitFDTD( );
FDTD = SetGaussExcite( FDTD, f0, fc );
BC = {'MUR' 'MUR' 'MUR' 'MUR' 'MUR' 'PML_8'}; % boundary conditions
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh
max_res = floor(c0 / (f0+fc) / unit / 20); % cell size: lambda/20
CSX = InitCSX();

% create helix mesh
mesh.x = SmoothMeshLines([-2*Monocone.a 0 2*Monocone.a], Helix.mesh_res);
% add the air-box
mesh.x = [mesh.x -SimBox(1)/2-gnd.radius  SimBox(1)/2+gnd.radius];
% create a smooth mesh between specified fixed mesh lines
mesh.x = SmoothMeshLines( mesh.x, max_res, 1.4);

% copy x-mesh to y-direction
mesh.y = mesh.x;

% create helix mesh in z-direction
mesh.z = SmoothMeshLines([0 feed.heigth 3*Monocone.a+feed.heigth],Helix.mesh_res);
% add the air-box
mesh.z = unique([mesh.z -SimBox(3)/3 max(mesh.z)+2*SimBox(3)/3]);
% create a smooth mesh between specified fixed mesh lines
mesh.z = SmoothMeshLines( mesh.z, max_res, 1.4 );

CSX = DefineRectGrid( CSX, unit, mesh );

%% create helix using the wire primitive
CSX = AddMetal( CSX, 'helix' ); % create a perfect electric conductor (PEC)


rc1 = round(93/2);
rc2 = round(102/2);
rc3 = round(108/2);
rc4 = round(111.4/2);
rc5 = round(112.5/2);
rc6 = round(110.7/2);
rc7 = round(105.2/2);
rc8 = round(95.38/2);
rc9 = round(79.55/2);
rf = round(7.5/2);
h1 = 43;
h2 = 54;
h3 = 63;
h4 = 76;
h5 = 87;
h6 = 97;
h7 = 107;
h8 = 117;
h9 = 127;
ht = 175;
clear p;
p(1,1) = feed.heigth; p(2,1) = rf;
p(1,2) = h1 +feed.heigth; p(2,2) = rc1;
p(1,3) = h2 +feed.heigth; p(2,3) = rc2;
p(1,4) = h3 +feed.heigth; p(2,4) = rc3;
p(1,5) = h4 +feed.heigth; p(2,5) = rc4;
p(1,6) = h5 +feed.heigth; p(2,6) = rc5;
p(1,7) = h6 +feed.heigth; p(2,5) = rc6;
p(1,8) = h7 +feed.heigth; p(2,6) = rc7;
p(1,9) = h8 +feed.heigth; p(2,5) = rc8;
p(1,10) = h9 +feed.heigth; p(2,6) = rc9;
p(1,11) = ht +feed.heigth; p(2,5) = 0;
CSX = AddRotPoly( CSX, 'helix', 0, 'y', p, 'z', [0,2*pi]);




%% create ground circular ground
CSX = AddMetal( CSX, 'gnd' ); % create a perfect electric conductor (PEC)
% add a box using cylindrical coordinates
start = [0          0    0];
stop  = [gnd.radius 2*pi 0];
CSX = AddBox(CSX,'gnd',10,start,stop,'CoordSystem',1);

%% apply the excitation & resist as a current source
start = [0 0 0];
stop  = [0 0 feed.heigth];
[CSX port] = AddLumpedPort(CSX, 5 ,1 ,feed.R, start, stop, [0 0 1], true);

%%nf2ff calc
start = [mesh.x(11)      mesh.y(11)     mesh.z(11)];
stop  = [mesh.x(end-10) mesh.y(end-10) mesh.z(end-10)];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop, 'OptResolution', lambda0/15);

%% prepare simulation folder
Sim_Path = 'tmp_Helical_Ant';
Sim_CSX = 'Helix_Ant.xml';

if (post_proc_only==0)
    [status, message, messageid] = rmdir( Sim_Path, 's' ); % clear previous directory
    [status, message, messageid] = mkdir( Sim_Path );      % create empty simulation folder

    %% write openEMS compatible xml-file
    WriteOpenEMS( [Sim_Path '/' Sim_CSX], FDTD, CSX );

    %% show the structure
    CSXGeomPlot( [Sim_Path '/' Sim_CSX] );

    %% run openEMS
    RunOpenEMS( Sim_Path, Sim_CSX);
end

%% postprocessing & do the plots
freq = linspace( f0-fc, f0+fc, 501 );
port = calcPort(port, Sim_Path, freq);

Zin = port.uf.tot ./ port.if.tot;
s11 = port.uf.ref ./ port.uf.inc;

% plot feed point impedance
figure
plot( freq/1e6, real(Zin), 'k-', 'Linewidth', 2 );
hold on
grid on
plot( freq/1e6, imag(Zin), 'r--', 'Linewidth', 2 );
title( 'feed point impedance' );
xlabel( 'frequency f / MHz' );
ylabel( 'impedance Z_{in} / Ohm' );
legend( 'real', 'imag' );

% plot reflection coefficient S11
figure
plot( freq/1e6, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
grid on
title( 'reflection coefficient S_{11}' );
xlabel( 'frequency f / MHz' );
ylabel( 'reflection coefficient |S_{11}|' );

% plot reflection coefficient S11
figure
plot( freq/1e6, (1+ abs(s11))./(1-abs(s11)), 'k-', 'Linewidth', 2 );
grid on
title( 'reflection coefficient S_{11}' );
xlabel( 'frequency f / MHz' );
ylabel( 'reflection coefficient |S_{11}|' );


drawnow

%% NFFF contour plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find resonance frequncy from s11
f_res = f0;

% get accepted antenna power at frequency f0
P_in_0 = interp1(freq, port.P_acc, f0);

% calculate the far field at phi=0 degrees and at phi=90 degrees
thetaRange = unique([0:0.5:180]);
phiRange = (0:2:360) - 180;
disp( 'calculating the 3D far field...' );

nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_res, thetaRange*pi/180, phiRange*pi/180,'Mode',0,'Outfile','3D_Pattern.h5','Verbose',1);
figure
plotFF3D(nf2ff);        % plot liear 3D far field
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
