
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
Monocone.width = 5
lambda0 = 4*Monocone.a;
f0 = (3*10^8)/(lambda0*unit);
fc = 0.5e9; % 20 dB corner frequency

Helix.mesh_res = 1;
gnd.radius = 1*lambda0;
% feeding
feed.R = 50;    %feed impedance
feed.heigth = 3;    %feed impedance

% size of the simulation box
SimBox = [2 2 3]*lambda0;

%% setup FDTD parameter & excitation function
FDTD = InitFDTD( );
FDTD = SetGaussExcite( FDTD, f0, fc );
BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'}; % boundary conditions
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh
max_res = floor(c0 / (f0+fc) / unit / 25); % cell size: lambda/20
CSX = InitCSX();

% create helix mesh
mesh.x = SmoothMeshLines([-round(Monocone.a*sin(Monocone.theta0)) 0 round(Monocone.a*sin(Monocone.theta0))], Helix.mesh_res);
% add the air-box
mesh.x = [mesh.x -SimBox(1)/2-gnd.radius SimBox(1)/2+gnd.radius];
% create a smooth mesh between specified fixed mesh lines
mesh.x = SmoothMeshLines( mesh.x, max_res, 1.4);

% copy x-mesh to y-direction
mesh.y = mesh.x;

% create helix mesh in z-direction
%mesh.z = [mesh.z SimBox(3)];

mesh.z = SmoothMeshLines([0 feed.heigth feed.heigth+round(Monocone.a)],Helix.mesh_res);

mesh.z = unique([mesh.z -SimBox(3)/4 max(mesh.z)+SimBox(3)/2 ]);
% create a smooth mesh between specified fixed mesh lines
mesh.z = SmoothMeshLines( mesh.z, max_res, 1.4 );

CSX = DefineRectGrid( CSX, unit, mesh );

%% create helix using the wire primitive
CSX = AddMetal( CSX, 'monocone' ); % create a perfect electric conductor (PEC)

clear p;
p(1,1) = feed.heigth; p(2,1) = 0;
p(1,2) = round(Monocone.a*cos(Monocone.theta0))+feed.heigth; p(2,2) = round(Monocone.a*sin(Monocone.theta0));
p(1,3) = round(Monocone.a*cos(7*Monocone.theta0/8))+feed.heigth; p(2,3) = round(Monocone.a*sin(7*Monocone.theta0/8));
p(1,4) = round(Monocone.a*cos(6*Monocone.theta0/8))+feed.heigth; p(2,4) = round(Monocone.a*sin(6*Monocone.theta0/8));
p(1,5) = round(Monocone.a*cos(5*Monocone.theta0/8))+feed.heigth; p(2,5) = round(Monocone.a*sin(5*Monocone.theta0/8));
p(1,6) = round(Monocone.a*cos(4*Monocone.theta0/8))+feed.heigth; p(2,6) = round(Monocone.a*sin(4*Monocone.theta0/8));
p(1,7) = round(Monocone.a*cos(3*Monocone.theta0/8))+feed.heigth; p(2,7) = round(Monocone.a*sin(3*Monocone.theta0/8));
p(1,8) = round(Monocone.a*cos(2*Monocone.theta0/8))+feed.heigth; p(2,8) = round(Monocone.a*sin(2*Monocone.theta0/8));
p(1,9) = Monocone.a+feed.heigth; p(2,9) = 0;
p(1,10) = Monocone.a-Monocone.width+feed.heigth; p(2,10) = 0;
p(1,11) = round((-2*Monocone.width+Monocone.a)*cos(2*Monocone.theta0/8))+Monocone.width+feed.heigth; p(2,11) = round((-2*Monocone.width+Monocone.a)*sin(2*Monocone.theta0/8));
p(1,12) = round((-2*Monocone.width+Monocone.a)*cos(3*Monocone.theta0/8))+Monocone.width+feed.heigth; p(2,12) = round((-2*Monocone.width+Monocone.a)*sin(3*Monocone.theta0/8));
p(1,13) = round((-2*Monocone.width+Monocone.a)*cos(4*Monocone.theta0/8))+Monocone.width+feed.heigth; p(2,13) = round((-2*Monocone.width+Monocone.a)*sin(4*Monocone.theta0/8));
p(1,14) = round((-2*Monocone.width+Monocone.a)*cos(5*Monocone.theta0/8))+Monocone.width+feed.heigth; p(2,14) = round((-2*Monocone.width+Monocone.a)*sin(5*Monocone.theta0/8));
p(1,15) = round((-2*Monocone.width+Monocone.a)*cos(6*Monocone.theta0/8))+Monocone.width+feed.heigth; p(2,15) = round((-2*Monocone.width+Monocone.a)*sin(6*Monocone.theta0/8));
p(1,16) = round((-2*Monocone.width+Monocone.a)*cos(7*Monocone.theta0/8))+Monocone.width+feed.heigth; p(2,16) = round((-2*Monocone.width+Monocone.a)*sin(7*Monocone.theta0/8));
p(1,17) = round((-2*Monocone.width+Monocone.a)*cos(Monocone.theta0))+Monocone.width+feed.heigth; p(2,17) = round((-2*Monocone.width+Monocone.a)*sin(Monocone.theta0));
p(1,18) = feed.heigth+Monocone.width; p(2,18) = 0;


%p(1,3) = round(Monocone.a*cos(Monocone.theta0)); p(2,3) = 0

CSX = AddRotPoly( CSX, 'monocone', 0, 'y', p, 'z', [0,2*pi]);

%% create ground circular ground
CSX = AddMetal( CSX, 'gnd' ); % create a perfect electric conductor (PEC)
% add a box using cylindrical coordinates
start = [4         0    feed.heigth];
stop  = [gnd.radius 2*pi 0];
CSX = AddBox(CSX,'gnd',10,start,stop,'CoordSystem',1);
start = [0         0    0];
stop  = [gnd.radius 2*pi -1];
CSX = AddBox(CSX,'gnd',10,start,stop,'CoordSystem',1);

%% apply the excitation & resist as a current source
start = [0 0 0];
stop  = [0 0 feed.heigth];
[CSX port] = AddLumpedPort(CSX, 5 ,1 ,feed.R, start, stop, [0 0 1], true);

%%nf2ff calc
start = [mesh.x(11)      mesh.y(11)     mesh.z(11)];
stop  = [mesh.x(end-10) mesh.y(end-10) mesh.z(end-10)];
[CSX nf2ff] = CreateNF2FFBox(CSX, 'nf2ff', start, stop, 'OptResolution', lambda0/15);

%% prepare simulatio3 folder
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

% plot reflection coefficient S11
figure
plot( freq/1e9, 20*log10(abs(s11)), 'k-', 'Linewidth', 2 );
grid on
title( 'reflection coefficient S_{11}' );
xlabel( 'frequency f / MHz' );
ylabel( 'reflection coefficient |S_{11}|' );

figure
plot( freq/1e9, real(60*log(abs(cot(Monocone.theta0/2)))*(1+s11)./(1-s11)), 'k-', 'Linewidth', 2 );
hold on
grid on
plot( freq/1e9, imag(60*log(abs(cot(Monocone.theta0/2)))*(1+s11)./(1-s11)), 'r--', 'Linewidth', 2 );
title( 'Antenna impedance' );
xlabel( 'frequency f / GHz' );
ylabel( 'impedance Z_{in} / Ohm' );
legend( 'real', 'imag' );

figure
plot( freq/1e9, real(Zin), 'k-', 'Linewidth', 2 );
hold on
grid on
plot( freq/1e9, imag(Zin), 'r--', 'Linewidth', 2 );
title( 'Antenna impedance' );
xlabel( 'frequency f / GHz' );
ylabel( 'impedance Z_{in} / Ohm' );
legend( 'real', 'imag' );



drawnow

%% NFFF contour plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find resonance frequncy from s11
f_res = f0;

% get accepted antenna power at frequency f0
P_in_0 = interp1(freq, port.P_acc, f0);

% calculate the far field at phi=0 degrees and at phi=90 degrees
thetaRange = unique([0:0.5:90 90:180]);
phiRange = (0:2:360) - 180;
disp( 'calculating the 3D far field...' );

nf2ff = CalcNF2FF(nf2ff, Sim_Path, f_res, thetaRange*pi/180, phiRange*pi/180,'Mode',0,'Outfile','3D_Pattern.h5','Verbose',1);

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
