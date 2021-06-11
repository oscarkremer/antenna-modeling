
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
f0 = 2.25e9; % center frequency, frequency of interest!
lambda0 = round(c0/f0/unit); % wavelength in mm
%fc = 1.75e9; % 20 dB corner frequency
fc = 0.25e9; % 20 dB corner frequency
feed.heigth = 2;
rf = round(7.5/2);
Monocone.a = 50;
Monocone.theta0 = 33.2836*pi/180;
Helix.mesh_res = 4;
Monocone.sphere_radius = round(rf-1+(Monocone.a)*tan(Monocone.theta0));
Monocone.sphere_center = round((Monocone.a)/cos(Monocone.theta0))+feed.heigth;
lenz.epsR = 2.1;
lenz.kappa = 1.5e-4;
nr= sqrt(lenz.epsR);
rho_g = lenz_project(Monocone.a*unit, Monocone.theta0, lenz.epsR)/unit
diff_rho_a = rho_g - round(Monocone.a*sin(Monocone.theta0));
rho_1 = round(Monocone.a*sin(Monocone.theta0)) + diff_rho_a/4;
rho_2 = round(Monocone.a*sin(Monocone.theta0)) + 2*diff_rho_a/4;
rho_3 = round(Monocone.a*sin(Monocone.theta0)) + 3*diff_rho_a/4;
rho_4 = round(Monocone.a*sin(Monocone.theta0)) + 8*diff_rho_a/10;
rho_5 = round(Monocone.a*sin(Monocone.theta0)) + 9*diff_rho_a/10;
z1 = sqrt((nr^2-1)/((nr+1)^2)*rho_g*rho_g - ((rho_1 - (rho_g/(nr+1)))/(nr/(sqrt(nr^2-1))))^2)
z2 = sqrt((nr^2-1)/((nr+1)^2)*rho_g*rho_g - ((rho_2 - (rho_g/(nr+1)))/(nr/(sqrt(nr^2-1))))^2)
z3 = sqrt((nr^2-1)/((nr+1)^2)*rho_g*rho_g - ((rho_3 - (rho_g/(nr+1)))/(nr/(sqrt(nr^2-1))))^2)
z4 = sqrt((nr^2-1)/((nr+1)^2)*rho_g*rho_g - ((rho_4 - (rho_g/(nr+1)))/(nr/(sqrt(nr^2-1))))^2)
z5 = sqrt((nr^2-1)/((nr+1)^2)*rho_g*rho_g - ((rho_5 - (rho_g/(nr+1)))/(nr/(sqrt(nr^2-1))))^2)

gnd.radius = 2*lambda0;

% feeding
feed.R = 50;    %feed impedance

% size of the simulation box
SimBox = [2 2 2]*2*lambda0;

%% setup FDTD parameter & excitation function
FDTD = InitFDTD( );
FDTD = SetGaussExcite( FDTD, f0, fc );
BC = {'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8' 'PML_8'}; % boundary conditions
FDTD = SetBoundaryCond( FDTD, BC );

%% setup CSXCAD geometry & mesh
max_res = floor(c0 / (f0+fc) / unit / 20); % cell size: lambda/20
CSX = InitCSX();

% create helix mesh
mesh.x = SmoothMeshLines([-rho_g-rf -Monocone.sphere_radius-rf 0 Monocone.sphere_radius+rf rho_g+rf], Helix.mesh_res);
% add the air-box
mesh.x = [mesh.x -SimBox(1)/2-gnd.radius  SimBox(1)/2+gnd.radius];
% create a smooth mesh between specified fixed mesh lines
mesh.x = SmoothMeshLines( mesh.x, max_res, 1.4);

% copy x-mesh to y-direction
mesh.y = mesh.x;

% create helix mesh in z-direction
mesh.z = SmoothMeshLines([0 feed.heigth round(Monocone.a*cos(Monocone.theta0))+feed.heigth Monocone.sphere_center Monocone.sphere_center+Monocone.sphere_radius], Helix.mesh_res);
% add the air-box
mesh.z = unique([mesh.z -SimBox(3)/2 max(mesh.z)+SimBox(3)/2]);
% create a smooth mesh between specified fixed mesh lines
mesh.z = SmoothMeshLines( mesh.z, max_res, 1.4 );

CSX = DefineRectGrid( CSX, unit, mesh );

%% create helix using the wire primitive
CSX = AddMetal( CSX, 'helix' ); % create a perfect electric conductor (PEC)

clear p;
rf = round(7.5/2);
p(1,1) = feed.heigth; p(2,1) = 0;
p(1,2) = feed.heigth; p(2,2) = rf;
p(1,3) = round(Monocone.a*cos(Monocone.theta0))+feed.heigth; p(2,3) = round(Monocone.a*sin(Monocone.theta0))+rf;
p(1,4) = Monocone.a+feed.heigth; p(2,4) = 0;


%p(1,3) = round(Monocone.a*cos(Monocone.theta0)); p(2,3) = 0

CSX = AddRotPoly( CSX, 'helix', 0, 'y', p, 'z', [0,2*pi]);


CSX = AddSphere(CSX, 'helix', 0, [0 0 Monocone.sphere_center], Monocone.sphere_radius);


CSX = AddMaterial( CSX, 'lenz' ); % create a perfect electric conductor (PEC)


clear p;
p(1,1) = 0; p(2,1) = 1;
p(1,2) = feed.heigth; p(2,2) = 1;
p(1,3) = feed.heigth; p(2,3) = rf;
p(1,4) = round(Monocone.a*cos(Monocone.theta0))+feed.heigth; p(2,4) = round(Monocone.a*sin(Monocone.theta0))+rf;
p(1,5) = z1+feed.heigth; p(2,5) = rho_1+rf;
p(1,6) = z2+feed.heigth; p(2,6) = rho_2+rf;
p(1,7) = z3+feed.heigth; p(2,7) = rho_3+rf;
p(1,8) = z4+feed.heigth; p(2,8) = rho_4+rf;
p(1,9) = z5+feed.heigth; p(2,9) = rho_5+rf;
p(1,10) = feed.heigth; p(2,10) = rho_g+rf;
p(1,11) = 0; p(2,11) = rho_g+rf;


CSX = SetMaterialProperty(CSX, 'lenz', 'Epsilon', lenz.epsR, 'Kappa', lenz.kappa);


CSX = AddRotPoly( CSX, 'lenz', 2, 'y', p, 'z', [0,2*pi]);



%% create ground circular ground
CSX = AddMetal( CSX, 'gnd' ); % create a perfect electric conductor (PEC)
% add a box using cylindrical coordinates
start = [-gnd.radius -gnd.radius 0];
stop  = [gnd.radius gnd.radius 0];
CSX = AddBox(CSX,'gnd',10,start, stop);

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
xlim ([500 4000])
ylim ([0 3])

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


