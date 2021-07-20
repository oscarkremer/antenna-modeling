clc;
clear all;
close all;

addpath('~/opt/openEMS/share/openEMS/matlab');
addpath('~/opt/openEMS/share/CSXCAD/matlab');
addpath('~/opt/openEMS/share/hyp2mat/matlab'); % hyp2mat package
addpath('~/opt/openEMS/share/CTB/matlab'); % circuit toolbox
data = load('antenna_75_with.mat');
data.nf2ff
freq = data.freq;
f_res = data.f0;
freq_points = size(freq)(2);

thetaRange = unique([0:0.5:180]);
phiRange = (0:2:360) - 180;
%freq_plot = [];
%gain = [];
P_in_0 = interp1(freq, data.port.P_acc, f_res);

nf2ff = CalcNF2FF(data.nf2ff, data.Sim_Path, f_res, thetaRange*pi/180, phiRange*pi/180,'Mode',0,'Outfile','3D_Pattern.h5','Verbose',1);
Dlog = 10*log10(nf2ff.Dmax);
efficiency_log = 10*log10(nf2ff.Prad/P_in_0);
%Gain = Dlog + efficiency_log;
%gain = [gain Gain]ff

