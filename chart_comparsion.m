clc;
clear all;
close all;

figure
data_50_with = load('antenna-modeling/mat_files/antenna_50_with.mat');
data_75_without = load('antenna-modeling/mat_files/antenna_75_without.mat');
data_75_with = load('antenna-modeling/mat_files/antenna_75_with.mat');
data_50_without = load('antenna-modeling/mat_files/antenna_50_without.mat');
% plot reflection coefficient S11
plot( data_50_with.freq/1e6, (1+ abs(data_50_with.s11))./(1-abs(data_50_with.s11)), 'k-', 'Linewidth', 2 );
grid on;
hold on;
plot( data_50_without.freq/1e6, (1+ abs(data_50_without.s11))./(1-abs(data_50_without.s11)), 'k-.', 'Linewidth', 2 );
plot( data_75_with.freq/1e6, (1+ abs(data_75_with.s11))./(1-abs(data_75_with.s11)), 'k--', 'Linewidth', 2 );
plot( data_75_without.freq/1e6, (1+ abs(data_75_without.s11))./(1-abs(data_75_without.s11)), 'k:', 'Linewidth', 2 );
legend('a=50mm, com dielétrico', 'a=50mm, sem dielétrico','a=75mm, com dielétrico', 'a=75mm, com dielétrico');
title( 'VSWR' );
xlabel( 'frequency f / MHz' );
xlim ([500 4000])
ylim ([0 3])
