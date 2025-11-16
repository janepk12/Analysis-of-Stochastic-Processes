clc;
clear all;
close all;
format long g;

%% STOCHASTIC PROCESSES
% Exercise 2 Task3 

temp_data = load('temperature.txt');
incl_data = load('inclination.txt');

figure(1)
subplot(211)
yyaxis left;
plot(temp_data(:,1), temp_data(:,2), '.-')
ylabel('temp [°C]')
hold on
yyaxis right;
plot(incl_data(:,1), incl_data(:,2), '.-')
xlabel('t [min]')
ylabel('inc [mgon]')

% t = incl_data(:,1);
% dt = t(2)-t(1);
[xcf, lags] = xcorr(temp_data(:,2)-mean(temp_data(:,2)), incl_data(:,2)-mean(incl_data(:,2)), 'normalized');
subplot(212)
plot(lags, xcf, '.')
