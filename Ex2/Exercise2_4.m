clc;
clear all;
close all;
format long g;

%% STOCHASTIC PROCESSES
% Exercise 2 Task4

% plot time series
tls_data = load("tls.txt");
ibis_data = load("ibis.txt");

t_tls = tls_data(:,1);
t_ibis = ibis_data(:,1);
y_tls = tls_data(:,2);
y_ibis = ibis_data(:,2);

figure
plot(t, y_tls, '.-')
hold on 
plot(t, y_ibis, '.-')
xlabel('Time [s]')
ylabel('Displacement [m]')

% determine time-shift
[xcf, lags] = xcorr(y_tls, y_ibis, 'normalized');

figure
plot(lags, xcf, '.')

% correct time-axis of TLS-measurements 
% t_tls = ...