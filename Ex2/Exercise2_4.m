clc;
clear;

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
plot(t_tls, y_tls, '.-')
hold on 
plot(t_ibis, y_ibis, '.-')
legend('TLS observations','IBIS observations')
title('Observation plots')
xlabel('Time [s]')
ylabel('Displacement [m]')

% determine time-shift

[xcf, lags] = xcorr(y_tls, y_ibis, 'normalized');
[maximum,delay] = max(abs(xcf));
shift_samples = lags(delay);


figure
plot(lags, xcf, '.')
xlabel('Shift / lags')
ylabel('Cross-corellation value')
title('Correllation function plot')

freq_est = length(t_tls) / max(t_tls) % 1/s
s_from_freq_est = 1/freq_est % s
total_shift_s = shift_samples * s_from_freq_est
sprintf('integer Shift until maximum corellation: %.i\n\ntime-adjusted shift until maximum corellation: %.3fs',shift_samples,total_shift_s);

% correct time-axis of TLS-measurements 
t_tls_corr = t_tls + 2*total_shift_s;

figure
hold on;
plot(t_tls,y_tls,'-')
plot(t_tls_corr,y_tls,'-')
plot(t_ibis, y_ibis,'-')
legend('TLS observations','shift-corrected TLS observations','original IBIS time series')
title('Observation plots (corrected)')
xlabel('Time [s]')
ylabel('Displacement [m]')
hold off;