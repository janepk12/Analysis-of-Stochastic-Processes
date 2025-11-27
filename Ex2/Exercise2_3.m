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
title('Observation data for tempreature and inclination')


[xcf, lags] = xcorr(temp_data(:,2)-mean(temp_data(:,2)), incl_data(:,2)-mean(incl_data(:,2)), 'normalized');
subplot(212)
plot(lags, xcf, '.')
title('Correllation function')
ylabel('Correllation value')
xlabel('Time steps/shifts')


[max, max_idx] = max(xcf)
lag_steps = floor(length(xcf)/2 - max_idx) % the difference between all observations and the 
timestep = lag_steps*5 % converting steps into time

% Adjust the time vector by subtracting the lag
%temp_corr = temp_data(:,1) + timestep;
incl_corr = incl_data(:,1) - timestep;



%% PLOTTING THE RESULTS
figure;
% plot(temp_corr, temp_data(:,2)-mean(temp_data(:,2)), 'b-')
plot(incl_corr, incl_data(:,2)-mean(incl_data(:,2)), '-','Color',[0.80, 0.25, 0.25])
hold on;

plot(incl_data(:,1), incl_data(:,2)-mean(incl_data(:,2)), '-', 'Color',[0.90, 0.60, 0.60])
plot(temp_data(:,1), temp_data(:,2)-mean(temp_data(:,2)), 'b--')
legend('Inclination corrected', 'Inclination original','Temperature original')
title('Corrected and uncorrected inclination and temperature data')
xlabel('Time [min]')
ylabel('Observation')
hold off;
