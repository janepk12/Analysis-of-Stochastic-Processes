clc;
clear all;
close all;
format long g;

%% STOCHASTIC PROCESSES
% Exercise 2 Task1

data = load('Exercise2-1.txt');

t = data(:,1);
y1 = data(:,2);
y2 = data(:,3);

figure
plot(t, y1, '.-')
hold on
plot(t, y2, '.-')
xlabel('t [s]')
ylabel('y')
legend('y2', 'y1')

% xcorr means we consider equally spaced time measurements of both time series
[xcf, lags] = xcorr(y1, y2, 'normalized'); 

dt = t(2)-t(1);

figure
plot(lags*dt, xcf, '.-')
xlabel('lags [s]')
ylabel('cross correlation function')

