clc;
clear all;
close all;
format long g;

%% STOCHASTIC PROCESSES
% Exercise 2 Task2
data2=load('Exercise2-2.txt');

t = data2(:,1);
y2_1 = data2(:,2);
y2_2 = data2(:,3);

% plot time series
figure
plot(t, y2_1, '.-')
hold on
plot(t, y2_2, '.-')
xlabel('t')
ylabel('y(t)')

% calc cross corr
[xcf, lags] = xcorr(y2_1, y2_2, 'normalized');
dt = t(2) - t(1);

% plot xcorr
figure
plot(lags*dt, xcf, '.-')
xlabel('lags')
ylabel('cross corr function')