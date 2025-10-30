%
%
% Exercise 1 of Analysis of Stochastic Processes
% Group 2
%


% Task 1
clc;
clear;
close all;  

% Loading in the data
data1 = load('/Users/user/Desktop/GGS_Master/GGS_WiSe25_26/Stochastic Processes/Analysis-of-Stochastic-Processes/Ex1/Exercise1-1.txt');
x1 = data1(:,1);
y1 = data1(:,2);
y2 = data1(:,3);
y3 = data1(:,4);

% Plotting everything in same window
hold on;
figure(1);
plot(x1,y1,'.')

plot(x1,y2,'.')

plot(x1,y3,'.')

legend('Time Series 1','Time Series 2','Time Series 3')


hold off;

