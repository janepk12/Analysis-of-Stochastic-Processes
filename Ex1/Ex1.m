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

%plot(x1,y1,'.')
cor1 = corr(x1,y1);

%plot(x1,y2,'.')
cor2 = corr(x1,y2);

%plot(x1,y3,'.')
cor3 = corr(x1,y3);

leg1 = sprintf('Time-Series1 -> corr: %.3f',cor1);
leg2= sprintf('Time-Series2 -> corr: %.3f',cor2);
leg3 = sprintf('Time-Series3 -> corr: %.3f',cor3);

title('Plots and correlation coefficient of 3 different time-series')
xlabel('x-values from column 1 from data')
ylabel('y-values from column 2,3,4 from data')
legend(leg1,leg2,leg3)

hold off;


% Task 1.5
% Aiming for these correlations:
% 1) 0.5  
% 2) -0.8 
% 3) 2

% Noise generation:
n = length(x1);
a = 2.0;
b = 5.0;
noise_level = 0.5;

x = linspace(0, 10, n);
y_linear = a * x + b;
noise = noise_level * randn(n,1); % Gaussian noise
yy1 = y_linear + noise;

%figure(2)
%cor11 = corr(x1,yy1)
%plot(x1 , yy1 , '.')




%% Task2

task2_1 = load("Ex1/Exercise1-2.txt");
task2_2 =load("Ex1/Exercise1-3.txt");


moving_mean1 = [];
moving_mean2 = [];
moving_std1 = [];
moving_std2 = [];
count = 0;
stepsize = 10;

if length(task2_1) == length(task2_1)
    obs_len = length(task2_2);
end


for i = 1:obs_len
    if mod(i,stepsize) == 0
        % take the previous 9 samples plus the current one to form a block of 10
        range1 = task2_1(i-9:i);
        range2 = task2_2(i-9:i);
        count = count + 1;
        moving_mean1(count) = mean(range1);
        moving_std1(count) = std(range1); 
        moving_mean2(count) = mean(range2);
        moving_std2(count) = std(range2);
    end
end


figure(3)
hold on;
%raw obs plots
plot(task2_1);
plot(task2_2);
title('Task 2.1 and Task 2.2');
legend('Timeseries 1', 'Timeseries 2');
hold off;

%moving mean plot:
figure(4)
hold on;
plot(moving_mean1);
plot(moving_mean2);
title('Moving mean for strides of 10');
legend('Stepwise mean of observations - Task2.1', 'Stepwise mean of observations - Task2.2');
hold off;

figure(5)
hold on;
plot(moving_std1);
plot(moving_std2);
title('Moving standard deviations for strides of 10');
legend('Stepwise mean of observations - Task2.1', 'Stepwise mean of observations - Task2.2');
hold off;


%% Task 3

acce = load("Exercise1-4/Accelerations_seconds.txt");
neigh = load("Exercise1-4/AlternatingNeighbors.txt");
sim = load("Exercise1-4/SimilarNeighbors.txt");
temp = load("Exercise1-4/Temperature_every5min.txt");

figure(6)


subplot(2,2,1)
plot(acce(:,1),acce(:,2))
title('Acceleration [s]');
xlabel('Time [s]')
ylabel('Data')
grid on

subplot(2,2,2)
plot(neigh(:,1),neigh(:,2))
title('Alternating neighbors');
xlabel('Time [s]')
ylabel('Data')
grid on

subplot(2,2,3)
plot(sim(:,1),sim(:,2))
title('Similar Neighbors')
xlabel('Time [s]')
ylabel('Data')
grid on

subplot(2,2,4)
plot(temp(:,1),temp(:,2))
title('Temperature [1/5min]')
xlabel('Time [s]')
ylabel('Data')
grid on


axis auto;
