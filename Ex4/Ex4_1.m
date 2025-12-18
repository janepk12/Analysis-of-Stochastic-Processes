%% TASK 4-1
clc; clear all; close all;
data2 = load('Exercise4-1.txt'); 
t = data2(:,1);
y = data2(:,2);

tau = [5,10,15,20,25]';
residuals =zeros(length(t),length(tau));
figure;
plot(t, y, '-', 'Color', [0.7 0.7 0.7], 'DisplayName', 'Original Observations');
hold on;

for i = 1:length(tau) 

    current_smooth = smooth_MA(y, tau(i));       
    legend_label = sprintf('MA filter size: %d', tau(i));    
    plot(t, current_smooth, '-', 'LineWidth', 1.5, 'DisplayName', legend_label);
    residuals(:,i) = y-current_smooth;
    title('Moving average ')
end
grid on;
legend show; 
hold off;

% 


figure;
%plot(t, y, '-', 'Color', [0.7 0.7 0.7], 'DisplayName', 'Original Observations');
hold on;
for i = 1:length(tau) 
    current_res = residuals(:,i);
    error = std(current_res);
    legend_label = sprintf('residual plot: %.i, std: %.2f', tau(i),error);    
    plot(t, current_res, '-', 'LineWidth', 1.5, 'DisplayName', legend_label);
    title('Residual Plot ')
    sprintf('Computed for iteration: %.1i',i)
end
grid on;
legend show; 
hold off;



%% Task 2 
% A moving average filter with window length 𝜏 = 25.
% o A Savitzky-Golay filter with window length 𝜏 = 25 and polynomial degree 𝑝 = 2.
% o A moving least squares filter with window length 𝜏 = 25, a polynomial degree 𝑝 = 2
% and Gaussian weighting.
% • Plot the resulting filtered series and interpret and discuss the results
%
clc; close all; clear all;

data1= load("Exercise4-2_1.txt");
data2= load("Exercise4-2_2.txt");
data3= load("Exercise4-2_3.txt");
data4= load("Exercise4-2_4.txt");


% RAW PLOT
figure;
subplot(2,2,1)
plot(data1(:,1), data1(:,2), '.-', 'MarkerSize', 6); grid on;
xlabel('t'); ylabel('x_1(t)');
title('Exercise4-2-1: Observation plot');


subplot(2,2,2)
plot(data2(:,1), data2(:,2), '.-', 'MarkerSize', 6); grid on;
xlabel('t'); ylabel('x_2(t)');
title('Exercise4-2-2: Observation plot');
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')

% xlim(padlim(data2(:,1))); ylim(padlim(data2(:,2)));

subplot(2,2,3)
plot(data3(:,1), data3(:,2), '.-', 'MarkerSize', 6); grid on;
xlabel('t'); ylabel('x_3(t)');
title('Exercise4-2-3: Observation plot');
% xlim(padlim(data2(:,1))); ylim(padlim(data2(:,2)));

subplot(2,2,4)
plot(data4(:,1), data4(:,2), '.-', 'MarkerSize', 6); grid on;
xlabel('t'); ylabel('x_4(t)');
title('Exercise4-2-4: Observation plot');
% xlim(padlim(data2(:,1))); ylim(padlim(data2(:,2)));



% MOVING AVERAGE FILTER
figure;
title('MOVING AVERAGE FILTER, size tau=25')
subplot(2,2,1)
plot(data1(:,1), smooth_MA(data1(:,2),25), '.-', 'MarkerSize', 6); grid on;
xlabel('t'); ylabel('x_1(t)');
title('Exercise4-2-1: Moving Average plot');


subplot(2,2,2)
plot(data2(:,1), smooth_MA(data2(:,2),25), '.-', 'MarkerSize', 6); grid on;
xlabel('t'); ylabel('x_2(t)');
title('Exercise4-2-2: Moving Average plot');
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')

% xlim(padlim(data2(:,1))); ylim(padlim(data2(:,2)));

subplot(2,2,3)
plot(data3(:,1), smooth_MA(data3(:,2),25), '.-', 'MarkerSize', 6); grid on;
xlabel('t'); ylabel('x_3(t)');
title('Exercise4-2-3: Moving Average plot');
% xlim(padlim(data2(:,1))); ylim(padlim(data2(:,2)));

subplot(2,2,4)
plot(data4(:,1), smooth_MA(data4(:,2),25), '.-', 'MarkerSize', 6); grid on;
xlabel('t'); ylabel('x_4(t)');
title('Exercise4-2-4: Moving Average plot');
% xlim(padlim(data2(:,1))); ylim(padlim(data2(:,2)));



% Savitzky Golay FILTER
figure;
subplot(2,2,1)
plot(data1(:,1), smooth_SG(data1(:,2),25,2), '.-', 'MarkerSize', 6); grid on;
xlabel('t'); ylabel('x_1(t)');
title('Exercise4-2-1: Savitzky Golay plot');


subplot(2,2,2)
plot(data2(:,1), smooth_SG(data2(:,2),25,2), '.-', 'MarkerSize', 6); grid on;
xlabel('t'); ylabel('x_2(t)');
title('Exercise4-2-2: Savitzky Golay plot');
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')

% xlim(padlim(data2(:,1))); ylim(padlim(data2(:,2)));

subplot(2,2,3)
plot(data3(:,1), smooth_SG(data3(:,2),25,2), '.-', 'MarkerSize', 6); grid on;
xlabel('t'); ylabel('x_3(t)');
title('Exercise4-2-3: Savitzky Golay plot');
% xlim(padlim(data2(:,1))); ylim(padlim(data2(:,2)));

subplot(2,2,4)
plot(data4(:,1), smooth_SG(data4(:,2),25,2), '.-', 'MarkerSize', 6); grid on;
xlabel('t'); ylabel('x_4(t)');
title('Exercise4-2-4: Savitzky Golay plot');
% xlim(padlim(data2(:,1))); ylim(padlim(data2(:,2)));



% Moving Least Squares FILTER
figure;
subplot(2,2,1)
plot(data1(:,1), smooth_MLS(data1(:,1), data1(:,2),25,2), '.-', 'MarkerSize', 6); grid on;
xlabel('t'); ylabel('x_1(t)');
title('Exercise4-2-1: Moving Least Squares plot');


subplot(2,2,2)
plot(data2(:,1), smooth_MLS(data2(:,1), data2(:,2),25,2), '.-', 'MarkerSize', 6); grid on;
xlabel('t'); ylabel('x_2(t)');
title('Exercise4-2-2: Moving Least Squares plot');
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')

% xlim(padlim(data2(:,1))); ylim(padlim(data2(:,2)));

subplot(2,2,3)
plot(data3(:,1), smooth_MLS(data3(:,1), data3(:,2),25,2), '.-', 'MarkerSize', 6); grid on;
xlabel('t'); ylabel('x_3(t)');
title('Exercise4-2-3: Moving Least Squares plot');
% xlim(padlim(data2(:,1))); ylim(padlim(data2(:,2)));

subplot(2,2,4)
plot(data4(:,1), smooth_MLS(data4(:,1), data4(:,2),25,2), '.-', 'MarkerSize', 6); grid on;
xlabel('t'); ylabel('x_4(t)');
title('Exercise4-2-4: Moving Least Squares plot');
% xlim(padlim(data2(:,1))); ylim(padlim(data2(:,2)));



%% Task 3

% x(t) = 1/2*t^2 + sin(6*pi*t) + N(0,0.1)

clc; clear all; close all;

task3 = load("Exercise4-3.txt");
t3 = task3(:,1);
x3 = task3(:,2);




figure;
plot(t3,x3,'.-','Color','r')



