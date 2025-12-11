%--------------------------------------------------------------------------
%   
%          ADJUSTMENT THEORY I
%    Exercise 9: Adjustment Calculation - part IV  
% 
%   Author         : Anastasia Pasioti
%   Version        : October 11, 2018
%   Last changes   : January 11, 2023
%
%--------------------------------------------------------------------------

clc;
clear all;
close all;

%--------------------------------------------------------------------------
%   Task 1
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   Observations and initial values for the unknowns
%--------------------------------------------------------------------------
%Load data


%Error-free values
data = load("Exercise3-3.txt");
t = data(:,1);
y = data(:,2);

figure;
plot(t,y,'.')
title('Task 3.1: Raw measurements')
xlabel('Time [s]')
ylabel('Observation [-]')
hold off

%Vector of observations
L = y;

%Number of observations
no_n = length(L);

%Initial values for the unknowns
a = 1; %
b = 1; %
c= 1;
d = 1;
e=1;
f = 1;

%Vector of initial values for the unknowns
X_0 = [a; b;c;d;e;f ];

%Number of unknowns
no_u = length(X_0);

%Redundancy
r = no_n - no_u;

%--------------------------------------------------------------------------
%  Stochastic model
%--------------------------------------------------------------------------

Q_LL = eye(no_n);

%Weight matrix
P = inv(Q_LL);

%--------------------------------------------------------------------------
%  Adjustment
%--------------------------------------------------------------------------
%break-off conditions
epsilon = 10^-12;
delta = 10^-12;
max_x_hat = Inf;
Check2 = Inf;

%Number of iterations
iteration = 0;

while max_x_hat > epsilon || Check2 > delta         
    
     %Observations as functions of the approximations for the unknowns
     %L_0 = a*t.^3 + b*t.^2 + c*t+d;
     %L_0 = a*t.^4 + b*t.^3 + c*t.^2+d*t+e;
     L_0 = a*t.^5 + b*t.^4 + c*t.^3+d*t.^2+e*t+f;
     
     %Vector of reduced observations
     l = L - L_0;
    
     %Design matrix with the elements from the Jacobian matrix J
     
     A = [t.^5, t.^4 , t.^3 , t.^2,t , ones(no_n, 1)];  
    
     %Normal matrix
     N = A' * P * A;
     
     %Vector of right hand side of normal equations
     n = A' * P * l;
    
     %Inversion of normal matrix / Cofactor matrix of the unknowns
     Q_xx = inv(N);
    
     %Solution of the normal equations
     x_hat = Q_xx * n;
      
     %Update
     X_hat = X_0 + x_hat;
     X_0 = X_hat;

     a = X_hat(1);
     b = X_hat(2);
     c = X_hat(3);
     d = X_hat(4);
     e = X_hat(5);
     f = X_hat(6);
    
     %Check 1
     max_x_hat = max(abs(x_hat));
     
     %Vector of residuals
     v = A * x_hat -l;
 
     %Vector of adjusted observations
     L_hat = L + v;
    
     %Objective function
     vTPv = v' * P * v;
    
     %Functional relationships without the observations
     %Psi = a*t.^3 +b*t.^2 + c*t + d;
     %Psi = a*t.^4 + b*t.^3 + c*t.^2+d*t+e;
     Psi = a*t.^5 + b*t.^4 + c*t.^3+d*t.^2+e*t + f;
     %Check 2
     Check2 = max(abs(L_hat- Psi)); 
    
     %Update number of iterations
     iteration = iteration+1;
  
end

if Check2<=delta
    disp('Everything is fine!')
else
    disp('Something is wrong.')
end

iteration

%--------------------------------------------------------------------------
%  Plot
%--------------------------------------------------------------------------
figure;
plot(t, y, '.', Markersize=15)
title('Task 3.1: Regressional computation')
xlabel('Time [s]')
ylabel('Observations [-]')
hold on 
txt = sprintf('%.2e x^5 + %.3f x^4 + %.3f x^3 +  %.3f x^2 +  %.3f x +  %.3f x ', a, b,c,d,e,f)

h_line = plot(t, a*t.^5+b*t.^4+c*t.^3+d*t.^2+e*t+f, 'DisplayName', txt, 'Color', 'r', 'LineWidth', 2);
legend(h_line, 'Location', 'best','FontSize',15);
hold off 

%--------------------------------------------------------------------------
%  Detrend the model
%--------------------------------------------------------------------------
figure;
y_new = y - Psi;
data_to_save = [t(:), y_new(:)];
writematrix(data_to_save, 'task3_1_data.txt', 'Delimiter', '\t'); % write to matrix to compute xcross

plot(t, y_new, '.', Markersize=15)
title('Task 3.1: Observation measurements - Detrended')
xlabel('Time [s]')
ylabel('Observations [-]')
hold on 


