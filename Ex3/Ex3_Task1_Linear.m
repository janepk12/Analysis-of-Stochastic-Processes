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
data = load("Exercise3-1.txt");
t = data(:,1);
y = data(:,2);

figure;
plot(t,y,'.')
title('Task 3.1: Acceleration measurements')
xlabel('Time [s]')
ylabel('Acceleration [g]')
hold off

%Vector of observations
L = y;

%Number of observations
no_n = length(L);

%Initial values for the unknowns
a = 1; %
b = 1; %


%Vector of initial values for the unknowns
X_0 = [a; b];

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
     L_0 = a*t + b;
     
     %Vector of reduced observations
     l = L - L_0;
    
     %Design matrix with the elements from the Jacobian matrix J
     
     A = [ t , ones(no_n, 1)];  
    
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
    
     %Check 1
     max_x_hat = max(abs(x_hat));
     
     %Vector of residuals
     v = A * x_hat -l;
 
     %Vector of adjusted observations
     L_hat = L + v;
    
     %Objective function
     vTPv = v' * P * v;
    
     %Functional relationships without the observations
     Psi = a * t+b;

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
ylabel('Acceleration [g]')
hold on 
txt = sprintf('%.2e x + %.3f', a, b)

h_line = plot(t, a*t+b, 'DisplayName', txt, 'Color', 'r', 'LineWidth', 2);
legend(h_line, 'Location', 'best');
hold off 

%--------------------------------------------------------------------------
%  Detrend the model
%--------------------------------------------------------------------------
figure;
y_new = y - Psi;
plot(t, y_new, '.', Markersize=15)
title('Task 3.1: Acceleration measurements - Detrended')
xlabel('Time [s]')
ylabel('Acceleration [g]')
hold on 


