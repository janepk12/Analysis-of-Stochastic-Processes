clc; clear;close all;
% Exercise 3
% Task 1 
data1 = load("Exercise3-1.txt");
data2 = load("Exercise3-2.txt");
t1 = data1(:,1);
y1 = data1(:,2);

t2 = data2(:,1);
y2 = data2(:,2);

figure;
plot(t1,y1,'.')
figure;
plot(t2,y2,'.')
%axis equal

