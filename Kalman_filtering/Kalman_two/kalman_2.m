clc
clear all;
close all;


data = load('pad.mat');
y = data.merenja;
t = data.t;
x = data.x_stvarno;

% podaci
T = t(2); %s
g = 9.81; %m/s^2
m = 1; %kg
k = 10; %kg/s
t1 = 1; %s
t2 = 1.4; %s

process_stdev_pos = .05;
process_stdev_vel = .1;
measure_stdev1 = 0.1;
measure_stdev2 = 0.5;
sigma  = 1;

% formiranje modela
A = [1 T; 0 (1-k*T/m)];
B = [0; -g*T];
C = [1 0];
Q = (B*B')*[process_stdev_pos^2 0; 0 process_stdev_vel^2];

R = measure_stdev1^2;
%R2 = measure_stdev2^2;

% inicijalizacija
x_est = [0;0];
p_est = [1 0;0 1]*sigma;

x_estimirano = zeros(2, length(y) + 1);
p_estimirano = zeros(2, length(y) + 1);
K_zapamceno = zeros(2, length(y));

x_estimirano(:,1) = x_est;
p_estimirano = diag(p_est);

for k = 1 :  length(y)
    % predikcija
    x_pred = A*x_est + B; % B je za slučaj kada upravljacki signal postoji
    p_pred = A*p_est*A' + Q;

    % estimacija
    K = p_pred*C'*(C*p_pred*C'+R)^(-1);
    x_est = x_pred + K*(y(k) - C*x_pred);
    p_est = (eye(2) - K*C)*p_pred;

    x_estimirano(:,k+1) = x_est;
    p_estimirano(:,k+1) = diag(p_est);
    K_zapamceno(:,k) = K;

end

% prikaz rezultata

figure(1);
plot(t,y);hold on;
plot(t,x_estimirano(1,1:end-1));hold on;
plot(t,x(1,:));
legend("Merenje","Estimacija","Stvarana vrednost");


figure(2);
plot(t,x_estimirano(2,1:end-1));hold on;
plot(t,x(2,:));
legend("Estimacija","Stvarana vrednost");





