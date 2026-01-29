clc;
clear all;
close all;

data = load('vozGPS.mat');
y = data.y;

%%
T = 1;
u = zeros(1, length(y));
u(1:30) = 1; % upravljacka sekvenca

A = [1 T;0 T];
C = [1 0];
B = [T^2/2;T];

R = 50^2;

Q1 = [T^4/4 T^3/2; T^3/2 T^2];
u1 = zeros(1, length(y));
u1(31:end) = 1;
%u1 = 1 + u1;


% inicijalizacija

x_est = [0;0];
P_est = zeros(2,2);

% vektori u kojim pamtimo

x_estimirano = zeros(2, length(y)+1);
P_estimirano = zeros(2, length(y)+1);
K_zapamceno = zeros(2,length(y));

x_estimirano(:,1) = x_est;
P_estimirano(:,1) = diag(P_est);

% iteracije

for k = 1:length(y)
    % predikcija
    x_pred = A*x_est + B*u(k);

    Q = Q1*u1(k);
    P_pred = A*P_est*A' + Q;

    %estimacija
    K = P_pred*C'*(C*P_pred*C'+R)^(-1);
    x_est = x_pred + K*(y(k)-C*x_pred);
    P_est = (eye(2,2) - K*C)*P_pred;

    % pamcenje vrednosti
    x_estimirano(:,k+1) = x_est;
    P_estimirano(:,k+1) = diag(P_est);
    K_zapamceno(:,k) = K;
    
end

% prikaz rezultata
n = 0:T:length(y)-T;
figure(1);
plot(n,y);hold on;
plot(n,x_estimirano(1,1:end-1));grid on;
title("Kalmanov filtar za pracenje pozicije voza");
legend("Merenja","Estimacija");
