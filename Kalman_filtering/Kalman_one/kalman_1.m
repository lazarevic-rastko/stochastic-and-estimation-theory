%% kalmanov filtar

data = load('kalman.mat');
y = data.y;
t = data.t;
x = data.x_stvarno;
% podaci
T = t(2);
A = [1 T;0 1];
B = [T^2/2;T];
C = [1 0];

sigmaV = 25;
sigmaA = 1;

R = sigmaV; % sum merenja
Q =[T^4/4 T^3/2;T^3/2 T^2]*sigmaA ; % sum modela

% inicijalizacija

sigma = 1;

% inicijalizacija
x_est = [0; 0];
P_est = [1 0;0 1]*sigma;

x_estimirano = zeros(2, length(y)+1);
p_estimirano = zeros(2, length(y)+1);
K_zapamceno = zeros(2,length(y));

x_estimirano(:,1) = x_est;
p_estimirano(:,1) = diag(P_est);


% iteracije

for k = 1:length(y)
    % predikacija
    x_pred = A*x_est;
    P_pred = A*P_est*A'+Q;

    % Kalmanovo pojačanje
    K = P_pred*C'*(C*P_pred*C'+R)^(-1);
    
    %estimacija

    x_est = x_pred + K*(y(k) - C*x_pred);
    P_est = (eye(2) - K*C)*P_pred;

    x_estimirano(:,k+1) = x_est;
    p_estimirano(:,k+1) = diag(P_est);
    K_zapamceno(:,k) = K;

end

% analiza rezultata

figure(1);
plot(t,y);hold on;
plot(t,x(1,:));hold on;
plot(t,x_estimirano(1,1:end-1));grid on;
title("Estimacija pozicije");
legend("Merenje","Stvarna vrednost","Estimacija");

figure(2);
plot(t,x(2,:));hold on;
plot(t,x_estimirano(2,1:end-1));grid on;
title("Estimacija brzine");
legend("Stvarna vrednost","Estimacija");

figure(3);
subplot(2,1,1)
plot(t,p_estimirano(1,1:end-1));
title("Varijansa estimacije pozicije");
subplot(2,1,2);
plot(t,p_estimirano(2,1:end-1));
title("Varijansa estimacije brzine");

figure(4);
plot(t,K_zapamceno);
title("Kalmanovo pojacanje");


