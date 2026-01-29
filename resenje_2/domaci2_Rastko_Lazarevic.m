% DRUGI DOMAĆI SAS
clc;
close all;
clear all;

data = load("podaci.mat");
x = data.x;
N = length(x);
a = data.A;
p = 4;

Az = tf(1,a,1);
poles = pole(Az);
figure
pzmap(Az);

% frekfencije od interesa
w_i =[angle(poles(1)) angle(poles(3))];

% Očekujemo 2 pika i to na učestanostima 0.7854rad/odb i 1.0472rad/odb.
% To su argumenti polova. Imamo dva para konj. kompleksnih polova.
% Dominatni pik očekujemo na učestanosti pola sa najvećim potegom
% i to je onda učestanost 0.784rad/odb


%% autokorelaciona metoda
[a_a, s_a] = autocorr_estimation(x,p);
% ugradjena funkcija
[a_aorg, s_a_org] = aryule(x, p);

%% modifikovana kovarijantna metoda
[a_mk, s_mk] = modcov_estimation(x,p);
% ugradjena funkcija
[a_mk_org, s_mk_org] = armcov(x, p);

%% poredjenje dobijenih koeficijenata 

% original
clc
disp("Original: "+num2str(a(1))+", "+num2str(a(2))+", "+num2str(a(3))+...
    ", "+num2str(a(4))+", "+num2str(a(5)));

% ACR
disp("Autokorelaciona metoda: ");
disp("Koeficijenti ACR (moja): "+num2str(a_a(1))+", "+num2str(a_a(2))+", "+num2str(a_a(3))+...
    ", "+num2str(a_a(4))+", "+num2str(a_a(5)));
disp("Varijansa ACR (moja): " + num2str(s_a));
disp("Koeficijenti ACR (org): "+num2str(a_aorg(1))+", "+num2str(a_aorg(2))+", "+num2str(a_aorg(3))+...
    ", "+num2str(a_aorg(4))+", "+num2str(a_aorg(5)));
disp("Varijansa ACR (org): " + num2str(s_a_org));


%%
 
% MODCOV
clc;
disp("Modifikovana kovarijantna metoda")
disp("Koeficijenti MC (moja): "+num2str(a_mk(1))+", "+num2str(a_mk(2))+", "+num2str(a_mk(3))+...
    ", "+num2str(a_mk(4))+", "+num2str(a_mk(5)));
disp("Varijansa MC (moja): " + num2str(s_mk));
disp("Koeficijenti MC (org): "+num2str(a_mk_org(1))+", "+num2str(a_mk_org(2))+", "+num2str(a_mk_org(3))+...
    ", "+num2str(a_mk_org(4))+", "+num2str(a_mk_org(5)));
disp("Varijansa MC (org): " + num2str(s_mk_org));


%% AIC(k)
pp = [2 3 4 5 6 7 8 9 10 11 12 14 16 18 20 22 24 26];
AIC = zeros(1,length(pp));
 
for i = 1:length(pp)
    %[~,sigma_i] = modcov_estimation(x,pp(i));
    [~,sigma_i] = autocorr_estimation(x,pp(i));
    %AIC(i) = N*log(sigma_i) + 2*pp(i);
    AIC(i) = N*log(sigma_i) + pp(i)*log(N);
end

disp("Estimirani red je 3.")
figure(1)
title("Drugi Akaikijev kriterijum");
plot(pp,AIC,'LineWidth',1.7);
grid on;
xlabel("k");ylabel("AIC(k)");


%% pordjenje razlicitih sgs
pest = 3;
f = linspace(0,0.5,1000);

% autokorelaciona
[est1, s1] = autocorr_estimation(x,pest);
[est1_org,s1_org] = aryule(x,pest);

[est2, s2] = modcov_estimation(x,pest);
[est2_org, s2_org] = armcov(x,pest);


% tacna
Pxx = 1*abs(1./(a(1)+a(2)*exp(-1j*2*pi*f)+a(3)*exp(-1j*4*pi*f)+a(4)*exp(-1j*6*pi*f)+a(5)*exp(-1j*8*pi*f))).^2;
% autokorelaciona
Pxx_e1 = sgs3(f,est1,s1);
Pxx_e1_org = sgs3(f,est1_org,s1_org);
% modifikovana kovarijantna
Pxx_e2 = sgs3(f,est2,s2);
Pxx_e2_org = sgs3(f,est2_org,s2_org);

figure(2);
title("SGS našeg procesa")
plot(f,10*log10(Pxx));hold on;

plot(f,10*log10(Pxx_e1));hold on;
plot(f,10*log10(Pxx_e1_org),'g--');hold on


plot(f,10*log10(Pxx_e2));hold on;
plot(f,10*log10(Pxx_e2_org),'b--');hold on;
title("SGS")
xlabel("f");
ylabel("Pxx(f)");
legend("Teorijska","ACR moja","ACR ugrađena","MC moja","MC ugrađena");

% cak ni u teorijskom prikazu ne vidimo manje dominantan pik
% on se prakticno utopi u prvi koji je prilicno dobro izrazen (visok i uzak)
% dok je manji širi i niži
% tako da i ima smisla što se ne primeti
% to opravdava i cinjenica o njihovoj blizini u spektru i odnos potega

