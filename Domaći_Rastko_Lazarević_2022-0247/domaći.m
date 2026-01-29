clc;
close all;
clear all;
%%
% Parametar
%2022/0247
q = mod(0+2+4+7,5);

% Učitavnje podataka
data = open('x3.csv');
x3 = data.x3;

%% Prikaz signala u vremenskom domenu

t=0:1:255;

plot(t,x3);grid on;
xlabel("n[odb]");
ylabel("x3(n)");

%% Frekvencijska analiza k-te realizacije procesa

k = 10;
x3_r = x3(k,:);

N = 2^(nextpow2(length(x3))+1);
X3_r = fft(x3_r,N)/length(x3_r);
X3_r(2:N/2+1) = 2*X3_r(2:N/2+1);

X3_r = X3_r(1:N/2+1);
f_fft = 0:1/N:1/2;

figure(1);
title("AFK")
plot(f_fft,abs(X3_r)); grid on;
xlabel("f - dgitalna učestanost");
ylabel("|X(f)|");

figure(2);
title("FFK")
plot(f_fft,unwrap(angle(X3_r))); grid on;
xlabel("f - dgitalna učestanost");
ylabel("arg{X(f)}");

%% PERIODOGRAM
f = -0.5:1/499:0.5;


rand_num = [10 35 26 43 7];
figure(1)
for num = rand_num
    Pper = per(x3(num,:),f);
    plot(f,Pper);grid on;hold on
end
xlabel("f - dgitalna učestanost");
ylabel("P(f)");
title("Periodogram signala: " + num2str(rand_num(1)) + ...
    ", " + num2str(rand_num(2)) + ", " + num2str(rand_num(3)) + ...
    ", " + num2str(rand_num(4)) + ", " + num2str(rand_num(5)));

legend("x3(" + num2str(rand_num(1)) + ")", ...
       "x3(" + num2str(rand_num(2)) + ")", ...
       "x3(" + num2str(rand_num(3)) + ")", ...
       "x3(" + num2str(rand_num(4)) + ")", ...
       "x3(" + num2str(rand_num(5)) + ")");

disp("Bitne odlike našeg signala su harmonici na 0.15 i 0.30.");
disp("Sinusoide na digitalnim učestanostima 0.15 i 0.30.")



%% Blackman-Tukey zatvaranje prozora
P_015 = [];
P_044 = [];
B = [];


%eksperiment
i = 1;
for M = 1:2:255
    break
    Pbt = bt_psd(x3(35,:),f,M,'bart');
    figure(i);
    plot(f,Pbt);grid on;
    xlabel('f[Hz]');ylabel('P(f)');
    title("Blackman - Tuckey " + num2str(M));
    i = i + 1;
    
    %B=[B M];
    %P_015 = [P_015 Pbt(1301)];
    %P_044 = [P_044 Pbt(1889)];

    answ = input('',"s");
    if answ == "d"
        %sledeći
        continue
    else 
        break
    end

end
    %figure(255);
    %plot(B,P_015);grid on;hold on;
    %plot(B,P_044);
    %legend("f = 0.15","f = 0.44");
    %xlabel("M - duzina prozorske funkcije")


% skiciranje
MO = [31 73 181];

figure(2);
for M = MO
    Pbt=bt_psd(x3(25,:),f,M,'bart');
    plot(f,Pbt);grid on;hold on;
    xlabel('f[Hz]');ylabel('P(f)');
end
title("Blackman - Tukey");
legend("M = "+num2str(MO(1)),"M = "+num2str(MO(2)),"M = "+num2str(MO(3)));


disp("Optimalan odbabir dužine prozora: M = 73");

%% Welch zatvaranje prozora

% Lmax = 256
% Startna merenja, koja sam modifikovao u zavinosti od 
% dinamike promene L i S
% 1) L = 65     S = 9
% 2) L = 65     S = 36
% 3) L = 65     S = 52
% ----------------------
% 4) L = 152    S = 22
% 5) L = 152    S = 83
% 6) L = 152    S = 122
% ----------------------
% 7) L = 205    S = 30
% 8) L = 205    S = 111
% 9) L = 205    S = 166

%53 32
LS=[120 22;120 83;120 100;65 9;58 36;65 55;205 30;205 111;205 166];

figure(3);
for i = 1:length(LS)
    L = LS(i,1);
    S= LS(i,2);

    w=hamming(L);
    Pwlc=welchp(x3(25,:),L,S,f,w);
    subplot(3,3,i)
    plot(f,Pwlc);grid on;
    xlabel("f");
    ylabel("P(f)");
    title("L = "+num2str(L)+" i S =  "+num2str(S));
end

% dodatna merenja
% L = 50    S = 5
% L = 50    S = 24
% L = 50    S = 38

LS=[50 5;50 24;50 38;65 36];
figure(10)
for i = 1:length(LS)
    L = LS(i,1);
    S= LS(i,2);

    w=hamming(L);
    Pwlc=welchp(x3(25,:),L,S,f,w);
    subplot(2,2,i)
    plot(f,Pwlc);grid on;
    xlabel("f");
    ylabel("P(f)");
    title("L = "+num2str(L)+" i S =  "+num2str(S));
end



disp("Odabirano: L = 65 i S = 36.");


%% Računanje varijanse Periodogram
Num_r = 50;
EX_per = zeros(1,length(f));
out_per = zeros(Num_r,length(f));

for i = 1 : Num_r
    out_per(i,:) = per(x3(i,:),f);
    EX_per = EX_per + out_per(i,:);
end
EX_per = EX_per/Num_r;

VarX_per = zeros(1,length(f));
for i = 1:length(f)
    VarX_per(i) = sum((out_per(:,i)-EX_per(i)).^2);
end

VarX_per = VarX_per/(Num_r-1);

%medianPer_r = median((out_per(:)-sum(out_per(:))/50/length(f)).^2);
medianPer_f = median(VarX_per);

%% Računanje varijanse Welch

w=hamming(65);

EX_welch = zeros(1,length(f));
out_welch = zeros(Num_r,length(f));

for i = 1:Num_r
    out_welch(i,:) = welchp(x3(i,:),65,36,f,w);
    EX_welch = EX_welch + out_welch(i,:);
end
EX_welch = EX_welch/Num_r;


VarX_welch = zeros(1,length(f));
for i = 1:length(f)
    VarX_welch(i) = sum((out_welch(:,i)-EX_welch(i)).^2);
end

VarX_welch = VarX_welch/(Num_r-1);

medianWelch_f = median(VarX_welch);

%% Računanje varijanse Blackman - Tukey

EX_bt = zeros(1,length(f));
out_bt = zeros(Num_r,length(f));

for i=1:Num_r
    out_bt(i,:) = bt_psd(x3(i,:),f,73,'bart');
    EX_bt= EX_bt + out_bt(i,:);
end
EX_bt = EX_bt/Num_r;

VarX_bt = zeros(1,length(f));
for i=1:length(f)
    VarX_bt(i) = sum((out_bt(:,i) - EX_bt(i)).^2);
end

VarX_bt = VarX_bt/(Num_r-1);

medianBT_f = median(VarX_bt);

%% poređenje varijanse

figure(5);
plot(f,VarX_per);hold on;
plot(f,VarX_bt);hold on;
plot(f,VarX_welch);hold on;grid on;
xlabel('f')
title("Varijansa");
legend("Periodogram","Blackman-Tuckey","Welch")
legend("Blackman-Tuckey","Welch")

%% Poređenje rezultata i dobijene varijanse Periodograma
var_analitical = zeros(1,length(f));

for i = 1:length(f)
    var_analitical(i) = EX_per(i)^2*(1+(sin(2*pi*256*f(i))/sin(2*pi*f(i))/256)^2);
end

figure(6);
plot(f,VarX_per);grid on;hold on;
plot(f,var_analitical);
legend("Numerički","Analitički");
title("Analitička i numerička vrednost varijanse periodograma")
xlabel("f");

figure(11)
plot(f,VarX_per);grid on;hold on;
title("Numerička vrednost varijanse periodograma")



%% Intervali poverenja Blackman - Tukey

M = 73;
N = 256;
bartlett(M);

ni = ceil(2*N/sum(bartlett(M).^2));

alpha = 0.05;

q_L = icdf('Chisquare', 1 - alpha/2, ni);
q_R = icdf('Chisquare', alpha/2, ni);

L = ni*out_bt(25,:)/q_L;
R = ni*out_bt(25,:)/q_R;

figure(7);
title("Interval poverenja");
plot(f,10*log(out_bt(25,:)));hold on;
plot(f,10*log(L));hold on;
plot(f,10*log(R));hold on;
plot(f,10*log(EX_per));
fill([f, fliplr(f)], [10*log(L), fliplr(10*log(R))], 'r', ...
    'FaceAlpha', 0.1, 'EdgeColor', 'none');
xlabel('f');
ylabel("Pf_db");
legend("BT estimacija", "Donja granica","Gornja granica","Aproksimacija Pxx(f)");
title("Estimacija Blackman - Tukey estimatora i interval poverenja");
 
%% Simulacija sa kraćim sekvencama signala

EX_per_2 = zeros(1,length(f));
out_per_2 = zeros(Num_r,length(f));

for i=1:Num_r
    out_per_2(i,:) = per(x3(i,1:64),f);
    EX_per_2 = EX_per_2 + out_per_2(i,:);
end
EX_per_2 = EX_per_2/Num_r;

VarX_per_2 = zeros(1,length(f));
for i=1:length(f)
    VarX_per_2(i) = sum((out_per_2(:,i)-EX_per_2(i)).^2);
end

VarX_per_2 = VarX_per_2/Num_r;

medianPer_f_2 = median(VarX_per_2);

figure(8);
plot(f,VarX_per);grid on;hold on;
plot(f,VarX_per_2);
xlabel("f");
ylabel("Pf")
title("Varijansa s.p. za različit broj odbiraka");
legend("256 - odbiraka","64 - odbirka");







%% Testiranje funkcija 1
Test = ispisiRedoveFajla('P01_test_sekvenca_1.csv');

%% Testiranje funckija 2

%[k,r]=akf(Test);

%f=-0.5:1/50:0.5;
%Pbt=per(Test,f);
%figure(11);

%plot(f,Pbt);







