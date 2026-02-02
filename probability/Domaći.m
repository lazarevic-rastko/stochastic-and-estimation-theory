%% DOMAĆI - 1 SSE %%
clc
clear all
close al

%% Nepoznati parametar A
NormSuma=0;
for k=12:18
    NormSuma=NormSuma+0.6^(abs(15-k));
end
A=1/NormSuma
%% FR
clc
Out=0;
for l=12:18
    Out=0;
    for k=12:l
        Out=Out+A*0.6^(abs(15-k));
    end
    disp(Out)
end

%% FMV
k=18;
Value=A*0.6^(abs(15-k));
disp(Value)

%% Analitičko računanje očekivanja
EX=0;
for k=12:18
    EX=EX+k*A*0.6^(abs(15-k));
end
disp(EX)
%% Analitičko računanje varijanse
VarX=0;
EXX=0;
for k=12:18
    EXX=EXX+(k^2)*A*0.6^(abs(15-k));
end
disp(EXX)
VarX=EXX-EX^2;
disp(VarX)

%% Simulacija eksperimenta - detekcija jedinice
N=1e4;
rand_vec=rand(N,1);
out=zeros(N,1);

for i=1:N
    tmp=rand_vec(i);
    if tmp<0.0644
        out(i,1)=12;
    elseif tmp<0.1718
         out(i,1)=13;
    elseif tmp<0.3508
            out(i,1)=14;
    elseif tmp<0.6492
            out(i,1)=15;
    elseif tmp<0.8282
            out(i,1)=16;
    elseif tmp<0.9356
            out(i,1)=17;
    else
        out(i,1)=18;
    end
end

figure(1)
histogram(out)
title("Eksperiment sa 10 000  ishoda ")
xlabel("Broj detektovanih elektrona")
ylabel("Broj ishoda")

%% Estimacija ocekivanja i varijanse

mx=sum(out)/N;
sigmax=sum((out-mx).^2)/(N-1);

%% Grafički prikazi
k=10:0.01:20;
y=zeros(size(k));
ye=zeros(size(k));

y(k<12)=0;
y(k>=12 & k<13)=0.0644;
y(k>=13 & k<14)=0.1718;
y(k>=14 & k<15)=0.3508;
y(k>=15 & k<16)=0.6492;
y(k>=16 & k<17)=0.8282;
y(k>=17 & k<18)=0.9356;
y(18<=k)=1;
figure(1)
plot(k,y,'r')
grid("on");hold("on")
ye(k<12)=0;
ye(k>=12 & k<13)=0.0622;
ye(k>=13 & k<14)=0.1677;
ye(k>=14 & k<15)=0.3459;
ye(k>=15 & k<16)=0.6477;
ye(k>=16 & k<17)=0.8321;
ye(k>=17 & k<18)=0.9372;
ye(18<=k)=1;

plot(k,ye,'g');hold("on");
legend("Analitički oblik","Estimacija");

%% fmv
u=zeros(size(k));
u(k==12)=0.0644;
u(k==13)=0.1074;
u(k==14)=0.1790;
u(k==15)=0.2983;
u(k==16)=0.1790;
u(k==17)=0.1074;
u(k==18)=0.0644;
figure(2);
plot(k,u);grid("on");hold("on");

ue=zeros(size(k));
ue(k==12)=0.0622;
ue(k==13)=0.1055;
ue(k==14)=0.1782;
ue(k==15)=0.3018;
ue(k==16)=0.1844;
ue(k==17)=0.1051;
ue(k==18)=0.0628;
stem(k,u);grid("on");hold("on");
title("Funkcija mase verovatnoće");
legend("Analitička vrednost","Estimacija")


%% ZADATAK 2 %%
clc
clear all
N=1e5;

t=-8:0.01:6;
fy=zeros(size(t));
fy(t>=-5 & t<0)=2/13;
fy(t>=0 & t<3)=1/13;

X=rand(N,1);

Y(X>=0 & X<10/13)=13/2*X(X>=0 & X<10/13)-5;
Y(X>=10/13 & X<=1)=13*X(X>=10/13 & X<=1)-10;

figure(1);
histogram(Y,'Normalization', 'pdf');grid("on");hold("on")
plot(t,fy, 'LineWidth', 2);hold("on");grid("on")
legend("Histogram-procena fgv","Analitički oblik fgv")

fz=zeros(size(t));
fz(t>=-2 & t<=2)=1/4;

Z=-2+rand(N,1)*4;

fw=zeros(size(t));
fw(-7<t & t<=-3)=(t(-7<t & t<=-3)+7)/26;
fw(-3<t & t<=-2)=2/13;
fw(-2<t & t<=1)=(6-t(-2<t & t<=1))/52;
fw(1<t & t<=2)=(7-2*t(1<t & t<=2))/52;
fw(2<t & t<=5)=(5-t(2<t & t<=5))/52;

W=Y'+Z;
figure(2)
histogram(W,'Normalization', 'pdf');hold on;grid on;

plot(t,fw,'LineWidth',2);hold on;
legend("Histogram-procena fgv","Analitički oblik fgv")

%% zadatak 3
clc;
clear all;
N=1e5;
X_1=randn(1,N);
X_2=randn(1,N);
X=[X_1;X_2];

Ay=[1/sqrt(3) sqrt(5/3);sqrt(3) 0];
Az=[0 sqrt(2);sqrt(3-2.2*2.2/2) -2.2/sqrt(2)]

Y=Ay*X;
Z=Az*X;


figure(1)
scatter(X_1,X_2,'red'); hold on; grid on;xlabel("X1");ylabel("X2")
title("Vektor X")
figure(2)
scatter(Y(1,:),Y(2,:),'green'); hold on; grid on;xlabel("Y1");ylabel("Y2")
title("Vektor Y")
figure(3);
scatter(Z(1,:),Z(2,:),'blue'); hold on; grid on;xlabel("Z1");ylabel("Z2")
title("Vektor Z")














