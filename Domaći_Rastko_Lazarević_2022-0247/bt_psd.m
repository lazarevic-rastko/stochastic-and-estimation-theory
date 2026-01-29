function Psdbt=bt_psd(x,f,M,wind)
% x - odbirci ulaznog signala
% f - frekvencijska osa
% M - dužina prozorske funkcije, M<N i M neparno
% wind - tip prozorske funkcije

if mod(M,2)==0
    disp("Dužina prozorske funkcije je neparan broj");
    return
end

N=length(x);
L=length(f);

[~,Rxx]=akf(x);

w=zeros(1,2*N-1);
switch wind
    case 'bart'
        w1=bartlett(M);
    case 'parzen'
        w1=parzenwin(M);
    otherwise
        w1=bartlett(M);
end

for k=N-(M-1)/2:N+(M-1)/2
    w(k)=w1(k-N+(M-1)/2+1);
end

Psdbt=zeros(1,L);
for i=1:L
    for k=1:2*N-1
        Psdbt(i)=Psdbt(i) + w(k)*Rxx(k)*exp(-1j*2*pi*f(i)*(k-N));
    end
end

end








