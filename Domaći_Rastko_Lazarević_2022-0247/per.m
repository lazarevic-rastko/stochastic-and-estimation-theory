function Pper=per(x,f)
% Periodogram
% x - odbirci ulaznog signala
% f - frekvencijska osa

M=length(f);
N=length(x);

Pper=zeros(1,M);

for k=1:M
    for n=1:N
        Pper(k)=Pper(k)+x(n)*exp(-1j*2*pi*f(k)*(n-1));
    end
    Pper(k)=abs(Pper(k))^2/N;
end

end

 




