function Pwelch=welchp(x,L,S,f,w)
% Welch metod
% x - odbirci ulaznog signala
% L - dužina segmenta
% S - preklapanje između segmenata
% f - frekvencijska osa
% w - prozorska funkcija dužine L

if length(w)~=L
    disp("Dužina prozorske funkcije ne odgovara traženoj dužini segmenta!");
    return
end

M=length(f);
N=length(x);

segments=my_buffer(x,L,S);
K=size(segments);
K=K(1);
Pwelch=zeros(1,M);
w=w';
for i=1:K

    segment=segments(i,:).*w;

    Pper=per(segment,f);
    Pwelch=Pwelch+Pper;
end

Pwelch=Pwelch/K;

end