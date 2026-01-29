function segments=my_buffer(vec, M, S)
% M - dužina podsegmenta
% S - preklapanje izmedju susednih segmenata

s=0;
segments=[];
if S>M
    error("Duzina preklapanja, veca od duzine prozora")
end
if M>length(vec)
    error("Duzina prozora, veca od duzine signala")
end
if M<=0
    error("Duzina prozora nije validna!");
end

if S<0
    error("Duzina preklapanja medju segmentima nije validna!");
end

%važno zbog ideje, na koji je funkcija realizovana
if S==0
    S=M;
else
    S=M-S;
end

while true
    row=zeros(1,M);
    k=1;
    for i=s+1:s+M
        if i>length(vec)
            break
        end
        row(k)=vec(i);
        k=k+1;
    end
    
segments=[segments;row];

if s+M>=length(vec)
    break
end
s=s+S;

end







