function MPper=mper(x,K,f)
% Usrednjeni periodogram
% x - odbirci ulaznog signala
% K - broj segmenata na koje usrednjavamo
% f - frekvencijska osa

N=length(x);
L=ceil(N/K);
w=rectwin(L);
MPper=welchp(x,L,0,f,w');

end