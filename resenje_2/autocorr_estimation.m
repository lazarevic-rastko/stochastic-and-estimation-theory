function [a,s] = autocorr_estimation(x, p)
% x - ulazni signal
% p - red modela
    N = length(x);
    rxx = zeros(2*p+1,1);
    
    for k = (p+1):(2*p+1)
        rxx(k) = sum(conj(x(1:(N-k+p+1))).*x(1+k-(p+1):N))/N;
    end
    rxx(p:-1:1) = conj(rxx(p+2:end));
    R = toeplitz(rxx(p+1:2*p),rxx((p+1):-1:2));
    
    a = -R\rxx(p+2:end);
    s = rxx(p+1) + sum(a.*rxx(p:-1:1));

    a = [1 a'];

end