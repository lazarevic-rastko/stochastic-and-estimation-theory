function [a,s] = modcov_estimation(x, p)
    N = length(x);
    cxx = zeros(p+1,p+1);

    for j = 0:p
        for k = 0:p
            cxx(j+1,k+1) = (sum(conj(x(p-j+1:N-j)) .* x(p-k+1:N-k)) ...
                     + sum(x(1+j:N-p+j) .* conj(x(1+k:N-p+k))) ) / (2*(N-p));
        end
    end

    a = -cxx(2:end,2:end)\cxx(2:end,1);
    a = a';
    s = cxx(1,1) + sum(a.*cxx(1,2:p+1));
    a = [1 a];

end