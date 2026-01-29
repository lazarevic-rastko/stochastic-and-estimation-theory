function Pxx = sgs3(f,a,sigma)
    Pxx = sigma*abs(1./(a(1)+a(2)*exp(-1j*2*pi*f)+a(3)*exp(-1j*4*pi*f)+a(4)*exp(-1j*6*pi*f))).^2;
end 