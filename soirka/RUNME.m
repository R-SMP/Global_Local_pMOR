clear; close all

load('beam')
K = A{1};
G = A{2};
M = A{3};

s0 = 1i*linspace(1,4,5);
s0 = [s0 conj(s0)];
[Kr, Gr, Mr, Br, Cr, snew] = SO_IRKA_SISO(K, G, M, B, C, s0, 'os', 1e-5, 20);

res_r = zeros(1,length(s));
for ii=1:length(s)
    res_r(ii) = Cr * ((s(ii)^2*Mr + s(ii)*Gr + Kr) \ Br);
end

loglog(abs(s),abs(res),'DisplayName','original')
hold on
plot(abs(s),abs(res_r),'--','DisplayName','reduced')
plot(imag(snew),1e-4*ones(length(snew),1),'d','DisplayName','expansion points')
hold off
legend
