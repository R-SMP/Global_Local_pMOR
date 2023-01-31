%% MOR

%Create beam model
[M, D, K, f, C] = fem_beam(1,10);

det(M)

det(D)

det(K)

o = sparse(1,length(f));
o(f~=0) = 1;
%I think C and o are the same

s = 1i*logspace(0,10000,100);
result_full = zeros(1,length(s));

C*((s(1)^2*M + s(1)*D + K) \ f)

tic
for ii=1:length(s)
    ii
    result_full(ii) = o * ((s(ii)^2*M + s(ii)*D + K) \ f);
end
ctime_full = toc
% 
% fig = figure('Name','Frequency response');
% set(fig,'defaulttextinterpreter','latex')
% semilogy(abs(s),abs(result_full))
% xlim([0 10000])
% ylim([1e-8 1e-1])
% legend('Response')
% title('Harmonic response of a clamped plate')
% ylabel('Displacement at load point')
% xlabel('Frequency (rad/s)')