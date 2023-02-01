clear all
close all

%Add SOIRKA to path
addpath('soirka')

%Create beam model
[M, D, K, f, C] = fem_beam(1,200);

%% Compute full solution

%I think C and o are the same or the very least they serve the same purpose,
% they get the DOF corresponding to the displacement at the tip.
% o = sparse(1,length(f));
% o(f~=0) = 1;

%s = 1i*linspace(0,10000,100);
s = 1i*logspace(0,4,100);
%Careful with logspace, it might lead to very high values if used
%incorrectly

result_full = zeros(1,length(s));

tic
for ii=1:length(s)
    ii;
    result_full(ii) = C * ((s(ii)^2*M + s(ii)*D + K) \ f);
end
ctime_full = toc;


%Plot full response
fig = figure('Name','Frequency response, full solution');
set(fig,'defaulttextinterpreter','latex')
semilogy(abs(s),abs(result_full))
xlim([0 10000])
ylim([1e-8 1e-1])
legend('Response')
title('Harmonic response of a clamped plate')
ylabel('Displacement at load point')
xlabel('Frequency (rad/s)')

%% Compute SO-IRKA solution using algorithm implementation from Exericse 3.2

%s = 1i*linspace(0,10000,100);
s = 1i*logspace(0,4,100);

%We can change the "r" here, we have chosen r = 20
s0 = 1i*linspace(0,10000,20); 
%s0 = 1i*logspace(0,4,5);

%IRKA iteration parameters
tol = 1e-10;
maxiter = 10;

tic

%2 options for the IRKA algorithm:

%1: Run IRKA algorithm given Exercise 3.2, bit less accurate than option 2
% "irka.m"

%[Mr, Dr, Kr, fr, Cr, sr] = irka(M, D, K, f, C, s0, tol, maxiter);

%2: Run SO-IRKA algorithm supplied by Sebastian, here r is doubled,
%but it takes as much time as option 2, making it a lot more accurate in
%the same amount of time.
% SO-IRKA-SISO.m"

[Kr, Dr, Mr, fr, Cr, sr] = SO_IRKA_SISO(K, D, M, f, C, s0, 'os', tol, maxiter);

ctime_irka_offline = toc;

tic
result_irka = zeros(1,length(s));
for ii=1:length(s)
    result_irka(ii) = Cr * ((s(ii)^2*Mr + s(ii)*Dr + Kr) \ fr);
    %Careful, s should be the initial frequency range, not the eigenvalues
    %obtained from IRKA
end
ctime_irka_online = toc;

%IRKA error calculations
error_irka = (result_full - result_irka)./result_full;
error_norm_irka = norm(result_full - result_irka)/norm(result_full);

%Plot reponse using IRKA
fig = figure('Name','Frequency response using IRKA');
set(fig,'defaulttextinterpreter','latex')
semilogy(abs(s),abs(result_irka))
xlim([0 10000])
ylim([1e-8 1e-1])
legend('Response')
title('Harmonic response of a clamped plate')
ylabel('Displacement at load point')
xlabel('Frequency (rad/s)')

%Plot overlaid results and error plot
fig = figure();
set(fig,'defaulttextinterpreter','latex')
subplot(2,1,1)
semilogy(abs(s),abs(result_full),...
    abs(s),abs(result_irka))
xlim([0 10000])
ylim([1e-8 1e-1])
legend('Full system', 'Reduced System (IRKA)')
title('Comparison of full and reduced models')
ylabel('Displacement at load point')
subplot(2,1,2)
semilogy(abs(s),abs(error_irka))
xlim([[0 10000]])
ylim([1e-12 100])
legend('Reduced system error (IRKA)',...
    'Location','northeast')
ylabel('Relative approximation error')
xlabel('Frequency (rad/s)')

% fprintf("Error norm for IRKA: %.3e\n", ...
%     error_norm_irka);


%% pMOR - Choose and initaliaze parameter space


%Parameter is length, switch to height later


%% Global pMOR










