clear all
close all

%Add SOIRKA to path
addpath('soirka')

%Create beam model
[M, D, K, f, C] = fem_beam(1,10);

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
xlim([0 10000])
ylim([1e-12 100])
legend('Reduced system error (IRKA)',...
    'Location','northeast')
ylabel('Relative approximation error')
xlabel('Frequency (rad/s)')

% fprintf("Error norm for IRKA: %.3e\n", ...
%     error_norm_irka);

%% Wth

tol = 1e-10;
maxiter = 10;
s0 = 1i*linspace(0,10000,20);

[M, D, K, f, C] = fem_beam(1,10);


[Kr, Gr, Mr, Br, Cr, s, Vr] = SO_IRKA_SISO(K, D, M, f, C, s0, 'os', tol, maxiter)


%% Global pMOR

% Choose and initaliaze parameter space
%Parameter is length, switch to height later

P = [0.214 0.5 0.75 1 1.25]; 
%For some reason anything below 0.214 as a length does not allow the 
% SO-IRKA algorithm to function properly.

tol = 1e-10;
maxiter = 10;
s0 = 1i*linspace(0,10000,5);

%Loop over the parameter space
k = length(P);

V = []; %Not sure if this is effective or if preallocating would be better

for i = 1:k
    
    %Create model matrices for that parameter
    [M, D, K, f, C] = fem_beam(P(i),10);

    %Obtain basis matrix for that parameter
    [~, ~, ~, ~, ~, ~, Vr] = SO_IRKA_SISO(K, D, M, f, C, s0, 'os', tol, maxiter);

    %Concatenate
    V = [V, Vr];
        
end

%Perform SVD to remove rank deficient components
[U, S, Vsd] = svd(V);

for i = 1:size(S,1)
    
    if S(i,i) < 1e-10

        cull = i;
        break

    end

end

S(S<1e-10) = 0;

V = U*S*Vsd';

%Final global basis %Not sure how I should be gettin it from here

%V(:,i:end) = [];

V = (U(:,1:i));

%[U, S, Vsd] = svd(V);

%Sample use case of global basis

%Compute reduced system matrices for that certain parameter value

[M, D, K, f, C] = fem_beam(1,10);

Mr = V' * M * V;
Dr = V' * D * V;
Kr = V' * K * V;
fr = V' * f;
Cr = C * V;

s = 1i*logspace(0,4,100); %Some trouble somewhere

result_global = zeros(1,length(s));
for ii=1:length(s)
    result_irka(ii) = Cr * ((s(ii)^2*Mr + s(ii)*Dr + Kr) \ fr);
    %Careful, s should be the initial frequency range, not the eigenvalues
    %obtained from IRKA
end

%Error calculations
error_irka = (result_full - result_global)./result_full;
error_norm_irka = norm(result_full - result_global)/norm(result_full);

%Plot reponse using GLOBAL BASIS
fig = figure('Name','Frequency response using Global MOR');
set(fig,'defaulttextinterpreter','latex')
hold on
semilogy(abs(s),abs(result_global))
xlim([0 10000])
ylim([1e-8 1e-1])
title('Global MOR')
ylabel('Displacement at load point')
xlabel('Frequency (rad/s)')

% %Plot reponse using GLOBAL BASIS
% fig = figure('Name','Frequency response using IRKA');
% set(fig,'defaulttextinterpreter','latex')
% semilogy(abs(s),abs(result_full))
% hold on
% semilogy(abs(s),abs(result_irka))
% hold on
% semilogy(abs(s),abs(result_global))
% xlim([0 10000])
% ylim([1e-8 1e-1])
% legend('Full', 'SO-IRKA', 'Global MOR')
% title('Harmonic response of a clamped plate')
% ylabel('Displacement at load point')
% xlabel('Frequency (rad/s)')





%% pMOR

% Choose and initaliaze parameter space
%Parameter is length, switch to height later

P = [0.214 0.5 0.75 1 1.25]; 
%For some reason anything below 0.214 as a length does not allow the 
% SO-IRKA algorithm to function properly.

tol = 1e-10;
maxiter = 10;
s0 = 1i*linspace(0,10000,5); %r = 5 then

%Loop over the parameter space
k = length(P);

V = []; %Not sure if this is effective or if preallocating would be better

%This matrices will store the reduced system matrices for each parameter
%value
Mp = [];
Dp = [];
Kp = [];
fp = [];
Cp = [];

for i = 1:k
    
    %Create model matrices for that parameter
    [M, D, K, f, C] = fem_beam(P(i),10);

    Mp = [M, Mr];


    %Obtain basis matrix for that parameter
    [~, ~, ~, ~, ~, ~, Vr] = SO_IRKA_SISO(K, D, M, f, C, s0, 'os', tol, maxiter);

    %Concatenate
    V = [V, Vr];

    Mp = [M, Mr];
    Dp = [D, Dr];
    Kp = [K, Kr];
    fp = [f, fr];
    Cp = [C, Cr];
        
end

%Perform SVD to obtain r most significant basis vectors
[U, S, Vsd] = svd(V);

r = length(s0);

R = U(:,1:r);

%Transform reduced system matrices to generalized coordinate system R

for i = 1:k
    
    T = R'*V(:,(1:r)) %Figure out proper index distribution

end















