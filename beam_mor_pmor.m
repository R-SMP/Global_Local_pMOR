clear
close all

%Add SOIRKA to path
addpath('soirka')
%% Inputs
numSamples = 1000; % number of frequency samples

%Create beam model
[M, D, K, f, C] = fem_beam(1,10);

%% Compute full solution

% C and o are the same and serve the same purpose,
% they get the DOF corresponding to the z-displacement at the tip.
% o = sparse(1,length(f));
% o(f~=0) = 1;

s = 1i*logspace(0,4,numSamples); % frequency sample values
%Careful with logspace, it might lead to very high values if used
%incorrectly

result_full = zeros(1,length(s));

tic
for ii=1:length(s)
    result_full(ii) = C * ((s(ii)^2*M + s(ii)*D + K) \ f);
end
ctime_full = toc;
fprintf("full:" + ctime_full);

%Plot full response
fig = figure('Name','Frequency response, full solution');
set(fig,'defaulttextinterpreter','latex')
semilogy(abs(s),abs(result_full))
xlim([0 10000])
ylim([1e-8 1e-1])
legend('Response')
title('Frequency response, full solution')
ylabel('Displacement at load point')
xlabel('Frequency (rad/s)')

%% Compute SO-IRKA solution using algorithm implementation from Exericse 3.2

s = 1i*logspace(0,4,numSamples);

%We can change the "r" here, we have chosen r = 20
r = 20;
s0 = 1i*linspace(0,10000,r); 
%s0 = 1i*logspace(0,4,5);

%IRKA iteration parameters
tol = 1e-10; % tolerance
maxiter = 10; % number of maximum iterations

tic

%2 options for the IRKA algorithm:

%1: Run IRKA algorithm given Exercise 3.2, bit less accurate than option 2
% "irka.m"

%[Mr, Dr, Kr, fr, Cr, sr] = irka(M, D, K, f, C, s0, tol, maxiter);

%2: Run SO-IRKA algorithm supplied by Sebastian, here r is doubled,
%but it takes as much time as option 2, making it a lot more accurate in
%the same amount of time.
% SO-IRKA-SISO.m"

% generate reduced-order beam model
[Kr, Dr, Mr, fr, Cr, sr] = SO_IRKA_SISO(K, D, M, f, C, s0, 'os', tol, maxiter);

ctime_irka_offline = toc;
fprintf("offline:" + ctime_irka_offline);

tic
result_irka = zeros(1,length(s));
for ii=1:length(s)
    result_irka(ii) = Cr * ((s(ii)^2*Mr + s(ii)*Dr + Kr) \ fr);
    %Careful, s should be the initial frequency range, not the eigenvalues
    %obtained from IRKA
end
ctime_irka_online = toc;
fprintf("irka:" + ctime_irka_online);

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
title('Frequency response using IRKA')
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


%% Global pMOR

% Choose and initaliaze parameter space
% Parameter is length, switch to height later

%P = [0.214 0.5 0.75 0.9 1.25]; %original
P = [0.214 0.3 0.4 0.5 0.6 0.75 0.9 1.25]; 
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

    %Concatenate reduced bases
    V = [V, Vr];
        
end

%Perform SVD to remove rank deficient components
[U, S, Vsd] = svd(V);

for i = 1:size(S,1)
    
    if S(i,i) < 1e-10

        redIndex = i; % index beyond which eigenvalues are small enough to be neglected 
        % equivalent to reduced dimension of reduced global basis
        break

    end

end
U = U(:, 1:redIndex);
S = S(1:redIndex, 1:redIndex);
Vsd = Vsd(1:redIndex, 1:redIndex);

V = U*S*Vsd'; % reduced global bases 

%Sample use case of global basis

%Compute reduced system matrices for that certain parameter value, 1
%chosen here

[M, D, K, f, C] = fem_beam(1,10);

Mr = V' * M * V;
Dr = V' * D * V;
Kr = V' * K * V;
fr = V' * f;
Cr = C * V;

s = 1i*logspace(0,4,numSamples); 

result_global = zeros(1,length(s));
for ii=1:length(s)
    lhs = s(ii)^2*Mr + s(ii)*Dr + Kr;
    result_global(ii) = Cr * (lhs \ fr);
    %Careful, s should be the initial frequency range, not the eigenvalues
    %obtained from IRKA
end

%Error calculations
% error_irka = (result_full - result_global)./result_full;
% error_norm_irka = norm(result_full - result_global)/norm(result_full);

%Plot reponse using GLOBAL BASIS
fig = figure('Name','Frequency response using Global MOR');
set(fig,'defaulttextinterpreter','latex')
semilogy(abs(s),abs(result_global))
xlim([0 10000])
ylim([1e-8 1e-1])
title('Global MOR')
ylabel('Displacement at load point')
xlabel('Frequency (rad/s)')

%Plot reponse using GLOBAL BASIS
fig = figure('Name','Frequency response comparison up to global pMOR');
set(fig,'defaulttextinterpreter','latex')
semilogy(abs(s),abs(result_full), 'LineWidth', 6)
hold on
semilogy(abs(s),abs(result_irka), 'LineWidth', 4)
hold on
semilogy(abs(s),abs(result_global), 'LineWidth', 2)
xlim([0 10000])
ylim([1e-8 1e-1])
legend('Full', 'SO-IRKA', 'Global MOR')
title('Frequency response comparison up to global pMOR')
ylabel('Displacement at load point')
xlabel('Frequency (rad/s)')





%% Local pMOR

% Choose and initaliaze parameter space
%Parameter is length, switch to height later

P = [0.214 0.5 0.75 0.9 1.25];
%For some reason anything below 0.214 as a length does not allow the 
% SO-IRKA algorithm to function properly.

tol = 1e-10;
maxiter = 10;
s0 = 1i*linspace(0,10000,5); %r = 5

%Loop over the parameter space
k = length(P);

V = []; %Not sure if this is effective or if preallocating would be better

%These matrices will store the reduced system matrices for all parameter
%values
Mp = [];
Dp = [];
Kp = [];
fp = [];
Cp = [];

for i = 1:k
    
    %Create model matrices for that parameter
    [M, D, K, f, C] = fem_beam(P(i),10);

    %Obtain basis matrix and reduced system matrices for that parameter
    [Kr, Dr, Mr, fr, Cr, sr, Vr] = SO_IRKA_SISO(K, D, M, f, C, s0, 'os', tol, maxiter);

    %Concatenate
    V = [V, Vr];

    Mp = [Mp, Mr];
    Dp = [Dp, Dr];
    Kp = [Kp, Kr];
    fp = [fp, fr];
    Cp = [Cp; Cr];
        
end

%Perform SVD to obtain r most significant basis vectors
[U, S, Vsd] = svd(V);

r = 2*length(s0); %Remember the SO-IRKA algorithm doubles r
%r = 24;
R = U(:,1:r);

%Transform reduced system matrices to generalized coordinate system R

Mgp = [];
Dgp = [];
Kgp = [];
fgp = [];
Cgp = [];

for i = 1:k
    
    Vk_indices = (i-1)*r + 1 : i*r;
    Vk  = V(:,Vk_indices);
    
    T = inv(R'*Vk); % transformation matrix

    M_tilde_k = T'*Mp(:,Vk_indices)*T;
    D_tilde_k = T'*Dp(:,Vk_indices)*T;
    K_tilde_k = T'*Kp(:,Vk_indices)*T;
    f_tilde_k = T'*fp(:,i);
    C_tilde_k = Cp(i,:);

    Mgp = [Mgp, M_tilde_k];
    Dgp = [Dgp, D_tilde_k];
    Kgp = [Kgp, K_tilde_k];
    fgp = [fgp, f_tilde_k];
    Cgp = [Cgp; C_tilde_k];

end


%Use cubic splines: get solution for each reduced and transformed matrix
%for each element in the sampled parameter space, then interpolate between the
%solutions.

%Get solution for frequency range for each parameter value
%Interpolate solution for each frequency range


%Array with one row for each the solution for each parameter value, columns
%are the frequency values
result_interpolate = zeros(k,length(s));
for j = 1:k
    for i=1:length(s)
        gk_indices = (j-1)*r + 1 : j*r; % indices corresponding to generalized k-th basis
        result_interpolate(j,i) = Cgp(j,:) * ((s(i)^2*Mgp(:,gk_indices) + s(i)*Dgp(:,gk_indices) + Kgp(:,gk_indices)) \ fgp(:,j));
        %Careful, s should be the initial frequency range, not the eigenvalues
        %obtained from IRKA
    end
end


%Parameter space for which we want interpolations, here we see the benefit
%of increasing offline costs.
p = [0.1 0.45 0.9 1 1.1 1.5];

result_local = zeros(length(p),length(s));
for j = 1:length(p)
    for i=1:length(s) %s should not be involved in interpolation
        
        result_local(:,i) = spline(P,result_interpolate(:,i),p);
        %result is interpolated for each frequency taking the parameter
        %space as X and the result for that generalized coordinate system
        %result as Y. Putting these together leads to complete interpolated
        %version for entire frequency range for that certain parameter in
        %from the p parameter range (range for which we want to obtain
        %interpolated results). We do this for each parameter in p.

    end
end


%Plot reponse using INTERPOLATION OF LOCAL MATRICES
fig = figure('Name','Frequency response using local pMOR');
set(fig,'defaulttextinterpreter','latex')
semilogy(abs(s),abs(result_local(4,:))) % Row 4 contains result for L = 1 m
xlim([0 10000])
ylim([1e-8 1e-1])
title('Local pMOR')
ylabel('Displacement at load point')
xlabel('Frequency (rad/s)')

%Plot reponse comparisons among all methods
fig = figure('Name','Frequency response comparison inc. local');
set(fig,'defaulttextinterpreter','latex')
semilogy(abs(s),abs(result_full), 'LineWidth', 10)
hold on
semilogy(abs(s),abs(result_irka), 'LineWidth', 7)
hold on
semilogy(abs(s),abs(result_global), 'LineWidth', 3)
hold on
semilogy(abs(s),abs(result_local(4,:)), 'LineWidth', 1)
xlim([0 10000])
ylim([1e-8 1e-1])
legend('Full', 'SO-IRKA', 'Global MOR', 'Local MOR')
title('Frequency response comparison inc. local pMOR')
ylabel('Displacement at load point')
xlabel('Frequency (rad/s)')

