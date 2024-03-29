clear
close all

%Add SOIRKA to path
addpath('soirka')
%% Inputs
numSamples = 1000; % number of frequency samples
numElements = 50; % number of beam FEM elements
height = 0.01;
length_beam = 1;

%Create beam model
[M, D, K, f, C] = fem_beam(length_beam,height,numElements);

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
fig = figure('Name','Frequency response using full solution');
set(fig,'defaulttextinterpreter','latex')
semilogy(abs(s),abs(result_full))
xlim([0 10000])
ylim([1e-12 1e-0])
title('Frequency response using full solution')
ylabel('Displacement at load point')
xlabel('Frequency (rad/s)')

%% Compute SO-IRKA solution using algorithm implementation from Exericse 3.2

s = 1i*logspace(0,4,numSamples);

%We can change the "r" here, we have chosen r = 5
r = 5;
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
fig = figure('Name','Frequency response using SO-IRKA');
set(fig,'defaulttextinterpreter','latex')
semilogy(abs(s),abs(result_irka))
xlim([0 10000])
ylim([1e-12 1e-0])
title('Frequency response using SO-IRKA')
ylabel('Displacement at load point')
xlabel('Frequency (rad/s)')

%Plot overlaid results and error plot
fig = figure();
set(fig,'defaulttextinterpreter','latex')
subplot(2,1,1)
semilogy(abs(s),abs(result_full),...
    abs(s),abs(result_irka))
xlim([0 10000])
ylim([1e-12 1e-0])
legend('Full system', 'Reduced System (SO-IRKA)')
title('Comparison of full and reduced models')
ylabel('Displacement at load point')
subplot(2,1,2)
semilogy(abs(s),abs(error_irka))
error_str = sprintf('%1d',error_norm_irka);
txt = {'Error norm: ',error_str};
text(500,0.01,txt)
xlim([0 10000])
ylim([1e-12 100])
legend('Reduced system error (SO-IRKA)',...
    'Location','northeast')
ylabel('Relative approximation error')
xlabel('Frequency (rad/s)')

% fprintf("Error norm for IRKA: %.3e\n", ...
%     error_norm_irka);


%% Global pMOR

% Choose and initaliaze parameter space
% Parameter is length, switch to height later

%P = [0.214 0.5 0.75 0.9 1.25]; %original
%P = [0.214 0.3 0.4 0.5 0.6 0.75 0.9 1.25];

P = [0.02 0.03 0.04 0.05 0.06 0.07 0.08 1.2];
%For some reason anything below 0.214 as a length does not allow the 
% SO-IRKA algorithm to function properly.

tol = 1e-10;
maxiter = 10;
s0 = 1i*linspace(0,10000,5);

tic

%Loop over the parameter space
k = length(P);

V = []; %Not sure if this is effective or if preallocating would be better

for i = 1:k
    
    %Create model matrices for that parameter
    [M, D, K, f, C] = fem_beam(length_beam,P(i),numElements);

    %Obtain basis matrix for that parameter
    [~, ~, ~, ~, ~, ~, Vr] = SO_IRKA_SISO(K, D, M, f, C, s0, 'os', tol, maxiter);

    %Concatenate reduced bases
    V = [V, Vr];
        
end

%Perform SVD to remove rank deficient components
[U, S, Vsd] = svd(V);

for i = 1:size(S,1)

    if S(i,i) < 10e-3

        r_index = i; % index beyond which eigenvalues are small enough to be neglected 
        % equivalent to reduced dimension of reduced global basis
        break

    end

end
U = U(:, 1:r_index);
S = S(1:r_index, 1:r_index);
Vsd = Vsd(1:r_index, 1:r_index);

V = U*S*Vsd'; % reduced global bases 

%Sample use case of global basis

%Compute reduced system matrices for that certain height value, 0.1
%chosen here

[M, D, K, f, C] = fem_beam(length_beam,height,numElements);

Mr = V' * M * V;
Dr = V' * D * V;
Kr = V' * K * V;
fr = V' * f;
Cr = C * V;

ctime_global_offline = toc;

s = 1i*logspace(0,4,numSamples); 

tic
result_global = zeros(1,length(s));
for ii=1:length(s)
    lhs = s(ii)^2*Mr + s(ii)*Dr + Kr;
    result_global(ii) = Cr * (lhs \ fr);
    %Careful, s should be the initial frequency range, not the eigenvalues
    %obtained from IRKA
end
ctime_global_online = toc;

%Error calculations
error_global = (result_full - result_global)./result_full;
error_norm_global = norm(result_full - result_global)/norm(result_full);

%Plot reponse using GLOBAL BASIS
fig = figure('Name','Frequency response using Global pMOR');
set(fig,'defaulttextinterpreter','latex')
semilogy(abs(s),abs(result_global))
xlim([0 10000])
ylim([1e-12 1e-0])
title('Frequency response using Global pMOR')
ylabel('Displacement at load point')
xlabel('Frequency (rad/s)')

%Overlaid plot of full, reduced, and GLOBAL BASIS
fig = figure('Name','Frequency response comparison including Global pMOR');
set(fig,'defaulttextinterpreter','latex')
semilogy(abs(s),abs(result_full), 'LineWidth', 6)
hold on
semilogy(abs(s),abs(result_irka), 'LineWidth', 4)
hold on
semilogy(abs(s),abs(result_global), 'LineWidth', 2)
xlim([0 10000])
ylim([1e-12 1e-0])
legend('Full', 'SO-IRKA', 'Global pMOR')
title('Frequency response comparison including Global pMOR')
ylabel('Displacement at load point')
xlabel('Frequency (rad/s)')

%Plot overlaid results and error plot for GLOBAL BASIS
fig = figure();
set(fig,'defaulttextinterpreter','latex')
subplot(2,1,1)
semilogy(abs(s),abs(result_full),...
    abs(s),abs(result_global))
xlim([0 10000])
ylim([1e-12 1e-0])
legend('Full system', 'Reduced System (GLOBAL pMOR)')
title('Comparison of full and reduced models')
ylabel('Displacement at load point')
subplot(2,1,2)
semilogy(abs(s),abs(error_global))
error_str = sprintf('%1d',error_norm_global);
txt = {'Error norm: ',error_str};
text(500,0.01,txt)
xlim([0 10000])
ylim([1e-12 100])
legend('Reduced system error (GLOBAL pMOR)',...
    'Location','northeast')
ylabel('Relative approximation error')
xlabel('Frequency (rad/s)')







%% Local pMOR

% Choose and initaliaze parameter space
%Parameter is length, switch to height later

%P = [0.214 0.5 0.75 0.9 1.25];
%P = [0.1 0.214 0.3 0.4 0.5 0.6 0.75 0.9 1.25];
P = [0.02 0.03 0.04 0.05 0.06 0.07 0.08 1.2];
%For some reason anything below 0.214 as a length does not allow the 
% SO-IRKA algorithm to function properly.

tol = 1e-10;
maxiter = 10;
s0 = 1i*linspace(0,10000,5); %r = 5

tic

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
    [M, D, K, f, C] = fem_beam(length_beam,P(i),numElements);

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
    
    T = inv(R' * Vk); % transformation matrix

    M_tilde_k = T' * Mp(:,Vk_indices)*T;
    D_tilde_k = T' * Dp(:,Vk_indices)*T;
    K_tilde_k = T' * Kp(:,Vk_indices)*T;
    f_tilde_k = T' * fp(:,i);
    C_tilde_k = Cp(i,:) * T;

    Mgp = [Mgp, M_tilde_k];
    Dgp = [Dgp, D_tilde_k];
    Kgp = [Kgp, K_tilde_k];
    fgp = [fgp, f_tilde_k];
    Cgp = [Cgp; C_tilde_k];

end

%---------  INTERPOLATION  ---------%
Mf = zeros(r);
Kf = zeros(r);
Df = zeros(r);
ff = zeros(r,1);
Cf = zeros(1,r);
p = height; % parameter for which interpolation is required, here H = 0.1m
numParam = length(P);

%Use cubic splines: interpolate between the reduced and transformed
%matrices to get the final system matrices corresponding to parameter 'p'
%This needs to be done element-wise

%Get solution for frequency range for that parameter value by solving the
%final system

for i = 1:r

    % My, Dy, Ky, fy, Cy are vectors with matrix elements corresponding to parameters in P 
    fy = fgp(i,:);
    ff(i) = spline(P,fy,p);
    Cy = Cgp(:,i);
    Cf(i) = spline(P',Cy,p);

    My = zeros(1,numParam);
    Ky = zeros(1,numParam);
    Dy = zeros(1,numParam);
    
    % matrices are symmetric so interpolation is required only for the lower
    % triangular part
    for j = 1:i 
       
        for iParam = 1:numParam
            colIndex = (iParam-1)*r + j;
            My(iParam) = Mgp(i,colIndex);
            Ky(iParam) = Kgp(i,colIndex);
            Dy(iParam) = Dgp(i,colIndex);
        end

        Mf(i,j) = spline(P,My,p);
        Kf(i,j) = spline(P,Ky,p);
        Df(i,j) = spline(P,Dy,p);
    
    end
end

% Obtain the final symmetric system matrices
Mf = Mf + transpose(tril(Mf, -1));
Kf = Kf + transpose(tril(Kf, -1));
Df = Df + transpose(tril(Df, -1));

ctime_local_offline = toc;

tic

% Solve the final system model 
result_local = zeros(1, length(s));
for i = 1:length(s)
   result_local(i) = Cf * ((s(i)^2*Mf + s(i)*Df + Kf) \ ff);
end

ctime_local_online = toc;

%Plot reponse using INTERPOLATION OF LOCAL MATRICES
fig = figure('Name','Frequency response using local pMOR');
set(fig,'defaulttextinterpreter','latex')
semilogy(abs(s),abs(result_local))
xlim([0 10000])
ylim([1e-12 1e-0])
title('Frequency response using Local pMOR')
ylabel('Displacement at load point')
xlabel('Frequency (rad/s)')

%Plot reponse comparisons among all methods
fig = figure('Name','Frequency response comparison inc. local');
set(fig,'defaulttextinterpreter','latex')
semilogy(abs(s),abs(result_full),'-')%'LineWidth', 7)
hold on
semilogy(abs(s),abs(result_irka),':')%'LineWidth', 5)
hold on
semilogy(abs(s),abs(result_global),'-.')%'LineWidth', 3)
hold on
semilogy(abs(s),abs(result_local),'--')%'LineWidth', 1)
xlim([0 10000])
ylim([1e-7 1e-1])
legend('Full', 'SO-IRKA', 'Global pMOR', 'Local pMOR')
title('Frequency response comparison inc. Local pMOR')
ylabel('Displacement at load point')
xlabel('Frequency (rad/s)')

%Error calculations
error_local = (result_full - result_local)./result_full;
error_norm_local = norm(result_full - result_local)/norm(result_full);

%Plot overlaid results and error plot for LOCAL BASIS
fig = figure();
set(fig,'defaulttextinterpreter','latex')
subplot(2,1,1)
semilogy(abs(s),abs(result_full),...
    abs(s),abs(result_local))
xlim([0 10000])
ylim([1e-7 1e-1])
legend('Full system', 'Reduced System (LOCAL pMOR)')
title('Comparison of full and reduced models')
ylabel('Displacement at load point')
subplot(2,1,2)
semilogy(abs(s),abs(error_local))
error_str = sprintf('%1d',error_norm_local);
txt = {'Error norm: ',error_str};
text(500,0.00001,txt)
xlim([0 10000])
ylim([1e-12 100])
legend('Reduced system error (LOCAL PMOR)',...
    'Location','southeast')
ylabel('Relative approximation error')
xlabel('Frequency (rad/s)')

 
% xlim([0 10000])
% ylim([1e-12 1e-0])
% %S = sprintf('%.2f',p(1:length(p)));
% S = sprintf('L = %0.2f*', p(1:length(p)));
% C = regexp(S, '*', 'split');
% legend(C{:})
% title('Frequency response comparison for different stiffnesses')
% ylabel('Displacement at load point')
% xlabel('Frequency (rad/s)')