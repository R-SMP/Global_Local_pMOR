function [Kr, Gr, Mr, Br, Cr, s] = SO_IRKA_SISO(K, G, M, B, C, s, proj, tol, maxiter)
%SO_IRKA_SISO_us SO-IRKA with adaptive order
%   Inputs: M, G, K, B, C : system
%           M : mass matrix
%           G : damping matrix
%           K : stiffness matrix
%           C : not sure what it is, but it comes into play when projection is two sided 
%           s : initial expansion points  %in IRKA they are assumed to be
%           imaginary
%           tol : convergence tolerance %use same error as in Exercise 3.2
%           proj : projection method (one-sided, two-sided) %Not sure about
%           this
%           maxiter : max number of iterations %10 was the number in Ex
%           3.2, but this could be a nice parametric study

if ~(strcmpi(proj, 'os') || strcmpi(proj, 'ts'))
    error('unknown projection method')
end

s = s(:);
r = length(s);
s = cplxpair(s); % sorts the elements along different dimensions of a complex array, grouping together complex conjugate pairs
Vr = zeros(size(M,1),r); %r is the dimension of the number of expansion points we want to use
Wr = zeros(size(M,1),r);

%Alternative
%V = zeros(length(B),length(s));

%Compute initial basis
fprintf('Starting initial basis computation...')
timer_initial = tic;

%Alternative
% for ii=1:length(s0)
%     V(:,ii) = (s0(ii)^2*M + s0(ii)*G + K) \ B;
% end

for ii = 1:r/2
    tmp1 = s(2*ii)^2*M + s(2*ii)*G + K; 
    x = tmp1\B;
    %The two previous steps are the initial reduction basis, could be just
    %one step
    Vr(:,2*ii-1) = real(x); %Not sure why we split this, an array can be of complex numbers
    Vr(:,2*ii) = imag(x);
    
    if strcmp(proj, 'ts')
        x = C / tmp1;
        Wr(:,2*ii-1) = real(x);
        Wr(:,2*ii) = imag(x);
    end
end

[Vr, ~] = qr(Vr,0); %Orthogonalize
if strcmp(proj, 'ts')
    [Wr, ~] = qr(Wr,0);
else
    Wr = Vr;
end
fprintf(' finished in %f s.\n', toc(timer_initial))


%Start iteration
iter = 1;
err = 1;
while (iter <= maxiter) && (err > tol)
    fprintf('Starting iteration i=%d with order r=%d...\n',iter,r)
    s_old = s;
    
    Mr = Wr' * M * Vr;
    Gr = Wr' * G * Vr;
    Kr = Wr' * K * Vr;
    
    lam = polyeig(Kr,Gr,Mr);
    
    % remove real and unstable eigenvalues %This is interesting
    lam(imag(lam) == 0) = [];
    lam(real(lam) > 0) = [];
    
    % sort for the eigenvalues next to the imaginary axis %Not quite sure
    % why we do it this way, 
    %  lam = sort(polyeig(Kr,Cr,Mr)) does the exact same thing, although it
    %  does skip the intermediate step of removing real and unstable
    %  eigenvalues
    [~,order] = sort(abs(imag(lam)),'ascend');
    
    % choose expansion point as mirror images of eigenvalues, new s is
    % assigned/created here
    s = -lam(order);
    s = s(1:r); %we've chosen the first r eigenvalues because the polyeig led to 2r eigenvalues

    % sort again
    s = cplxpair(s);

    %Take the new s points and calculate new basis

    % expand basis
    Vr = zeros(size(M,1),r);
    Wr = zeros(size(M,1),r);
    for ii = 1:r/2
        tmp1 = s(ii*2)^2*M + s(ii*2)*G + K;
        x = tmp1\B;
        Vr(:,2*ii-1) = real(x);
        Vr(:,2*ii) = imag(x);
        
        if strcmp(proj, 'ts')
            x = C / tmp1;
            Wr(:,2*ii-1) = real(x);
            Wr(:,2*ii) = imag(x);
        end
        
    end
    
    % orthogonalize
    [Vr, ~] = qr(Vr,0);
    if strcmp(proj, 'ts')
        [Wr, ~] = qr(Wr,0);
    else
        Wr = Vr;
    end
    
    % convergence
    err = norm(sort(s,'ComparisonMethod','abs') - ...
        sort(s_old,'ComparisonMethod','abs'))/norm(s_old);
    fprintf('\terr = %e.\n',err)
    
    iter = iter + 1;
    
end

Mr = Wr' * M * Vr;
Gr = Wr' * G * Vr;
Kr = Wr' * K * Vr;
Br = Wr' * B;
Cr = C * Vr; %Not necessary unless in two sided projection
% hold off

end

