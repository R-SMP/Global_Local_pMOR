function [Mr, Dr, Kr, fr, Cr, s0] = irka(M, D, K, f, C, s0, tol, maxiter)
%irka Simple implementation of IRKA algorithm as seen in Exercise 3.2
%   Inputs: M, D, K, f, C, s0, tol, maxiter : system
%           M : mass matrix
%           D : damping matrix
%           K : stiffness matrix
%           C : solution matrix, picks out tip displacement DOF 
%           s0 : initial expansion points  %in IRKA they are assumed to be
%           imaginary
%           tol : convergence tolerance 
%           maxiter : max number of iterations

%fprintf("r and s0_old");
r = length(s0);
s0_old = s0;

%compute the initial reduction basis
%fprintf("compute the initial reduction basis");
V = zeros(length(f),length(s0));
for ii=1:length(s0)
    V(:,ii) = (s0(ii)^2*M + s0(ii)*C + K) \ f;
end
[V, ~] = qr(V,0);


% start IRKA loop
%fprintf("start IRKA loop")
iter = 1;
err = 1;
while (iter < maxiter) && (err > tol)
    fprintf('Starting iteration i=%d with order r=%d...\n',iter,r)
    Mr = V' * M * V;
    Dr = V' * D * V;
    Kr = V' * K * V;
    
    lam = sort(polyeig(Kr,Dr,Mr));
    s0 = -lam(1:r);
    
    err = norm(sort(abs(s0)) - sort(abs(s0_old)))/norm(s0_old);
    iter = iter + 1;
    s0_old = s0;
    
    for ii=1:length(s0)
        V(:,ii) = (s0(ii)^2*M + s0(ii)*C + K) \ f;
    end
    [V, ~] = qr(V,0); 
end

% project the system
fprintf("project the system")
Mr = V' * M * V;
Dr = V' * D * V;
Kr = V' * K * V;
fr = V' * f;
o = sparse(1,length(f));
o(f~=0) = 1; 
Cr = o * V;

end