%% Static condensation example

K = [2 -1 0;
    -1 2 -1;
    0 -1 1];

f = [1; 0; 0];

x = K\f


%% Dynamic condensation example

K = [2 -1 0;
    -1 2 -1;
    0 -1 1];

M = [1 0 0;
    0 1 0;
    0 0 1];

f = [1; 0; 0];

xa = [1 2];

n = length(K);
na = length (xa);
ni = n - na;

% Kaa = zeros(length(xa));
% Kii = zeros(length(K)-length(xa));
% 
% Kai = zeros(length(xa),length(K)-length(xa));
% Kia = zeros(length(K)-length(xa),length(xa));

Kaa = K(1:na, 1:na);
Kii = K(na+1:n, na+1:n);
Kai = K(1:na, na+1:n);
Kia = K(na+1:n, 1:na);

V = [ones(na); Kii\Kia];

Vh = V'; %conjugate transpose aka adjoint


syms w

Mr = Vh*M*V;

Kr = Vh*K*V;

fr = Vh*f;


eqn = det((-w^2)*Mr + w*Kr) == 0;

S = solve(eqn,w)

double(S)

det(double(S(3)^2)*Mr + double(S(3))*Kr)
det(double(S(4)^2)*Mr + double(S(4))*Kr)

det(double(S(3)^2)*Mr + double(S(3))*Kr)
det(double(S(4)^2)*Mr + double(S(4))*Kr)

u1 = linsolve(double(S(3)^2)*Mr + double(S(3))*Kr, zeros(na,1))
u2 = linsolve(double(S(4)^2)*Mr + double(S(4))*Kr, zeros(na,1))













% s_vals = linspace(0,2,100);
% 
% response = zeros(2,100);
% 
% for i = 1:100
% 
%     response(1:2, i) = ((s_vals(i)^2)*Mr + Kr)\fr;
% 
% end
% 
% atan(0)











