function [h,y]=peremneumann1(N)

a=0; b=1;
alpha=1; beta=exp(1);

%%Racs
h = (b-a)/(N+1);

%%A_h_hullam
e = ones(N+2,1);
A_h_hullam=(1/h^2)*spdiags([e -2*e e],[-1,0,1],N+2,N+2);
% Neumann perem (elso)
A_h_hullam(1,1)=-1/h;
A_h_hullam(1,2)=1/h;
% Dirichlet perem
A_h_hullam(N+2,N+1)=0;
A_h_hullam(N+2,N+2)=1;

%% f_hullam
x = linspace(a,b,N+2)';
f = exp(x); f(1)=alpha; f(N+2)=beta;

%% LAER mo-sa
y=A_h_hullam\f;








