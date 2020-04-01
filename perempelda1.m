function [h,y]=perempelda1(N)


a=0; b=1; %intervallum
alpha=0; beta=1; %peremértékek
h=(b-a)/(N+1);
% A_h
e=ones(N,1);
A_h=(1/h^2)*spdiags([e -2*e-h^2*e e],[-1,0,1],N,N);
%f jobboldal
f=zeros(N,1); f(1)=-alpha/h^2; f(N)=-beta/h^2;
% Numerikus és pontos
y=A_h\f;

