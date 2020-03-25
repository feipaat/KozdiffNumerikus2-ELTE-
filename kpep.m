function [h,x,y]=kpep(inta,intb,alpha,beta,N)
%% Kétpontos peremérték-feladat megoldása
%
%     a(x)u''(x)+b(x)u'(x)+c(x)u(x)=f(x)
%     u(a)=\alpha u(b)=\beta
%

%% Bemenõ paramatérek listája:

%     inta intervallum kezdete
%     intb intervallum vége
%     N    intervallumok száma


%% Elõkészületek

% Lépésköz

h=(intb-inta)/(N+1);

%% Az a,b,c,f függvények összerakása
x=linspace(inta,intb,N+2)';
e=ones(N+2,1);
nulla=zeros(N+2,1);
a=feval(@(x)(x.^2+e),x);
b=feval(@(x)(-2*x),x);
c=feval(@(x)(2*e),x);
f=feval(@(x)((x.^2+e).^2),x);
%% A_h összerakása és a perem
A_h=(1/h^2)*spdiags([a-0.5*h*b -2*a+h^2*c a+0.5*h*b],[-1,0,1],N+2,N+2);
A_h(1,1)=1; A_h(1,2)=0;
A_h(N+2,N+1)=0; A_h(N+2,N+2)=1;
f(1)=alpha; f(N+2)=beta;
y=A_h\f;

