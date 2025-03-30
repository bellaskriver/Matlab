function test
close all
format compact

%Skapar vektorn nvec
nvec=ceil(logspace(1,2,20));

%Beräknar konditionstalet för systemet beroende på dess storlek
for i=1:20;
konditionstal(i)=cond(matris(nvec(i))); %konditionstalet ökar med systemets
storlek
end
konditionstal

%Ritar figur med axlar i formatet 1 2 3
figure(1)
plot(log10(nvec),log10(konditionstal),'*');
xlabel('log10(x)')
ylabel('log10(y)')
grid on

%Ritar figur med axlar i formatet 10^1 10^2 10^3
figure(2)
loglog(nvec,konditionstal,'*');
xlabel('x')
ylabel('y')
grid on

%Funktionen på formen y=a*x^p
p=(log(konditionstal(20))-log(konditionstal(1)))/(log(nvec(20))-log(nvec(1)))
a=konditionstal(20)/nvec(20)^p

%Beräkna hur stora ekvationssystem är värda att lösa
%med gränser på 0,1 (egenvald) och relativ nogrannhet 0,001
x=(0.1/a/0.0001)^(1/p)

%matrisen D
function D=matris(n)
sup=ones(n-1,1);
cen=-2*ones(n,1);
sub=ones(n-1,1);
D=diag(sub,-1)+diag(cen,0)+diag(sup,1);

function problem3
close all
lambda=[-10 -100 -1000];
n=99;
u=zeros(n+2,3);
%Forloop för att skapa matrisen D och kolonnvektorn b
for i=1:3
D=matris(n,lambda(i));
b=[1; zeros(n-2,1); 1];
u(2:end-1,i)=D\b;
u(:,i)=[1; u(2:end-1,i);1];
end
u;

%Rita upp diagrammet
xvec=linspace(0,1,n+2);
plot(xvec, u(:,1),'+b')
hold on
plot(xvec, u(:,2),'*g')
hold on
plot(xvec, u(:,3),'or')
hold on
legend('λ=-10','λ=-100','λ=-1000')

%Ekvationssystemet
function D=matris(n,lambda)
h=1/(n+1);
cen=(2-lambda*h^2)*ones(n,1);
sup=-1*ones(n-1,1);
sub=sup;
D=diag(cen)+diag(sup,-1)+diag(sub,1);