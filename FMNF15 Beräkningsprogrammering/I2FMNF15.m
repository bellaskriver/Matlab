%% Problem 1
function egenunder
close all

%Huvudfunktion
r=1;
for n=1:10;
[x]=ufunc(n,r);
resultvec(n)=(x);
end

%Underfunktion
function [x]=ufunc(n,r)
x=n*r*cos(n);

%% Probelm 2
function nonlin
%Ger värden för alla metoder
dist=0.01;
for k=1:5
xk=pi/2+(k)*pi+dist; %gissningar på rötter
xk2=pi/2+(k+1)*pi-dist;
nollor(k) = mybisect(xk,xk2)
[iter(k),newtonnollor(k)]=g(xk2); %metodanrop för att beräkna nollställe med gissning b
fzero1(k) = fzero(@f, xk2);
end

%Kurvan
n=5;
x=linspace(0,2*pi,n);
y=tan(x)-x;

%Tabeller
disp(' ')
diary nonlin.txt
disp(['Point no.','Intervallhalvering ','Fzero'])
disp(' ')
for i=1:n
entry1=sprintf('%2u',i);
entry2=sprintf('%8.5f' ,nollor(i));
entry3=sprintf('%8.5f',newtonnollor(i));

entry4=sprintf('%2u',iter(i));
entry5=sprintf('%8.5f',fzero1(i));
disp([' ',entry1,' ',entry2, entry4,' ', entry5])
end
diary off

%Intervallhalvering
function myzero=mybisect(a,b)
for k=1:53
c=(a+b)/2;
if sign(f(c)) == sign(f(a))
a=c;
else
b=c;
end
end
bisecnollstalle=(a+b)/2;
myzero=bisecnollstalle;
if abs(f(bisecnollstalle)) > 1e-10
disp('Varning: intervallhalveringsmetoden har kanske inte konvergerat')
end

%Newtons metod
function [iter,newtnollstalle]=g(xk)
tol=10*eps; % eps är maskinepsilon i MATLAB
dxk=xk; % tar oss förbi "while", första varvet
iter=0; % varvräknare
while abs(dxk/xk)>tol && iter<20
dxk=-f(xk)/fprim(xk);
xk=xk+dxk;
iter=iter+1;
end
iter;
newtnollstalle=xk;

%Funktionen f(x)=tan(x)-x
function fout=f(x)
fout=tan(x)-x;

%Derivatan av funktionen f(x)
function fpout=fprim(x)
fpout=1/(cos(x))^2-1;

%% Problem 3
 
function matvec
format compact
n=input('Ange n, n skall vara större än 2500: ');
cen=rand(n,1);
sub=rand(n-1,1);
sup=rand(n-1,1);
A=diag(sub,-1)+diag(cen, 0)+diag(sup,1);
x=rand(n,1);
disp('Setup done')

% *** Algoritm 1 ***
tic
b1=A*x;
t1=toc

% *** Algoritm 2 ***
tic
b2=zeros(n,1);
for i=1:n
for j=1:n
b2(i)=b2(i)+A(i,j)*x(j);
end
end
t2=toc
check2=norm(b2-b1)

% *** Algoritm 3 ***
tic
b3=zeros(n,1);
b3(n)=A(n,n-1)*x(n-1)+A(n,n)*x(n);
b3(1)=A(1,1)*x(1)+A(1,2)*x(2);
for i=2:n-1
for j=i-1:i+1
b3(i)=b3(i)+A(i,j)*x(j);
end
end
t3=toc
check3=norm(b3-b1)