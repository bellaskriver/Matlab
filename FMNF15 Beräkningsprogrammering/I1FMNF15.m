function problem1
for m=1:2
v=pi/2*m
for n=4:5
varv=[m n]
u=10^n+v
end
end

%% v1
function mysumav(n)
% Kollar om n är ett positivt heltal
if n<1
error('n är för litet')
end
if n-floor(n)~=0
error('n är inte ett heltal')
end
% Beräknar summan och medelvärdet
s1=0;
for k=1:n
s1=s1+k^2;
end
s2=s1/n
% Visar svaret
disp('Summan blir: ')
disp(s1)
disp('Medelvärdet blir: ')
disp(s2)

%% v2

function mysumav(n)
% Beräknar summan och medelvärdet
s1=1:n;
s2=s1.^2
s2=sum(s1.^2)
s3=s2/n
% Visar summan och medelvärdet
disp('Summan blir: ')
disp(s2)
disp('Medelvärdet blir: ')
disp(s3)

%%

function hjarta
xbody=[ 0 1 2 0 -2 -1 0];
ybody=[0 1 0 -3 0 1 0];
plot(xbody,ybody,'k','Linewidth', 1.5)
hold on
fill(xbody,ybody,'r')

%%´

function hjarttapet
% underfunktion - tapet med hjärtan
close all
for i=-2:2
for j=-2:2
tapet(i,j,0.158)
end
end
axis equal
axis([-2.5 2.5 -2.5 2.5])
title('Glad alla hjärtans dag!')
function tapet(xc,yc,alpha)
xbody=[ 0 1 2 0 -2 -1 0];
ybody=[0 1 0 -3 0 1 0];
plot(xc+alpha*xbody,yc+alpha*ybody,'k', 'Linewidth', 1.5)
hold on
fill (xc+alpha*xbody,yc+alpha*ybody,'r')

%% v1

function twofunc (n)
% Skapar en vektor med tal från -pi/2 till pi/2, antalet punkter är n
xvec=linspace(-pi/2,pi/2,n);
% Sätter funktionerna till sig själva
fvec=xvec;
gvec=tan(xvec);
% Plottar f(x)
plot(linspace(0,20,n),linspace(0,20,n),'b --') ;
hold on
% Plottar g(x)
for k=0:6
xvec1=pi*k+xvec ;
plot(xvec1, gvec,'r-')
end
% Graf tillbehör
title('Nice plot')
xlabel('x-values')
ylabel('function values')
legend('Function f(x) = x','Function g(x) = tan(x)')
axis([0 20 -20 20])
grid on;

%% v2

function twofunc(n)
%Skapar en vektor med tal från -pi/2 till pi/2, antalet punkter är n
%Som inte inkluderar de odefinerade punkterna -pi/2 och pi/2
xvec=linspace(-pi/2+0.0000069,pi/2-0.00000069,n);
% Sätter funktionerna till sig själva
fvec=xvec;
gvec=tan(xvec);
hvec=tan(xvec)-xvec;
% Plottar h(x)
for k=0:6
xvec1=pi*k+xvec ;
hvec = tan(xvec1)-xvec1 ;
plot(xvec1, hvec,'r')
hold on
end
% Graft tillbehör
title('Nice plot')
xlabel('x-values')
ylabel('function values')
legend('Function f(x) = x','Function g(x) = tan(x)')
axis([0 20 -20 20])
grid on;