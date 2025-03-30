%Huvudfunktion
function projektuppgift
close all
[L, E, I, q, nfack, npl]=myinit;
mycheck(L,E,I,q,nfack)
f=L./(3.*E.*I);
g=q.*L.^3./(24.*E.*I);
A=diag(0.5*f(2:nfack-1),-1)+diag(f(1:nfack-1)+f(2:nfack),0)+diag(0.5*f(2:nfack-1),1);
b=g(1:nfack-1)+g(2:nfack);
z=A\b;
z=[0;z;0];
kondensationstal=cond(A)
[MA, MB, RA, RB, Mbeam, Mextreme, Vbeam, xbeam] = fortsatthuvudprogram(L, E, I, q, nfack, npl, z);

function [MA, MB, RA, RB, Mbeam, Mextreme, Vbeam, xbeam] = fortsatthuvudprogram(L, E, I, q, nfack, npl, z)
%Steg 5
MA=z(1:nfack);
MB=z(2:nfack+1);

%Steg 6
RA=(-q.*L)/2-(MA-MB)./L;
RB=(-q.*L)/2+(MA-MB)./L;

%Mbeam
Mbeam=[];
Mextreme=[];

for k=1:nfack
x=linspace(0,L(k),npl)';
Mloc=MA(k)+RA(k)*x+(q(k)*x.^2)/2;
Mextreme=[Mextreme myextreme(Mloc,npl)];
Mbeam=[Mbeam;Mloc(1:npl-1)];
end

Mbeam=[Mbeam;Mloc(npl)];

%Vbeam
for k=1:nfack
x=linspace(0,L(k),npl);
V=-RA(k)-q(k).*x;
Vbeam(k)=V(k);
Vbeam=[Vbeam V(2:npl)];
end

%Ändpunkter (xbeam)
i=0;
xbeam=0;
for k=1:nfack
x=linspace(i,i+L(k),npl);
xbeam=[xbeam x(2:npl)];
i=i+L(k);
end

%Kallar på plot och tabell
myplot(xbeam,Mbeam,Vbeam)
mytable(nfack,RA,RB,MA,MB,Mextreme,Mbeam,x)

%Max moment och plats Mbeam
[maxmomentMbeam, Plats_iMaxMbeam] = max(Mbeam);
x_platsMaxMbeam = xbeam(Plats_iMaxMbeam);
[minstaMbeam, Plats_iMinMbeam] = min(Mbeam);
x_platsMinMbeam = xbeam(Plats_iMinMbeam);
maxmomentMbeam
x_platsMaxMbeam
minstaMbeam
x_platsMinMbeam

function extrv=myextreme(Mloc,np1)
extrv=0;
for i=2:np1-1
if Mloc(i)>Mloc(i-1) && Mloc(i)>Mloc(i+1)
extrv=Mloc(i);
end
end

%Värden
function [L, E, I, q, nfack, npl]=myinit
nfack=5;
npl=1000;
L=[1 1 1 1 1]';
E=[1 1 1 1 1]';
I=[1 1 1 1 1]';
q=[-1 -1 0 -1 0]';

%Kontroll
function mycheck(L, E, I, q, nfack)
if length(L)==length(E)
disp('L och E är samma längd')
else
error('L och E är olika långa')
end
if length(I)==length(q)
disp('I och q är samma längd')
else
error('I och q är olika långa')
end
if L(1:nfack)<0
error('Alla längder är inte positiva')
else
disp('Alla längder är positiva')
end

%Plottar
function myplot(xbeam,Mbeam,Vbeam)
figure(1)
plot(xbeam, Mbeam, 'b-')
grid on
title('Bending moment diagram')
xlabel('position')
ylabel('bending moment')
axis ij
figure(2)
plot(xbeam, Vbeam ,'-b')
grid on
title('Shear force diagram')
xlabel('position')
ylabel('shear force')
axis ij

%Tabell
function mytable(nfack,RA,RB,MA,MB,Mextreme,Mbeam,x)
disp(' ')
disp(['Fack ','Vänster stödreaktion ','Höger stödrektion'])
disp(' ')
for i=1:nfack
disp([int2str(i),' ', num2str(RA(i)),' ', num2str(RB(i))])
end
upplag=zeros(1,nfack+1);
upplag(1)=RA(1);
upplag(nfack+1)=RB(nfack);
upplag(2:nfack)=RA(2:nfack)+RB(1:nfack-1);
disp(' ')
disp(['Stöd ','Upplagskrafter'])
disp(' ')
for i=1:nfack+1
disp([int2str(i),' ', num2str(upplag(i))])
end
disp(' ')
disp(['Fack ','Vänster stödmoment ','Höger stödmoment'])
for i=1:nfack
disp([int2str(i),' ',num2str(MA(i)),' ', num2str(MB(i))])
end
Mextreme=Mextreme(Mextreme~=0);
disp(' ')
disp(['Fack ','Maximalt fältmoment'])
for i=1:length(Mextreme)
disp([int2str(i),' ', num2str(Mextreme(i))])
end
disp(' ')