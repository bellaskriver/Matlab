%Fel i rad 227

%----- Beräkna värmeflödet -------------------------------------
format compact

%Topologi
Edof=[1 1 2;
2 2 3;
3 3 4;
4 4 5;
5 5 6;
6 6 7];

%Styvhetsmatris K och vektor f
K=zeros(6);
f=zeros(6,1);

%Element egenskaper
ep1=[1/0.04];
ep2=[0.58/0.12];
ep3=[0.04/0.1];
ep4=[0.15/0.07];
ep5=[1/0.13];

%Elementstyvhetsmatriser
Ke1=spring1e(ep1);
Ke2=spring1e(ep2);
Ke3=spring1e(ep3);
Ke4=spring1e(ep4);
Ke5=spring1e(ep5);

%Assemblering
K=assem(Edof(1,:),K,Ke1);
K=assem(Edof(2,:),K,Ke2);
K=assem(Edof(3,:),K,Ke3);
K=assem(Edof(4,:),K,Ke4);
K=assem(Edof(5,:),K,Ke5);

%Ekvationssystem
bc=[1 -10; 6 22];
[a,r]=solveq(K,f,bc);
Temperatur=a

%----- Beräkna ånghalt -------------------------------------
%Topologi
DEdof=[1 1 2;
2 2 3;
3 3 4];

%Styvhetsmatris K (Dk) och vektor f (Df)
DK=zeros(4);
Df=zeros(4,1);

%Element egenskaper Dep
Dep1=[50e-7/0.12];
Dep2=[175e-7/0.1];
Dep3=[50e-7/0.07];

%Elementstyhetsmatriser DKe
DKe1=spring1e(Dep1);
DKe2=spring1e(Dep2);
DKe3=spring1e(Dep3);

%Assemblering i DK
DK=assem(DEdof(1,:),DK,DKe1);
DK=assem(DEdof(2,:),DK,DKe2);
DK=assem(DEdof(3,:),DK,DKe3);

%Ekvationssystem
Dbc=[1 1.926; 4 6.793];
[Da,Dr]=solveq(DK,Df,Dbc);
Anghalt=Da

%---- Med mättnadsånghalt -------------------------------------------------------
%Räkning börjar från vänster med 1
%För -30<T<0
T1=a(3)
Da(2)

%Ånghalt för skiktet mellan tegel och mineralull
CDB=10.16*((1.486+(T1/100))^12.3)/(273.15+T1)

%Beräknad ånghalt större än mättnadsånghalt
%För 0<T<30
T2=a(4)
Da(3)

%Ånghalt för skiktet mellan mineralull och gasbetong
CDA=625.67*((1.098+(T2(1)/100))^8.02)/(273.15+T2(1))

%Ej större än mättnadsånghalt
%Gör om men sätter ånghalten till mättnadsånghalten
Dbc=[1 1.926; 2 CDB; 4 6.793];
[Da2,Dr]=solveq(DK,Df,Dbc);
kondenspersekundigram=-Dr(2)
kondensperveckaikilogram=-(Dr(2)/1000)*60*60*24*7

%Kod för 9-6a och 9-6c
close all
clear

%----- Topology -------------------------------------------------
Edof=[1 4 5 6 1 2 3;
2 4 5 6 7 8 9;];

%----- Stiffness matrix K and load vector f ---------------------
K=zeros(9); f=zeros(9,1);

%----- Element stiffness and element load matrices -------------
E=210e9;
A1=3e-3; A2=4.8e-3;
I1=9.6e-6; I2=19.2e-6;
ep1=[E A1 I1]; ep2=[E A2 I2];
ex1=[0 0]; ex2=[0 4.8];
ey1=[3 0]; ey2=[3 3];
eq1=[0 0]; eq2=[0 -100e3];
Ke1=beam2e(ex1,ey1,ep1);
[Ke2,fe2]=beam2e(ex2,ey2,ep2,eq2);

%----- Assemble Ke into K ---------------------------------------
K=assem(Edof(1,:),K,Ke1);
[K,f]=assem(Edof(2,:),K,Ke2, f, fe2);

%----- Solve the system of equations and compute reactions ------
%Fast inspänning vid c
bc=[1 0; 2 0; 3 0; 7 0; 8 0; 9 0;];

%Rulllager vid C
%bc=[1 0; 2 0; 3 0; 8 0];
[a,r]=solveq(K,f,bc);

%----- Section forces -------------------------------------------
Ed=extract_ed(Edof,a);
[es1,edi1]=beam2s(ex1,ey1,ep1,Ed(1,:),eq1,21);
[es2,edi2]=beam2s(ex2,ey2,ep1,Ed(2,:),eq2,21);
Inspanningsmoment=r(3)
Forskjutning=a(4)

%----- Draw deformed frame ---------------------------------------
figure(1)
plotpar=[2 1 0];
eldraw2(ex1,ey1,plotpar);
eldraw2(ex2,ey2,plotpar);
sfac=scalfact2(ex2,ey2,edi2,0.1);
plotpar=[1 2 1];
dispbeam2(ex1,ey1,edi1,plotpar,sfac);
dispbeam2(ex2,ey2,edi2,plotpar,sfac);
axis([-1.5 7.5 -0.5 5.5]);
scalgraph2(sfac,[1e-2 0.5 0]);
title('Displacements')

%----- Draw normal force diagram --------------------------------
figure(2)
plotpar=[2 1];
sfac=scalfact2(ex1,ey1,es1(:,1),0.2);
secforce2(ex1,ey1,es1(:,1),plotpar,sfac);
secforce2(ex2,ey2,es2(:,1),plotpar,sfac);
axis([-1.5 7.5 -0.5 5.5]);
scalgraph2(sfac,[3e4 1.5 0]);
title('Normal force')

%----- Draw shear force diagram ---------------------------------
figure(3)
plotpar=[2 1];
sfac=scalfact2(ex2,ey2,es2(:,2),0.2);
secforce2(ex1,ey1,es1(:,2),plotpar,sfac);
secforce2(ex2,ey2,es2(:,2),plotpar,sfac);
axis([-1.5 7.5 -0.5 5.5]);
scalgraph2(sfac,[3e4 0.5 0]);
title('Shear force')

%----- Draw moment diagram --------------------------------------
figure(4)
plotpar=[2 1];
sfac=scalfact2(ex2,ey2,es2(:,3),0.2);
secforce2(ex1,ey1,es1(:,3),plotpar,sfac);
secforce2(ex2,ey2,es2(:,3),plotpar,sfac);
axis([-1.5 7.5 -0.5 5.5]);
scalgraph2(sfac,[3e4 0.5 0]);
title('Moment')

% uppg. 9-6b, 9-6d, 9-7 & 9-8

close all
echo on

% ----- Topology ----------------------------------------------------------
Edof=[1 4 5 6 1 2 3;
2 4 5 6 7 8 9];

% ----- Element properties and global coordinates -------------------------
E=210e9;
A1=3e-3; A2=4.8e-3;
I1=9.6e-6; I2=19.2e-6;
ep1=[E A1 I1]; ep2=[E A2 I2];
eq1=[0]; eq2=[-100e3];
ex1=[0 0]; ex2=[0 4.8];
ey1=[3 0]; ey2=[3 3];
%uppg. 9-8 a & b
%eq2=[-500e3]

% ----- Initial axial forces ----------------------------------------------
QX1=0.0001; QX2=0;
QX01=1;

% ----- Iteration for convergence -----------------------------------------
format compact
eps=1e-6;
n=0;
while(abs((QX1-QX01)/QX01)>eps)
n=n+1;
K=zeros(9,9);
f=zeros(9,1);
[Ke1]=beam2ge(ex1,ey1,ep1,QX1);
[Ke2,fe2]=beam2ge(ex2,ey2,ep2,QX2,eq2);
K=assem(Edof(1,:),K,Ke1);
[K,f]=assem(Edof(2,:),K,Ke2,f,fe2);


% Fast inspänd vid C
% bc=[1 0; 2 0; 3 0; 7 0; 8 0; 9 0];
% Rullager vid C
bc=[1 0; 2 0; 3 0; 8 0];
[a,r]=solveq(K,f,bc);
Ed=extract_ed(Edof,a);
QX01=QX1;
[es1,QX1,edi1]=beam2gs(ex1,ey1,ep1,Ed(1,:),QX1,eq1,21);
[es2,QX2,edi2]=beam2gs(ex2,ey2,ep2,Ed(2,:),QX2,eq2,21);


if(n==1);
edi10=edi1;
edi20=edi2;
% K0=K läggs till till uppg. 9-7.b för att beräkna
% knäckningssäkerheten
K0=K;
end;
if(n>20)
break
end
end
disp('The solution does not converge')
forskjutning=a(4)
Inspanningsmoment=r(3)
axiell_kraft=QX1
% Tillhörande uppg. 9-7.b
b=bc(:,1);
[lambda,phi]=eigen(K,K0,b);
nmods=size(lambda);
one=ones(nmods);
alpha=one./(one-lambda);
knackningssakerhet=alpha(1)
phi(:,1);

%----- Draw deformed frame ---------------------------------------
figure(1)
plotpar=[3 1 0];
eldraw2(ex1,ey1,plotpar);
eldraw2(ex2,ey2,plotpar);
sfac=scalfact2(ex2,ey2,edi2,0.1);
plotpar=[1 2 0];
dispbeam2(ex1,ey1,edi1,plotpar,sfac);
dispbeam2(ex2,ey2,edi2,plotpar,sfac);
plotpar=[2 4 0];
dispbeam2(ex1,ey1,edi10,plotpar,sfac);
dispbeam2(ex2,ey2,edi20,plotpar,sfac);
axis([-1.5 7.5 -0.5 5.5]);
scalgraph2(sfac,[0.1 0.5 0]);
title('Displacements')

%----- Draw normal force diagram --------------------------------
figure(2)
plotpar=[2 1];
sfac=scalfact2(ex1,ey1,es1(:,1),0.2);
secforce2(ex1,ey1,es1(:,1),plotpar,sfac);
secforce2(ex2,ey2,es2(:,1),plotpar,sfac);
axis([-1.5 7.5 -0.5 5.5]);
scalgraph2(sfac,[2e5 1.5 0]);
title('Normal force')

%----- Draw shear force diagram ---------------------------------
figure(3)
plotpar=[2 1];
sfac=scalfact2(ex2,ey2,es2(:,2),0.2);
secforce2(ex1,ey1,es1(:,2),plotpar,sfac);
secforce2(ex2,ey2,es2(:,2),plotpar,sfac);
axis([-1.5 7.5 -0.5 5.5]);
scalgraph2(sfac,[2e5 0.5 0]);
title('Shear force')

%----- Draw moment diagram --------------------------------------
figure(4)
plotpar=[2 1];
sfac=scalfact2(ex2,ey2,es2(:,3),0.2);
secforce2(ex1,ey1,es1(:,3),plotpar,sfac);
secforce2(ex2,ey2,es2(:,3),plotpar,sfac);
axis([-1.5 7.5 -0.5 5.5]);
scalgraph2(sfac,[2e5 0.5 0]);
title('Moment')

%------------------------ end ---------------------------------------------
echo off

% uppg. 9-8c)
%----------------------------------------------------------------
% PURPOSE
% Buckling analysis of a plane frame.
%--------------------------------------------------------------------------
echo off
% ----- Topology ----------------------------------------------------------
Edof=[1 4 5 6 1 2 3;
2 7 8 9 4 5 6;
3 7 8 9 10 11 12];

% ----- Element properties and global coordinates -------------------------
E=210e9;
A1=3e-3; A2=3e-3; A3=4.8e-3;
I1=9.6e-6; I2=9.6e-6; I3=19.2e-6;
ep1=[E A1 I1]; ep2=[E A2 I2]; ep3=[E A3 I3];
eq1=[0]; eq2=[0]; eq3=[-500e3];
ex1=[0 0]; ex2=[0 0]; ex3=[0 4.8];
ey1=[1.5 0]; ey2=[3 1.5]; ey3=[3 3];

% ----- Initial axial forces ----------------------------------------------
QX1=0.0001; QX2=0; QX3=0;
QX01=1;

% ----- Iteration for convergence -----------------------------------------
eps=1e-6;
n=0;
while(abs((QX1-QX01)/QX01)>eps)
n=n+1;
K=zeros(12,12);
f=zeros(12,1);
[Ke1]=beam2ge(ex1,ey1,ep1,QX1);
[Ke2]=beam2ge(ex2,ey2,ep2,QX2);
[Ke3,fe3]=beam2ge(ex3,ey3,ep3,QX3,eq3);
K=assem(Edof(1,:),K,Ke1);
K=assem(Edof(2,:),K,Ke2);
[K,f]=assem(Edof(3,:),K,Ke3,f,fe3);
if n==1;
K0=K;
end;
bc=[1 0;2 0;3 0;10 0;11 0;12 0];
[a,r]=solveq(K,f,bc);
Ed=extract_ed(Edof,a);
QX01=QX1;
[es1,QX1,edi1]=beam2gs(ex1,ey1,ep1,Ed(1,:),QX1,eq1,21);
[es2,QX2,edi2]=beam2gs(ex2,ey2,ep2,Ed(2,:),QX2,eq2,21);
[es3,QX3,edi3]=beam2gs(ex3,ey3,ep3,Ed(3,:),QX3,eq3,21);
if(n>20)
break
end
end
disp('The solution does not converge')

% ----- Buckling analysis -------------------------------------------------
echo on
b=bc(:,1);
[lambda,phi]=eigen(K,K0,b);
nmods=size(lambda);
one=ones(nmods);
alpha=one./(one-lambda);
alpha(1)
phi(:,1);