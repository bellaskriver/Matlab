echo on

%----- Topology matrix Edof -------------------------------------
Edof=[1 7 8 5 6;
2 7 8 5 6;
3 1 2 7 8];

%----- Stiffness matrix K and load vector f ---------------------
K=zeros(8);
f=zeros(8,1); f(7)=40e3; f(8)=-30e3;

%----- Element properties ---------------------------------------
E=2.0e11;
% A1=6.0e-4; A2=3.0e-4; A3=10.0e-4;
A=0.25e-3
% ep1=[E A1]; ep2=[E A2]; ep3=[E A3];
ep=[E A];

%----- Element coordinates --------------------------------------
ex1=[3 6]; ex2=[3 3]; ex3=[0 3];
ey1=[0 4]; ey2=[0 4]; ey3=[4 0];

%----- Element stiffness matrices ------------------------------
Ke1=bar2e(ex1,ey1,ep)
Ke2=bar2e(ex2,ey2,ep)
Ke3=bar2e(ex3,ey3,ep)

%----- Assemble Ke into K ---------------------------------------
K=assem(Edof(1,:),K,Ke1);
K=assem(Edof(2,:),K,Ke2);
K=assem(Edof(3,:),K,Ke3)

%----- Solve the system of equations ----------------------------
bc= [1 0;2 0;3 0;4 0;5 0;6 0];
[a,r]=solveq(K,f,bc)

%----- Element forces -------------------------------------------
ed1=extract_ed(Edof(1,:),a);
N1=bar2s(ex1,ey1,ep,ed1)
ed2=extract_ed(Edof(2,:),a);
N2=bar2s(ex2,ey2,ep,ed2)
ed3=extract_ed(Edof(3,:),a);
N3=bar2s(ex3,ey3,ep,ed3)

%----- Draw deformed truss ---------------------------------------
figure(1)
plotpar=[2 1 0];
eldraw2(ex1,ey1,plotpar);
eldraw2(ex2,ey2,plotpar);
eldraw2(ex3,ey3,plotpar);
sfac=scalfact2(ex1,ey1,ed1,0.1);

plotpar=[1 2 1];
eldisp2(ex1,ey1,ed1,plotpar,sfac);
eldisp2(ex2,ey2,ed2,plotpar,sfac);
eldisp2(ex3,ey3,ed3,plotpar,sfac);

axis([-1 7 -0.7 4.5]);
scalgraph2(sfac,[1e-3 0 -0.3]);
title('Displacements')

%----- Draw normal force diagram --------------------------------
figure(2)
plotpar=[2 1];
sfac=scalfact2(ex1,ey1,N2(:,1),0.1);
secforce2(ex1,ey1,N1(:,1),plotpar,sfac);
secforce2(ex2,ey2,N2(:,1),plotpar,sfac);
secforce2(ex3,ey3,N3(:,1),plotpar,sfac);
axis([-1 7 -1 4.5]);
scalgraph2(sfac,[5e4 0 -0.3]);
title('Normal force')

%---------------------------- end -------------------------------

echo off