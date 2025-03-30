%Fel i rad 34

echo on

%----- Element stiffness and element load matrices -------------
E1=210e9; E2=210e9;
A1=10e-3; A2=1e-3;
I1=2e-4;

ep1=[E1 A1 I1]; ep2=[E2 A2];
eq2=[0 0]; eq1=[0 -12e3];
ex1=[0 4]; ey1=[0 0];
ex2=[4 8]; ey2=[0 0];
ex3=[8 12]; ey3=[0 0];
ex4=[0 4]; ey4=[0 4];
ex5=[4 4]; ey5=[0 4];
ex6=[4 8]; ey6=[0 4];
ex7=[4 8]; ey7=[4 4];
ex8=[8 8]; ey8=[0 4];
ex9=[8 12]; ey9=[4 0];

[Ke1,fe1]=beam2e(ex1,ey1,ep1,eq1);
[Ke2,fe2]=beam2e(ex2,ey2,ep1,eq1);
[Ke3,fe3]=beam2e(ex3,ey3,ep1,eq1);

Ke4=bar2e(ex4,ey4,ep2);
Ke5=bar2e(ex5,ey5,ep2);
Ke6=bar2e(ex6,ey6,ep2);
Ke7=bar2e(ex7,ey7,ep2);
Ke8=bar2e(ex8,ey8,ep2);
Ke9=bar2e(ex9,ey9,ep2);

%----- Assemble Ke into K ---------------------------------------
[K,f]=assem(Edof1(1,:),K,Ke1,f,fe1);
[K,f]=assem(Edof1(2,:),K,Ke2,f,fe2);
[K,f]=assem(Edof1(3,:),K,Ke3,f,fe3);

K=assem(Edof2(1,:),K,Ke4);
K=assem(Edof2(2,:),K,Ke5);
K=assem(Edof2(3,:),K,Ke6);
K=assem(Edof2(4,:),K,Ke7);
K=assem(Edof2(5,:),K,Ke8);
K=assem(Edof2(6,:),K,Ke9);

%----- Solve the system of equations and compute reactions ------
bc=[1 0; 2 0; 11 0;];
[a,r]=solveq(K,f,bc)

%----- Section forces -------------------------------------------
Ed1=extract_ed(Edof1,a);
Ed2=extract_ed(Edof2,a);

es1=beam2s(ex1,ey1,ep1,Ed1(1,:),eq1,11)
es2=beam2s(ex2,ey2,ep1,Ed1(2,:),eq1,11)
es3=beam2s(ex3,ey3,ep1,Ed1(3,:),eq1,11)

es4=bar2s(ex4,ey4,ep2,Ed2(1,:))
es5=bar2s(ex5,ey5,ep2,Ed2(2,:))
es6=bar2s(ex6,ey6,ep2,Ed2(3,:))
es7=bar2s(ex7,ey7,ep2,Ed2(4,:))
es8=bar2s(ex8,ey8,ep2,Ed2(5,:))
es9=bar2s(ex9,ey9,ep2,Ed2(6,:))

%----- Max moment ----------------------------------------
b1=max(es1(:,3));
b2=max(es2(:,3));
b3=max(es3(:,3));
B=[b1 b2 b3];
B_max=max(abs(B))

%----- Moment diagram ----------------------------------------
figure(1)
plotpar=[2 1];
sfac=scalfact2(ex3,ey3,es3(:,3),0.2);
secforce2(ex1,ey1,es1(:,3),plotpar,sfac);
secforce2(ex2,ey2,es2(:,3),plotpar,sfac);
secforce2(ex3,ey3,es3(:,3),plotpar,sfac);
title('Moment')
figure(1)

%----- Max normal forces ----------------------------------------
n1=max(es1(:,1));
n2=max(es2(:,1));
n3=max(es3(:,1));
n4=max(es4);
n5=max(es5);
n6=max(es6);
n7=max(es7);
n8=max(es8);
n9=max(es9);
N=[n1 n2 n3 n4 n5 n6 n7 n8 n9];
N
max=max(abs(N))

%------------------------ end -----------------------------------

echo off