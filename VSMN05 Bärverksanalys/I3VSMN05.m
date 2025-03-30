echo on

%----- Topology -------------------------------------------------
Edof=[1 4 5 6 1 2 3;
2 4 5 6 7 8 9]

%----- Stiffness matrix K and load vector f ---------------------
K=zeros(9); f=zeros(9,1); f(6)=15e+3;

%----- Element stiffness and element load matrices -------------
E=210e9;
A1=3e-3; A2=4.8e-3;
I1=9.6e-6; I2=19.2e-6;
ep1=[E A1 I1]; ep2=[E A2 I2];
ex1=[0 0]; ex2=[0 4.8];
ey1=[3 0]; ey2=[3 3];
eq1=[0 10e+3]; eq2=[0 -20e+3];
[Ke1,fe1]=beam2e(ex1,ey1,ep1,eq1)
[Ke2,fe2]=beam2e(ex2,ey2,ep2,eq2)

%----- Assemble Ke into K ---------------------------------------
[K,f]=assem(Edof(1,:),K,Ke1,f,fe1);
[K,f]=assem(Edof(2,:),K,Ke2,f,fe2)

%----- Solve the system of equations and compute reactions ------
bc=[1 0; 2 0; 3 0; 7 0; 8 0; 9 0]
[a,r]=solveq(K,f,bc)

%----- Section forces -------------------------------------------
Ed=extract_ed(Edof,a);
[es1,edi1]=beam2s(ex1,ey1,ep1,Ed(1,:),eq1,21);
[es2,edi2]=beam2s(ex2,ey2,ep1,Ed(2,:),eq2,21);

%----- Draw deformed frame ---------------------------------------
figure(1)
plotpar=[2 1 0];
eldraw2(ex1,ey1,plotpar);
eldraw2(ex2,ey2,plotpar);
sfac=scalfact2(ex1,ey1,edi1,0.1);
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

%------------------------ end -----------------------------------

echo off