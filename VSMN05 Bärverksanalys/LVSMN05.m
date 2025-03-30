%----------------------- Topology ---------------------------------
Edof1=[1 1 2 3 4 5 6;
2 4 5 6 7 8 9;
3 7 8 9 10 11 12;
4 10 11 12 13 14 15;
5 13 14 15 16 17 18;
20 19 20 29 21 22 30;
21 21 22 30 23 24 31;
22 23 24 31 25 26 32;
23 25 26 32 27 28 33];

Edof2=[6 19 20 1 2;
7 19 20 4 5;
8 19 20 7 8;
9 21 22 4 5;
10 21 22 7 8;
11 21 22 10 11;
12 23 24 7 8;
13 23 24 10 11;
14 23 24 13 14;
15 25 26 10 11;
16 25 26 13 14;
17 25 26 16 17;
18 27 28 13 14;
19 27 28 16 17];

%----- Stiffness matrix K and load vector f ---------------------
s=[-125 -25];
utbojning=[];
for i=1:2
K=zeros(33); f=zeros(33,1);
f(28)=s(i);

%----- Element stiffness and element load matrices -------------
E1=12e9; E2=1e9;
A1=(10*15)*(10^-6); A2=(5*15)*(10^-6);
I1=0.010*(0.015^3)/12;
ep1=[E1 A1 I1]; ep4=[E2 A2];
eq1=[0 0]; eq2=[0 0]; eq3=[0 0];

ex1=[0 0.1]; ey1=[0 0];
ex2=[0.1 0.3]; ey2=[0 0];
ex3=[0.3 0.5]; ey3=[0 0];
ex4=[0.5 0.7]; ey4=[0 0];
ex5=[0.7 0.9]; ey5=[0 0];
ex6=[0.1 0]; ey6=[0.15 0];
ex7=[0.1 0.1]; ey7=[0.15 0];
ex8=[0.1 0.3]; ey8=[0.15 0];
ex9=[0.3 0.1]; ey9=[0.15 4];
ex10=[0.3 0.3]; ey10=[0.15 0];
ex11=[0.3 0.5]; ey11=[0.15 0];
ex12=[0.5 0.3]; ey12=[0.15 0];
ex13=[0.5 0.5]; ey13=[0.15 0];
ex14=[0.5 0.7]; ey14=[0.15 0];
ex15=[0.7 0.5]; ey15=[0.15 0];
ex16=[0.7 0.7]; ey16=[0.15 0];
ex17=[0.7 0.9]; ey17=[0.15 0];
ex18=[0.9 0.7]; ey18=[0.15 0];
ex19=[0.9 0.9]; ey19=[0.15 0];
ex20=[0.1 0.3]; ey20=[0 0];
ex21=[0.3 0.5]; ey21=[0 0];
ex22=[0.5 0.7]; ey22=[0 0];
ex23=[0.7 0.9]; ey23=[0 0];

Ke1=beam2e(ex1,ey1,ep1);
Ke2=beam2e(ex2,ey2,ep1);
Ke3=beam2e(ex3,ey3,ep1);
Ke4=beam2e(ex4,ey4,ep1);
Ke5=beam2e(ex5,ey5,ep1);
Ke20=beam2e(ex20,ey20,ep1);
Ke21=beam2e(ex21,ey21,ep1);
Ke22=beam2e(ex22,ey22,ep1);
[Ke23,fe23]=beam2e(ex23,ey23,ep1,eq2);

Ke6=bar2e(ex6,ey6,ep4);
Ke7=bar2e(ex7,ey7,ep4);
Ke8=bar2e(ex8,ey8,ep4);
Ke9=bar2e(ex9,ey9,ep4);
Ke10=bar2e(ex10,ey10,ep4);
Ke11=bar2e(ex11,ey11,ep4);
Ke12=bar2e(ex12,ey12,ep4);
Ke13=bar2e(ex13,ey13,ep4);
Ke14=bar2e(ex14,ey14,ep4);
Ke15=bar2e(ex15,ey15,ep4);
Ke16=bar2e(ex16,ey16,ep4);
Ke17=bar2e(ex17,ey17,ep4);
Ke18=bar2e(ex18,ey18,ep4);
Ke19=bar2e(ex19,ey19,ep4);

%----- Assemble Ke into K ---------------------------------------
K=assem(Edof1(1,:),K,Ke1);
K=assem(Edof1(2,:),K,Ke2);
K=assem(Edof1(3,:),K,Ke3);
K=assem(Edof1(4,:),K,Ke4);
K=assem(Edof1(5,:),K,Ke5);
K=assem(Edof1(6,:),K,Ke20);
K=assem(Edof1(7,:),K,Ke21);
K=assem(Edof1(8,:),K,Ke22);
[K,f]=assem(Edof1(9,:),K,Ke23,f,fe23);

K=assem(Edof2(1,:),K,Ke6);
K=assem(Edof2(2,:),K,Ke7);
K=assem(Edof2(3,:),K,Ke8);
K=assem(Edof2(4,:),K,Ke9);
K=assem(Edof2(5,:),K,Ke10);
K=assem(Edof2(6,:),K,Ke11);
K=assem(Edof2(7,:),K,Ke12);
K=assem(Edof2(8,:),K,Ke13);
K=assem(Edof2(9,:),K,Ke14);
K=assem(Edof2(10,:),K,Ke15);
K=assem(Edof2(11,:),K,Ke15);
K=assem(Edof2(12,:),K,Ke17);
K=assem(Edof2(13,:),K,Ke18);
K=assem(Edof2(14,:),K,Ke19);

%----- Solve the system of equations and compute reactions ------
bc=[2 0; 16 0; 27 0];
[a,r]=solveq(K,f,bc);
utbojning(i)=a(17);
end
styvhet=(s(2)-s(1))/(utbojning(2)-utbojning(1));
disp(['Styvheten i bron beräknas till ', num2str(styvhet),' kg/s2'])