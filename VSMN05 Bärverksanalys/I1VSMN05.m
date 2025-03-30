echo on
format compact

%----- Topology matrix Edof -------------------------------------
Edof=[1 1 2;
2 2 3;
3 3 4];

%----- Stiffness matrix K and load vector f ---------------------
K=zeros(4,4)
f=zeros(4,1); f(3)=1

%----- Element stiffness matrices ------------------------------
k=1; ep1=k; ep2=2*k; ep3=k;
Ke1=spring1e(ep1)
Ke2=spring1e(ep2)
Ke3=spring1e(ep3)

% ----- Assemble Ke into K ---------------------------------------
K=assem(Edof(1,:),K,Ke1)
K=assem(Edof(2,:),K,Ke2)
K=assem(Edof(3,:),K,Ke3)

% ----- Solve the system of equations ----------------------------
bc= [1 0; 4 0];
[a,r]=solveq(K,f,bc)

%----- Element forces -------------------------------------------
ed1=extract_ed(Edof(1,:),a)
ed2=extract_ed(Edof(2,:),a)
ed3=extract_ed(Edof(3,:),a)

es1=spring1s(ep1,ed1)
es2=spring1s(ep2,ed2)
es3=spring1s(ep3,ed3)
%---------------------------- end -------------------------------

echo off