close all
format compact
format short

%Assignment 1
%Del a - förskjutningar
%Topologi
Edof=[1 1 2 3 4;
2 9 10 11 12;
3 1 2 11 12;
4 3 4 11 12;
5 3 4 5 6;
6 11 12 13 14;
7 3 4 13 14;
8 5 6 13 14;
9 5 6 7 8;
10 13 14 15 16;
11 5 6 15 16;
12 7 8 15 16;
13 15 16 17 18;
14 7 8 17 18];

%K matriser
K=zeros(18);
f=zeros(18,1); f(18)=-20e3;

%Element egenskaper
A=2124*10^-6;
E=2.1*10^11;
ep=[E A];

%Koordinater
ex=[0 1.5;
0 1.5;
0 1.5;
1.5 1.5;
1.5 3;
1.5 3;
1.5 3;
3 3;
3 4.5;
3 4.5;
3 4.5;
4.5 4.5;
4.5 6;
4.5 6;];
ey=[1.5 1.5;
0 0;
1.5 0;
1.5 0;
1.5 1.5;
0 0;
1.5 0;
1.5 0;
1.5 1.5;
0 0;
1.5 0;
1.5 0;
0 0;
1.5 0;];

%K matris assemblering
for i=1:14
Ke=bar2e(ex(i,:),ey(i,:),ep);
K=assem(Edof(i,:),K,Ke);
end;

%Upplagsvilkor och förskjutningar
bc= [1 0;2 0;9 0; 10 0];
[a,r]=solveq(K,f,bc);
ForskjutningX=a(17)
ForskjutningY=a(18)

%Krafter
ed=extract_ed(Edof,a);
for i=1:14
es=bar2s(ex(i,:),ey(i,:),ep,ed(i,:));
N(i,:)=es(1);
end
N;

%Del b - maximal spänning
%Spänningen
sigma=N./A;
[minS placemin]=min(sigma)
[maxS placemax]=max(sigma)

%Del c - upplagskrafter vid nedre stödet
SupportforceX=r(9)
SupportforceY=r(10)

%Del d – rullager vid nedre upplaget
format compact

%Topologi
Edof=[1 1 2 3 4;
2 9 10 11 12;
3 1 2 11 12;
4 3 4 11 12;
5 3 4 5 6;
6 11 12 13 14;
7 3 4 13 14;
8 5 6 13 14;
9 5 6 7 8;
10 13 14 15 16;
11 5 6 15 16;
12 7 8 15 16;
13 15 16 17 18;
14 7 8 17 18;
15 1 2 9 10];

%K matriser
K=zeros(18);
f=zeros(18,1); f(18)=-20e3;

%Element egenskaper
A=2124*10^-6;
E=2.1*10^11;
ep=[E A];

%Koordinater
ex=[0 1.5;
0 1.5;
0 1.5;
1.5 1.5;
1.5 3;
1.5 3;
1.5 3;
3 3;
3 4.5;
3 4.5;
3 4.5;
4.5 4.5;
4.5 6;
4.5 6;
0 0];
ey=[1.5 1.5;
0 0;
1.5 0;
1.5 0;
1.5 1.5;
0 0;
1.5 0;
1.5 0;
1.5 1.5;
0 0;
1.5 0;
1.5 0;
0 0;
1.5 0;
1.5 0];

%K matris assemblering
for i=1:15
Ke=bar2e(ex(i,:),ey(i,:),ep);
K=assem(Edof(i,:),K,Ke);
end;

%Upplagsvilkor och förskjutningar
bc= [1 0;2 0;9 0];
[a,r]=solveq(K,f,bc);
Forskjutningar=a

%Krafter
ed=extract_ed(Edof,a);
for i=1:15
es=bar2s(ex(i,:),ey(i,:),ep,ed(i,:));
N(i,:)=es(1);
end
N;

%Plot
figure(1)
eldraw2(ex, ey, [1 2 1]);
hold on
eldisp2(ex, ey, ed, [1 4 1]);
