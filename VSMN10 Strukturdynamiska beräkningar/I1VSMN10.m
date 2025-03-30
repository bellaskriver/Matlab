%% Project 1 VSMN10 - Hugo Persson & Isabella McAuley Skriver (Group 6)
clear all; close all; format compact; format short;

%% Experiment
R1=readmatrix('resonance_1.txt'); % 2501x8 [Tid, Sensor nederst --> Överst]
R2=readmatrix('resonance_2.txt'); % 19360x8
R3=readmatrix('resonance_3.txt'); % 19360x8
R4=readmatrix('resonance_4.txt'); % 19360x8
R5=readmatrix('resonance_5.txt'); % 19360x8
R6=readmatrix('resonance_6.txt'); % 19360x8
T1=readmatrix('transient.txt'); % 19360x8

for i=1:7
[freqT(:,i),YfftT(:,i)]=fftoperator(T1(:,1),T1(:,i+1));
end

% Impulse load -> Estimate eigenfrequencies
figure(1)
clf
plot(freqT,YfftT)
axis([0 30 0 0.06])
sgtitle('Impulse load')
xlabel('Frequency [Hz]')
wnexp=[2.48 7.49 12.35 16.32 19.52 21.59]; %Estimate damping
a=1./sqrt(2).*[0.0402 0.0203 0.0169 0.0281 0.0194 0.0492]; %Frequency divided by square root of two
ay=[];

%Half-power bandwidth plotting
xplot=[];
for j=1:6
xplot(:,j)=linspace(wnexp(j)-1,wnexp(j)+1,100);
ay(:,j)=a(j).*ones(length(xplot(:,j)),1);
end

% Plotting the heigth amplitude/sqrt(2) for all the natural frequencies
figure(1)
hold on
plot(xplot,ay,'Color','k')

% Half-power bandwidth calculation
w12=[2.46 2.54; 7.46 7.55; 12.31 12.41; 16.28 16.34; 19.48 19.55; 21.55 21.62]; %w1 and w2 from figure 1
zeta=zeros(1,6);
for j=1:6
zeta(j)=(w12(j,2)-w12(j,1))/(2*wnexp(j)); %Estimated damping
end

% Plotting vector
y1=linspace(0,6,7)';
x1=zeros(7,1);
y2=y1.*ones(7,7);
R=[10*R1(100,2:8); R2(100,2:8); 0.3*R3(100,2:8); 0.2*R4(100,2:8); 0.3*R5(100,2:8); 0.2*R6(100,2:8)];

% Harmonic excitation -> Mode shapes
figure(2)
clf
sgtitle('Mode shapes from experimental data')
for j=1:6
subplot(2,3,j)
plot(R(j,:),y1,0.7+R(j,:),y1,'Color','b')
hold on
plot(x1,y1,x1+0.7,y1,'LineStyle','--','Color','k')
plot([R(j,:); 0.7+R(j,:)],[y1 y1]','Color','k')
title('Mode shape', num2str(j))
axis([-0.2 0.9 0 6])
yticks([1 2 3 4 5 6])
end

%% Material data
% Steel
Es=210e9; %Young-modulus Pa
rho_s=8000; %Density kg/m^3
A1=2*15e-6; %Cross section area type 1 m^2
A2=3*10e-6; %Cross section area type 2 m^2
I1=15e-12*2^3/12; %Moment of inertia type 1 m^4
I2=10e-12*3^3/12; %Moment of inertia type 2 m^4
Atot=A1+A2; %Total area m^2
Itot=I1+I2; %Total moment of inertia m^4
Vs=Atot*(0.245+0.029); %Volume m^3
ms=rho_s*Vs; %Mass kg
eps=[Es Atot Itot rho_s*Atot];

% Wood
Ew=11e9; %Young-modulus Pa
rho_w=590; %Density kg/m^3
Aw=0.2*0.029; %Cross section area m^2
Iw=0.2*0.029^3/12; %Moment of inertia m^4
Vw=0.7*0.2*0.029; %Volume m^3
mw=rho_w*Vw; %Mass kg
epw=[Ew Aw Iw rho_w*Aw];

%% Shear building analysis
% Mass and stiffness matrix
M=mw*eye(6); %Assumption: only the mass of the wood is considered
l=0.245;
k=24*Es*Itot/l^3; %Spring constant
v=-1*ones(1,5);
K=diag(v,1)+diag(v,-1)+2*eye(6);
K(6,6)=1;
K=k*K; %Stiffness matrix

% Solving eigenproblem without damping
[L,X]=eigen(K,M);

% Eigenfrequencies
wn=sqrt(L) %Eigenfrequencies [rad/s]
fn=wn/(2*pi) %Eigenfrequencies [Hz]
% Plotting eigenmodes
x=0.5*[zeros(1,6);X];
x=[0.8*x(:,1) 0.5*x(:,2) 0.3*x(:,3) 0.3*x(:,4) 0.3*x(:,5) 0.3*x(:,6)];

figure(3)
clf
sgtitle('Mode shapes from shear building analysis')
for j=1:6
subplot(2,3,j)
plot(x(:,j),y1,0.7+x(:,j),y1,'Color','b')
hold on
plot(x1,y1,0.7+x1,y1,'LineStyle','--','Color','k')
plot([x(:,j) 0.7+x(:,j)]',[y1 y1]','Color','k')
title('Mode shape', num2str(j))
axis([-0.2 0.9 0 6])
yticks([1 2 3 4 5 6])
end

% Harmonic load (Earthquake forces)
ag=1;
p0=-mw*ag*ones(6,1); %Harmonic loading
w=linspace(0,150,10000);
u0_sba=zeros(6,length(w));

% Calculating Kd and u0 for 0<w<150 rad/s
for j=1:length(w)
Kd=K-w(j)^2*M;
u0_sba(:,j)=solveq(Kd,p0);
end

% Plot vectors
x2=linspace(-50,50,10);
for j=1:length(wn)
wnplot(:,j)=wn(j).*ones(length(x2),1);
end

% Plotting frequency sweep
figure(4)
clf
sgtitle('Resonance frequencies with earthquake forces')
for j=1:6
subplot(6,1,j)
plot(w,u0_sba(j,:),'Color','b')
hold on
plot(wnplot,x2,'LineStyle','-','Color','k','LineWidth',0.1)
ylabel(['Dof ', num2str(j)])
xticks(wn)
axis([0 140 -0.1 0.1])
grid on
end
xlabel('\omega [rad/s]')

%% FE analysis
% Topology
Edof=[1 1 2 3 7 8 9
2 4 5 6 10 11 12
3 7 8 9 13 14 15
4 10 11 12 16 17 18
5 13 14 15 19 20 21
6 16 17 18 22 23 24
7 19 20 21 25 26 27
8 22 23 24 28 29 30
9 25 26 27 31 32 33
10 28 29 30 34 35 36
11 31 32 33 37 38 39
12 34 35 36 40 41 42
13 10 11 12 7 8 9
14 16 17 18 13 14 15
15 22 23 24 19 20 21
16 28 29 30 25 26 27
17 34 35 36 31 32 33
18 40 41 42 37 38 39];

L=0.7;
h=0.245;

Coord=[L 0; 0 0; L h; 0 h; L 2*h; 0 2*h; L 3*h; 0 3*h; L 4*h; 0 4*h; L 5*h; 0 5*h; L 6*h; 0 6*h];
Dof=[1 2 3; 4 5 6; 7 8 9; 10 11 12; 13 14 15; 16 17 18; 19 20 21; 22 23 24; 25 26 27; 28 29 30; 31 32 33; 34 35 36; 37 38 39; 40 41 42];

[Ex,Ey]=coordxtr(Edof,Coord,Dof,2);

% Element matrices
K=zeros(42);
M=zeros(42);

for j=1:12
[k,m,c]=beam2de(Ex(j,:),Ey(j,:),eps);
K=assem(Edof(j,:),K,k);
M=assem(Edof(j,:),M,m);
end

for j=13:18
[k,m,c]=beam2de(Ex(j,:),Ey(j,:),epw);
K=assem(Edof(j,:),K,k);
M=assem(Edof(j,:),M,m);
end

%Plotting FE model
figure(5)
clf
eldraw2(Ex,Ey,[1 2 2],Edof);
sgtitle('2D Frame structure')
axis([-1 L+1 0 7*h])

% Eigenvalue analysis
b=[1 2 3 4 5 6]; %Boundary conditions
[L,X]=eigen(K,M,b);
X6= X(1:6:end,:); %Choosing relevant dofs
fn_FEM=sqrt(L)/(2*pi);
fn_FEM=fn_FEM(1:6)
wn_FEM=fn_FEM*2*pi

% Scaling the vector
X6=[0.5*X6(:,1) 0.4*X6(:,2) 0.3*X6(:,3) 0.2*X6(:,4) 0.2*X6(:,5) 0.2*X6(:,6)];

% Mode shapes
figure(6)
clf
sgtitle('Mode shapes from undamped FE analysis')
for j=1:6
subplot(2,3,j)
plot(X6(:,j),y1,0.7+X6(:,j),y1,'Color','b')
hold on
plot(x1,y1,0.7+x1,y1,'LineStyle','--','Color','k')
plot([X6(:,j) 0.7+X6(:,j)]',[y1 y1]','Color','k')
title('Mode shape', num2str(j))
axis([-0.2 0.9 0 6])
yticks([1 2 3 4 5 6])
end
u_undamped=[];

% Frequency sweep undamped
p=ones(3*length(Dof),1);
wu=linspace(0,150,1000);
bc=[b; zeros(1,length(b))]';

for j=1:length(wu)
Kud=K-M*(wu(j))^2;
uu=solveq(Kud,p,bc);
u_undamped(:,j) = uu;
end

u_undamped=u_undamped(1:6:end,:);

for j=1:length(wn)
wnplot(:,j)=wn_FEM(j).*ones(length(x2),1);
end

figure(7)
clf
sgtitle('Resonance frequencies without damping')
for j=1:6
subplot(6,1,j)
plot(wu,u_undamped(j+1,:),'Color','b')
hold on
plot(wnplot,x2,'LineStyle','-','Color','k','LineWidth',0.1)
ylabel(['Dof ', num2str(j)])
axis([0 140 -0.01 0.01])
xticks(wn_FEM)
end

xlabel('\omega [rad/s]')
nu=length(wu);

u1u=[min(u_undamped(2,0.08*nu:0.12*nu)), min(u_undamped(3,0.08*nu:0.12*nu)), min(u_undamped(4,0.08*nu:0.12*nu)), min(u_undamped(5,0.08*nu:0.12*nu)), min(u_undamped(6,0.08*nu:0.12*nu)), min(u_undamped(7,0.08*nu:0.12*nu))];
u2u=[min(u_undamped(2,0.29*nu:0.33*nu)), min(u_undamped(3,0.29*nu:0.33*nu)),min(u_undamped(4,0.29*nu:0.33*nu)), min(u_undamped(5,0.29*nu:0.33*nu)), max(u_undamped(6,0.29*nu:0.33*nu)), max(u_undamped(7,0.29*nu:0.33*nu))];
u3u=[min(u_undamped(2,0.48*nu:0.52*nu)), min(u_undamped(3,0.48*nu:0.52*nu)), max(u_undamped(4,0.48*nu:0.52*nu)), max(u_undamped(5,0.48*nu:0.52*nu)), max(u_undamped(6,0.48*nu:0.52*nu)), min(u_undamped(7,0.48*nu:0.52*nu))];
u4u=[min(u_undamped(2,0.64*nu:0.68*nu)), max(u_undamped(3,0.64*nu:0.68*nu)),max(u_undamped(4,0.64*nu:0.68*nu)), min(u_undamped(5,0.64*nu:0.68*nu)), min(u_undamped(6,0.64*nu:0.68*nu)), max(u_undamped(7,0.64*nu:0.68*nu))];
u5u=[min(u_undamped(2,0.76*nu:0.80*nu)), max(u_undamped(3,0.76*nu:0.80*nu)), min(u_undamped(4,0.76*nu:0.80*nu)), min(u_undamped(5,0.76*nu:0.80*nu)), max(u_undamped(6,0.76*nu:0.80*nu)), min(u_undamped(7,0.76*nu:0.80*nu))];
u6u=[min(u_undamped(2,0.84*nu:0.88*nu)), max(u_undamped(3,0.84*nu:0.88*nu)), min(u_undamped(4,0.84*nu:0.88*nu)), max(u_undamped(5,0.84*nu:0.88*nu)), min(u_undamped(6,0.84*nu:0.88*nu)), max(u_undamped(7,0.84*nu:0.88*nu))];

% Scaling the amplitudes for the plot
uplot=[0.1*u1u; 4*u2u; 10*u3u; 10*u4u; 20*u5u; 30*u6u];
uplot=[zeros(6,1) uplot]';
figure(8)
clf
sgtitle('Mode shapes without damping using frequency sweep')
for j=1:6
subplot(2,3,j)
plot(uplot(:,j),y1,0.7+uplot(:,j),y1,'Color','b')
hold on
plot(x1,y1,0.7+x1,y1,'LineStyle','--','Color','k')
plot([uplot(:,j) 0.7+uplot(:,j)]',[y1 y1]','Color','k')
title('Mode shape', num2str(j))
axis([-0.2 0.9 0 6])
yticks([1 2 3 4 5 6])
end

% Damping matrix
w=linspace(0,150,10000);
damp=mean(zeta');
for j=1:length(w)
C=2*damp*K/w(j); %Damping matrix
Kfem=K+1i*w(j)*C-M*(w(j))^2;
u=solveq(Kfem,p,bc);
u_damp(:,j) = u;
end
u6=imag(u_damp(1:6:end,:)); %Choosing relevant dofs
u6(:,1)=[0 0 0 0 0 0 0]'; %Creating "boundry conditions"

figure(9)
clf
sgtitle('Resonance frequencies with damping')
for j=1:6
subplot(6,1,j)
plot(w,u6(j+1,:),'Color','b')
hold on
plot(wnplot,x2,'LineStyle','-','Color','k','LineWidth',0.1)
ylabel(['Dof ', num2str(j)])
axis([0 140 -0.0005 0.0005])
xticks(wn_FEM)
end
xlabel('\omega [rad/s]')

% Picking extreme values for each natural frequency
n=length(w);
u1=[min(u6(2,0.08*n:0.12*n)) min(u6(3,0.08*n:0.12*n)) min(u6(4,0.08*n:0.12*n)) min(u6(5,0.08*n:0.12*n)) min(u6(6,0.08*n:0.12*n)) min(u6(7,0.08*n:0.12*n))];
u2=[min(u6(2,0.29*n:0.33*n)) min(u6(3,0.29*n:0.33*n)) min(u6(4,0.29*n:0.33*n)) min(u6(5,0.29*n:0.33*n)) max(u6(6,0.29*n:0.33*n)) max(u6(7,0.29*n:0.33*n))];
u3=[min(u6(2,0.48*n:0.52*n)) min(u6(3,0.48*n:0.52*n)) max(u6(4,0.48*n:0.52*n)) max(u6(5,0.48*n:0.52*n)) max(u6(6,0.48*n:0.52*n)) min(u6(7,0.48*n:0.52*n))];
u4=[min(u6(2,0.64*n:0.68*n)) max(u6(3,0.64*n:0.68*n)) max(u6(4,0.64*n:0.68*n)) min(u6(5,0.64*n:0.68*n)) min(u6(6,0.64*n:0.68*n)) max(u6(7,0.64*n:0.68*n))];
u5=[min(u6(2,0.76*n:0.80*n)) max(u6(3,0.76*n:0.80*n)) min(u6(4,0.76*n:0.80*n)) min(u6(5,0.76*n:0.80*n)) max(u6(6,0.76*n:0.80*n)) min(u6(7,0.76*n:0.80*n))];
u6=[min(u6(2,0.84*n:0.88*n)) max(u6(3,0.84*n:0.88*n)) min(u6(4,0.84*n:0.88*n)) max(u6(5,0.84*n:0.88*n)) min(u6(6,0.84*n:0.88*n)) max(u6(7,0.84*n:0.88*n))];

% Scaling the amplitudes for the plot
u=[0.3*u1; 4*u2; 20*u3; 50*u4; 100*u5; 300*u6];
u=[zeros(6,1) u];
figure(10)
clf
sgtitle('Mode shapes with damping using frequency sweep')
for j=1:6
subplot(2,3,j)
plot(u(j,:),y1,0.7+u(j,:),y1,'Color','b')
hold on
plot(x1,y1,0.7+x1,y1,'LineStyle','--','Color','k')
plot([u(j,:); 0.7+u(j,:)],[y1 y1]','Color','k')
title('Mode shape', num2str(j))
axis([-0.2 0.9 0 6])
yticks([1 2 3 4 5 6])
end