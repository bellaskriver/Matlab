% Hand-in 2 VSMN10 - Isabella McAuley Skriver

close all; clear all; format short; format compact;

% Part a)
k=100000; %Stiffness N/m
m=1000; %Mass kg
K=k*[2 -1 0; -1 2 -1; 0 -1 1]; %Stiffness matrix
M=m*[1 0 0; 0 1 0; 0 0 1]; %Mass matrix
[L,X]=eigen(K,M); %Solve for eigenvalues and eigenmodes
L; %Eigenvalues for nodes 1,2,3
X; %Eigenmodes 1,2,3 for nodes 1,2,3, first coloum is mode 1
w=sqrt(L) %Natural frequency (circular) in rad/s for nodes 1,2,3
fn=w./(2*pi) %Natural frequency (cyclic) in Hz for nodes 1,2,3

%For plotting a straight line to compare with the eigenmodes
y=[0;1;2;3];
x=zeros(4,1);

%Illustrating the eigenmodes
for i=1:3
figure(1)
subplot(1,3,i)
S=[0 0 0; X];
plot(S(:,i),y,'-b', "LineWidth",2 ) %Mode shapes of the structure
hold on
plot(S(:,i),y, '*r',"LineWidth",2) %DOF of the structure
hold on
subplot(1,3,i)
plot(x,y,'--b',"LineWidth",2) %Shape of the initial structure
axis([-0.025 0.025 0 3.5])
end

subplot(1,3,1)
title('Eigenmode 1')
subplot(1,3,2)
title('Eigenmode 2')
subplot(1,3,3)
title('Eigenmode 3')

% Del b)
%For the harmonic force p0sinwt, sinwt can be disregarded
%This dynamic system is a time independent system of equations
p0=1; %Force amplitude in N
f=p0*ones(3,1); %Force vector
w=linspace(0, 20, 1000); %Vector omega 0<w<20
for i=1:length(w)
K1=K-w(i)^2.*M; %From the time independant equation
u0(:,i)=solveq(K1,f); %Displacement/steady state response
end

%Illustrating the steady state response
figure(2)
subplot(3,1,1)
plot(w,u0(1,:), "LineWidth",2)
title('Steady state respons for DOF 1')
xlabel('Eigenfrequencies \omega')
ylabel('Displacement {u_{1,0}}')

subplot(3,1,2)
plot(w,u0(2,:),"LineWidth",2)
title('Steady state respons for DOF 2')
xlabel('Eigenfrequencies \omega')
ylabel('Displacement {u_{2,0}}')

subplot(3,1,3)
plot(w,u0(3,:),"LineWidth",2)
title('Steady state respons for DOF 3')
xlabel('Eigenfrequencies \omega')
ylabel('Displacement {u_{3,0}}')