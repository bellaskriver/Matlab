% Hand-in 1 VSMN10 - Isabella McAuley Skriver

close all; clear all; format short; format compact
d=[0.01; 0.1; 0.2; 0.7; 1]; %Dampening ratio
t=linspace(0,3,100)'; %omega/omega_n are the values along the x-axis

for k=1:5 %Looping over all dampening ratios
Rd(:,k)=abs(1./(1-(t).^2+1i*2.*d(k).*t)); %Complex deformation response factor
phi(:,k)=rad2deg(angle(1-(t.^2)+1i*2*d(k).*t)); %Phase-shift

figure(1) %Complex deformation response factor
plot(t, Rd, "LineWidth",2)
xlabel('Frequency ratio \omega / \omega_{n}')
ylabel('Deformation response factor, \it{R_d = u_0 / (u_{st})_0}')
grid on
axis([0 3 0 5])
title('Displacement amplitude')

figure(2) %Phase-shift
plot(t, phi,"LineWidth",2)
xlabel('Frequency ratio \omega / \omega_{n}')
ylabel('Phase angle \phi')
grid on
axis([0 3 0 180])
title('Phase-angle')
end

figure(1) %Complex deformation response factor
legend('\zeta=0.01','\zeta=0.1','\zeta=0.2','\zeta=0.7','\zeta=1')

figure(2) %Phase-shift
legend('\zeta=0.01','\zeta=0.1','\zeta=0.2','\zeta=0.7','\zeta=1')

%Calculation for 12 Hz and 3% dampening ratio
L=4; %Beam length in m
mm=750; %Motor mass in kg
mb=16; %Beam mass in kg/m
EI=1.5e6; %Stiffnes multiplied by moment of inertia in Nm^2
e=0.002; %Excentricity in mm
mr=250; %Rotating mass in kg
f=12; %Frequency in Hz
k=2*48*EI/L^3; %Beam stiffness from table (simply supported), se example in lecture 2
mtotal=2*0.5*mb*L+mm; %Total mass
wn=sqrt(k/mtotal); %Natural angular frequency
w=f*2*pi; %Loading angular frequency
p0=mr*e*w^2; %Force
u0=(p0/k)*abs(1/((1-(w/wn)^2)+1i*2*0.03*(w/wn))); %Vibration amplitude, dampening ratio 3%=0.003
disp(['The vibration amplitude is ', num2str(u0), ' m, which is equivalent to 1.2 mm'])