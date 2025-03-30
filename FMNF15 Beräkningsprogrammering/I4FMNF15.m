function main
p2i4();
p4i4();

function p2i4
close all

%Alt 1
P=[10.00; 15.00; 20.00; 25.00;];
A=[ones(length(P),1) P];
L=[ 11.60; 11.85; 12.25; 12.75;];
c=A\L

%Alt 2
P=linspace(10, 25, 4);
L=[11.60 11.85 12.25 12.75];
p=polyfit(P,L,1)

function p4i4
close all

%F(t)=c1+c2*t+c3*t.^2
F1=[13.18; 15.78; 17.97; 18.38; 16.63; 14.07;];
t1=[91; 121; 152; 182; 213; 244;];
A1=[t1.^0 t1 t1.^2];
c1=A1\F1

% 6 juni = dag 157
d1=[157];
cm1=c1;
Langd1=cm1(1)+cm1(2)*d1+cm1(3)*d1^2

% F(t)=c1+c2*cos(w*t2)+c3*sin(w*t2)
F2=[6.13; 8.02; 10.42; F1; 11.43; 8.73; 6.55;];
t2=[1; 32; 60; t1; 274; 305; 335;];
w=(2*pi)/365;
A2=[t2.^0 cos(w*t2) sin(w*t2)];
c2=A2\F2

% 6 juni = dag 157
d1=[157];
cm2=c2;
Langd2=cm2(1)+cm2(2)*cos(w*d1)+cm2(3)*sin(w*d1)
residual=norm(F2-A2*c2)

plot (t2,F2,'*', t2, F2,'-b')
xlabel('Dagnummer')
ylabel('Dagens längd i timmar')
title ('Dagslängd 1:a varje månad')
grid on