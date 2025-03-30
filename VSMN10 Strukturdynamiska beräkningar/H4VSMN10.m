% Hand-in 4 VSMN10 - Isabella McAuley Skriver

clear all; close all;

%The values of k, p0, m don't affect the shape. They are set to 1.
k=1;
p0=1;
m=1;
syms w %Creates variables

%Natural eigenfrequencies
wn1=sqrt(k/m);
wn2=sqrt(3*k/m);

%Part a
u1a=k*p0/(m^2*(wn1^2-w^2)*(wn2^2-w^2));
u2a=(2*k-w^2*m)*p0/(m^2*(wn1^2-w^2)*(wn2^2-w^2));

%Part b
u1b=p0/(2*m)*((1/(wn1^2-w^2)-1/(wn2^2-w^2)));
u2b=p0/(2*m)*(1/(wn1^2-w^2)+(1/(wn2^2-w^2)));

%Plot, use fplot and syms to plot the functions
figure;
subplot(2,1,1);
fplot(u1a, [-2, 2],'blue','LineWidth',5);
hold on;
fplot(u1b, [-2, 2],'red','LineWidth',2);
hold off;
title('Part a');
legend('u1a', 'u1b');
grid on;
subplot(2,1,2);
fplot(u2a, [-2, 2],'blue','LineWidth',5);
hold on;
fplot(u2b, [-2, 2],'red','LineWidth',2);
hold off;
title('Part b');
legend('u2a', 'u2b');
grid on;