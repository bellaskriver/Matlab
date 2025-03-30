function i5u1
close all

%alla decimaler
format long
exact=0.6576698563283957;
a=0;
b=0.8;
integr=[]
relfel=[]
nvec=[]

for m=2:22;
n=2.^m+1;
h=(b-a)/(n-1);
xvec=linspace(a,b,n)';
gvec=g(xvec);
w=h*[0.5; ones(n-2,1); 0.5];
I=gvec'*w
relerr=abs(exact-I)/abs(exact)
mlI=integral(@g,a,b)
mlIerr=abs(exact-mlI)/abs(exact)
integr=[integr;I]
relfel=[relfel;relerr]
nvec=[nvec;n]
end

tic
gvec=g(xvec);
I=gvec'*w;
t1=toc
tic
mlI=integral(@g,a,b);
mlIerr=abs(exact-mlI)/abs(exact);
t2=toc
nogrannhet=relerr-mlIerr

%log-log graf
loglog(nvec,relfel,'*')
title('Log-log plot of function av n och relerr')
xlabel('nvec')
ylabel('relfel')
grid on

function gout=g(x)
gout=exp(-x.^2);

%Problem 2
function i5u2
a=9.38
r=[9.38 9.90 10.42 10.94 11.46 11.98 12.50 13.02 13.54 14.06 14.58]'
b=14.58
n=length(r)
h=(b-a)/(n-1)
w=h*[0.5; ones(n-2,1); 0.5]
fvec=f(r)
gvec=g(r)
I1=fvec'*w
I2=gvec'*w
Tmed=I1/I2
wi=ones(n,1);
wi(2:2:10)=4;
wi(3:2:9)=2;
wi=h/3*wi;
I=fvec'*wi
Tmed1=I/I2
skillnad=Tmed-Tmed1

function fout=f(x)
T=[338 423 474 506 557 573 601 622 651 661 671]'
fout=T.*x;

function gout=g(x)
gout=x;