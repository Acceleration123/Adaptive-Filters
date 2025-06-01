% IIR Adaptive Line Enhancement
clear;clc
tic
s=input("enter the parameter s:");
w0=input("enter the parameter w:");
sigma2=input("enter the variance of v(n):");
a=input("enter the constant a:");
k=0:10000000;w=k*2*pi/10000000;
W=(1-s)*(w0-exp(-1i*w))./(1-(1+s)*w0*exp(-1i*w)+s*exp(-2*1i*w));
subplot(3,1,1);
plot(w,abs(W));grid on;
title("|W(exp(j\omega))|");
% the performance function I
subplot(3,1,2);
wo=-1:0.01:1;
W1=(1-s)*(wo-exp(-1i*acos(w0)))./(1-(1+s)*wo*exp(-1i*acos(w0))+s*exp(-2*1i*acos(w0)));
pf1=2*sigma2/(1+s)+(a^2/2)*abs(1-exp(-1i*acos(w0))*W1).^2;
plot(wo,pf1,'r');grid on;
title('The Performance Function I');
% the performance function II for adjusting s
subplot(3,1,3);
pf2=(1+s)/(1-s)*a^2/2*abs(W1).^2+sigma2;
plot(wo,pf2,'b');grid on;
title('The performance Function II');
toc
% IIR ALE(algorithm 1)
clear;clc;
tic
a=sqrt(2);
s=0.75;w0=0;
sigma2=1;
miuw=0.025;mius=0.0005;% the step-size parameters 
theta=input('enter the desired frequency:');
n=0:30000;
x1=zeros(1,length(n));
thetaj=-theta;
for j=0:length(n)-1
thetaj=thetaj+theta;
if thetaj>2*pi
    thetaj=thetaj-2*pi;
end
x1(j+1)=a*sin(thetaj)+wgn(1,1,sigma2);
end
x1(1:3)=0;
y=[0,0,0];
alpha=y;beta=y;
wtest=zeros(1,30000);
for i=3:25000
y(3)=y(2);
y(2)=y(1);
y(1)=(1+s)*w0*y(2)-s*y(3)+(1-s)*(w0*x1(i)-x1(i-1));
e=x1(i+1)-y(1);
alpha(1)=(1+s)*w0*alpha(2)-s*alpha(3)+(1+s)*y(2)+(1-s)*x1(i);
beta(1)=(1+s)*w0*beta(2)-s*beta(3)-(w0*(x1(i)-y(2))-(x1(i-1)-y(3)));
alpha(2:3)=alpha(1:2);
beta(2:3)=beta(1:2);
w0=w0+2*miuw*(1-s)^3*e*alpha(1);
if abs(w0)>0.999
		w0=0.999*sign(w0);
end
s=s+(2*mius/(1-s)^2)*(y(1)*y(1)+(1-s*s)*y(1)*beta(1));
if s<0.25
    s=0.25;
elseif s>0.9
    s=0.9;
end
wtest(i-2)=w0;
end
plot(wtest);
toc
% IIR ALE(algorithm 2)
clear;clc;
tic
a=sqrt(2);w0=0;
s=0.75;theta0=0;
sigma2=0.5;
mius=0.0005;% the step-size parameters 
theta=input('enter the desired frequency:');
n=0:30000;
x1=a*sin(theta*n)+wgn(1,length(n),sigma2);
x1(1:3)=0;
y=[0,0,0];
alpha=y;beta=y;
thetatest=zeros(1,30000);
for i=3:25000
w0=cos(theta0);
y(3)=y(2);
y(2)=y(1);
y(1)=(1+s)*w0*y(2)-s*y(3)+(1-s)*(w0*x1(i)-x1(i-1));
miutheta=0.05*(1-s)^3;    
e=x1(i+1)-y(1);
alpha(1)=(1+s)*w0*alpha(2)-s*alpha(3)-sin(theta0)*((1+s)*y(2)+(1-s)*x1(i));
beta(1)=(1+s)*w0*beta(2)-s*beta(3)-(w0*(x1(i)-y(2))-(x1(i-1)-y(3)));
alpha(2:3)=alpha(1:2);
beta(2:3)=beta(1:2);
theta0=theta0+2*miutheta*e*alpha(1);
if theta0<0.05
	theta0=0.05;
end
s=s+2*mius*((y(1))^2/(1-s)^2+((1+s)/(1-s))*y(1)*beta(1));
if s<0.25
    s=0.25;
elseif s>0.9
    s=0.9;
end
wtest(i-2)=w0;
end
plot(wtest);
toc