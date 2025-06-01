clear;clc;
tic
misa=input("enter the misadjustment:");
MSE1=zeros(1,20000);
MSE2=zeros(1,20000);
N=15;% the length of the filter
wo=[ones(1,7),-1*ones(1,8)];% the impulse response of the filter

h1=[0.35,1,-0.35];nh1=0:length(h1)-1;
h2=[0.35,1,0.35];nh2=0:length(h2)-1;% the impulse responses of coloring filters
k=0:1000;w=k*2*pi/1000;
H1=h1*exp(-1i*nh1'*w);H2=h2*exp(-1i*nh2'*w);
figure(1);
plot(w,(abs(H1)).^2,'r');grid on;% .^2! not ^2 directly!
hold on;
plot(w,abs((H2)).^2,'g');% the power spectral density of the two input sequences 
title('the power spectral density of the two input sequences ');

for m=1:100
v=randn(1,20000);nv=0:length(v)-1;% generate the Gaussian white noise sequence
x1=filter(h1,1,v);nx1=0:length(x1)-1;
x2=filter(h2,1,v);nx2=0:length(x2)-1;
R1=corlnm2(x1,N);% calculate the correlation matrix of the input x
R2=corlnm2(x2,N);
lambda1=eig(R1);lambda2=eig(R2);
spread1=max(lambda1)/min(lambda1);
spread2=max(lambda2)/min(lambda2);
tr1=trace(R1);
tr2=trace(R2);
miu1=misa/(N*(h1*h1'));miu2=misa/(N*(h2*h2'));
x11=zeros(1,N);x22=zeros(1,N);
w1=zeros(1,N);w2=zeros(1,N);
for j=0:10000
x11=x1(N+j:-1:j+1);
x22=x2(N+j:-1:j+1);
eo=sqrt(0.001)*randn(1,1);
d1=wo*x11'+eo;d2=wo*x22'+eo;% generate the desired output
y1=w1*x11';
y2=w2*x22';
e1=d1-y1;MSE1(j+1)=MSE1(j+1)+e1^2;
e2=d2-y2;MSE2(j+1)=MSE2(j+1)+e2^2;
w1=w1+2*miu1*e1*x11;
w2=w2+2*miu2*e2*x22;
end
end

disp(w1);
disp(w2);
MSE1=MSE1/100;MSE2=MSE2/100;
figure(2);
semilogy(MSE1);
hold on;
semilogy(MSE2);
legend(["H1(z)","H2(z)"]);
toc
