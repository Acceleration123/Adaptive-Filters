% Recursive Least-Squares For Modeling
clear;clc;
tic
sigma2=input("enter the variance of the noise:");
Mis=input("enter the misadjustment:");
delta=0.0001;
n=1:30000;
h1=[0.35,1,-0.35];
h2=[0.35,1,0.35];
wo=[ones(1,8),-1*ones(1,7)]';
N=length(wo);
lambda=(1-Mis/N)/(1+Mis/N);
MSE1=zeros(1,30000);MSE2=zeros(1,30000);
for j=1:100
x1=zeros(N,30000);x2=zeros(N,30000);% different colunm for different time
for m=1:N
v=randn(1,length(n));
x1(m,:)=filter(h1,1,v);
x2(m,:)=filter(h2,1,v);
end
psiinv1=(1/delta)*eye(N);
psiinv2=(1/delta)*eye(N);
w1=zeros(N,1);w2=zeros(N,1);
for i=1:300
x1tdl=x1(:,i);
x2tdl=x2(:,i);
u1=psiinv1*x1tdl; 
u2=psiinv2*x2tdl; 
k1=u1/(lambda+x1tdl'*u1); % the calculation of the gain vector
k2=u2/(lambda+x2tdl'*u2);
y1=w1'*x1tdl;y2=w2'*x2tdl;
d1=wo'*x1tdl+sqrt(sigma2)*randn(1,1);   
d2=wo'*x2tdl+sqrt(sigma2)*randn(1,1);
e1=d1-y1;e2=d2-y2;
MSE1(i)=MSE1(i)+e1^2;
MSE2(i)=MSE2(i)+e2^2;
w1=w1+e1*k1;w2=w2+e2*k2;
psiinv1=(1/lambda)*(psiinv1-k1*(x1tdl'*psiinv1));
psiinv2=(1/lambda)*(psiinv2-k2*(x2tdl'*psiinv2));
end
end
MSE1=MSE1/100;
MSE2=MSE2/100;
figure(1);
semilogy(MSE1);
xlabel("the number of the iterations");ylabel("MSE");
title("Learning Curve");
hold on;
semilogy(MSE2);
hold off;
toc