clear;clc;
tic
n=1:30000;
a=sqrt(2);sigma2=0.001;
w0=input("enter the desired frequency:"); 
N=input("enter the length of the filter:");
R=toeplitz((a^2/2)*cos(w0*(0:N-1)))+sigma2*eye(N);% generate the correlation matrix eye:generate a unit matrix
M=input("enter the misadjustment:");
miu=M/trace(R);
MSE=zeros(1,30000);
ytest=zeros(1,30000);
for j=1:100
    w=zeros(1,N);
    x=a*sin(w0*n+2*pi*randn)+sqrt(sigma2)*randn(1,length(n));
for i=1:20000
    x1=x(i+N-1:-1:i);
    d=x(i+N);
    y=w*x1';
    ytest(i)=y;
    e=d-y;MSE(i)=MSE(i)+e^2;
    w=w+2*miu*e*x1;
end
end
MSE=MSE/100;
semilogy(MSE);
xlabel("the number of iterations");ylabel("MSE");
toc
