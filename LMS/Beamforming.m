%% Beamforming(LMS algorithm)
clear;clc;
tic
n=1:30000;
sigma2_alpha=input("enter the variance of the sequence alpha:");
sigma2_beta=input("enter the variance of the sequence beta:");
M=input("enter the misadjustment:");
theta=input("enter the angle:");
w0=input("enter the frequency:");
phi=pi*sin(theta);
R=((sigma2_alpha+sigma2_beta)/2)*eye(2);
miu=M/trace(R);
MSE=zeros(1,30000);
k=0:10000;theta0=k*2*pi/10000;
G=zeros(1,length(k));
for j=1:100
    alpha=sqrt(sigma2_alpha)*randn(1,length(n));
    beta=sqrt(sigma2_beta)*randn(1,length(n));
    phi1=2*pi*rand(1,length(n));
    phi2=2*pi*rand(1,length(n));
    sa=alpha.*cos(w0*n+phi1)+beta.*cos(w0*n+phi2-phi);
    x1=alpha.*cos(w0*n+phi1)+beta.*cos(w0*n+phi2);
    x2=alpha.*sin(w0*n+phi1)+beta.*sin(w0*n+phi2);
    w=zeros(1,2);
    for i=1:20000
        d=sa(i);
        xinput=[x1(i),x2(i)];
        y=w*xinput';
        e=d-y;MSE(i)=MSE(i)+e^2;
        w=w+2*miu*e*xinput;
    end
    G=G+(cos(pi*sin(theta0))-w(1)).^2+(sin(pi*sin(theta0))-w(2)).^2;
end
MSE=MSE/100;
subplot(2,1,1);
semilogy(MSE);
xlabel("the number of iterations");ylabel("MSE");
G=G/100;
subplot(2,1,2);
polarplot(theta0,G);
toc