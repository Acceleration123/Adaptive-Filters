clear;clc;
tic
sigma2=1e-3;
h1=[0.35,1,-0.35];nh1=0:length(h1)-1;
h2=[0.35,1,0.35];nh2=0:length(h2)-1;% the impulse responses of coloring filters
n=1:30000;
N=input("enter the length of the euqalizer:");
delta=input("enter the delay:");
MSE1=zeros(1,30000);
MSE2=zeros(1,30000);
for j=1:100
    s=sign(randn(1,length(n)));% use sign to generate the sequence s(n)
    x1=filter(h1,1,s)+sqrt(sigma2)*randn(1,length(s));
    x2=filter(h2,1,s)+sqrt(sigma2)*randn(1,length(s));
    R1=corlnm2(x1,N);
    R2=corlnm2(x2,N);
    miu1=1500/trace(R1);
    miu2=2000/trace(R2);
    w1=zeros(1,N);
    w2=zeros(1,N);
for i=N+1:20000+N
    d=s(i);
    x11=x1(i+delta:-1:i+delta-N+1);
    x22=x2(i+delta:-1:i+delta-N+1);
    y1=w1*x11';
    y2=w2*x22';
    e1=d-y1;MSE1(i)=MSE1(i)+e1^2;
    e2=d-y2;MSE2(i)=MSE2(i)+e2^2;
    w1=w1+2*miu1*e1*x11;
    w2=w2+2*miu2*e2*x22;
end
end
MSE1=MSE1/100;
MSE2=MSE2/100;
semilogy(MSE1);hold on;
semilogy(MSE2);
xlabel(" the number of iterations");ylabel("MSE");
legend(["H1(z)","H2(z)"]);
toc
