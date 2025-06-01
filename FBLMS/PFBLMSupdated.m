% Constrained PFBLMS Algorithm For Modeling(M=p*L)
tic
clear;clc;
runs=input("enter the number of the experiments:");
num=input("enter the number of the blocks:");
sigma2=input("enter the variance of the plant noise:");
misa=input("enter the adjustment:");
p=input("enter the parameter p:");
h=[0.1,-0.2,-0.3,0.4,0.4,-0.2,-0.1];
L=64; % the length of every block
wo=rand(1,1985);wo=wo/sqrt(sum(wo.^2)); % normalization of wo
N=length(wo); % the length of the filter
beta=0.95;
P=5;
miup0=misa/P; % the step size parameter
M=p*L;
MSEPFB=zeros(1,10000);
epsilon=0.1;
for j=1:runs
    v=randn(1,200000);
    xp=filter(h,1,v);
    dp=filter(wo,1,xp)+sqrt(sigma2)*randn(1,length(xp));
    xfm=zeros((P-1)*p+1,L+M);
    for i=1:(P-1)*p+1
        xfm(i,:)=fft(xp((i-1)*L+1:i*L+M));
    end
    wfm=zeros(P,L+M);    
    sigma2p=zeros(1,L+M); miup=zeros(1,L+M);
for k=1:num
% filtering
    t=zeros(1,L+M);
    for i=1:P
        t=t+wfm(i,:).*xfm((P-1)*p+1+p-p*i,:);
    end
yp=ifft(t);yp=yp(end-L+1:end);
% error estimation    
d=dp((P-1)*L+1+M+(k-1)*L:P*L+M+(k-1)*L);
e=d-yp;MSEPFB(k)=MSEPFB(k)+sum(e.^2)/L;
% step normalization
    for i=1:L+M
        sigma2p(i)=beta*sigma2p(i)+(1-beta)*abs(xfm(end,i))^2;
        miup(i)=miup0/(sigma2p(i)+epsilon); % add epsilon to avoid that the number is too close to 0
    end
% adaptation of the tap-weight vector & constraint of the tap-weight vector
epf=fft([zeros(1,M),e]);
    for i=1:P
        wfm(i,:)=wfm(i,:)+2*miup.*conj(xfm((P-1)*p+1+p-p*i,:)).*epf;
        %{
        iwf=ifft(wfm(i,:));
        iwf=iwf(1:M);
        wfm(i,:)=fft([iwf,zeros(1,L)]);
        %}
    end
% update xfm
    for i=1:(P-1)*p
        xfm(i,:)=xfm(i+1,:);
    end
xfm(end,:)=fft(xp((P-1)*L+1+k*L:P*L+M+k*L));
end
end
MSEPFB=MSEPFB/runs;
semilogy(MSEPFB);
hold on;
xlabel("the number of the blocks");ylabel("MSE");
title("Learning Curve");
hold off;
toc