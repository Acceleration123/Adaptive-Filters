tic
clear;clc;
runs=input("enter the number of the experiments:");
num=input("enter the number of the blocks:");
sigma2=input("enter the variance of the plant noise:");
misa=input("enter the adjustment:");
h=[0.1,-0.2,-0.3,0.4,0.4,-0.2,-0.1];
L=64; % the length of every block
wo=rand(1,1985);wo=wo/sqrt(sum(wo.^2)); % normalization of wo
N=length(wo); % the length of the filter
beta=0.95;
miu0=misa*(N+L-1)/N;% the step size parameter
MSEFB=zeros(1,10000);
MSEFBv=zeros(1,10000);
for j=1:runs
    v=randn(1,200000);
    x=filter(h,1,v);
    d=filter(wo,1,x)+sqrt(sigma2)*randn(1,length(x));
    dv=filter(wo,1,v)+sqrt(sigma2)*randn(1,length(v));
    w=zeros(1,N);w1=zeros(1,N);
    sigma2es=zeros(1,N+L-1);sigma2esv=zeros(1,N+L-1);
    miu=zeros(1,N+L-1);miuv=zeros(1,N+L-1);
    wf=fft([w,zeros(1,L-1)]);w1f=fft([w1,zeros(1,L-1)]); % the tap-weight vector needs to be extended
for k=1:num
    xex=x((k-1)*L+1:k*L+N-1);vex=v((k-1)*L+1:k*L+N-1);
    dex=d((k-1)*L+N:k*L+N-1);dvex=d((k-1)*L+N:k*L+N-1);
% filtering    
    xf=fft(xex);vf=fft(vex);
    y=ifft(xf.*wf);y=y(end-L+1:end);
    yv=ifft(vf.*w1f);yv=yv(end-L+1:end);
% error estimation
    e=dex-y;MSEFB(k)=MSEFB(k)+mean(e.^2);
    ev=dvex-yv;MSEFBv(k)=MSEFBv(k)+mean(ev.^2);
% step normalization
    for i=1:N+L-1
        sigma2es(i)=beta*sigma2es(i)+(1-beta)*abs(xf(i))^2;
        sigma2esv(i)=beta*sigma2esv(i)+(1-beta)*abs(vf(i))^2;
        miu(i)=miu0/sigma2es(i);
        miuv(i)=miu0/sigma2esv(i);
    end
% adaptation of the tap-weight vector
    ef=fft([zeros(1,N-1),e]);evf=fft([zeros(1,N-1),ev]);
    wf=wf+2*miu.*conj(xf).*ef;w1f=w1f+2*miuv.*conj(vf).*evf;
% constraint of the tap-weight vector
    iwf=ifft(wf);iw1f=ifft(w1f);
    wf=fft([iwf(1:N),zeros(1,L-1)]);
    w1f=fft([iw1f(1:N),zeros(1,L-1)]);
end
end
MSEFB=MSEFB/runs;
MSEFBv=MSEFBv/runs;
semilogy(MSEFB);
hold on;
semilogy(MSEFBv);
xlabel("the number of the blocks");ylabel("MSE");
title("Learning Curve");
legend(["colored","white"]);
hold off;
toc