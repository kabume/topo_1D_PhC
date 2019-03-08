%zww²Î¿¼ÎÄÏ×25
%This program could draw the figure 3. h==>l, D==>Lambda
clc,clear
L=100;D=2*L;h0=L/2;
Delta=-0.5;
h=(1+Delta)*L;h1=L-h/2;h2=h;ne=2;c=3e8;
w=0:0.001:0.6;
k=w*1e7*ne/c;
Tc1=1+1j*tan(k*h0)/2;Tc2=1j*tan(k*h0)/2;Tc3=-1j*tan(k*h0)/2;Tc4=1-1j*tan(k*h0)/2;
num=length(w);
t=zeros(1,num);r=zeros(1,num);
for n=1:num
    Tc=[Tc1(n),Tc2(n);Tc3(n),Tc4(n)];
    Th1=[exp(1j*k(n)*h1),0;0,exp(-1j*k(n)*h1)];
    Th2=[exp(1j*k(n)*h2),0;0,exp(-1j*k(n)*h2)];
    T=Th1*Tc*Th2*Tc*Th1;
    yeta=acos(abs(T(1,1)))/1i;
    Z0=T(1,2)/sinh(yeta);
    if real(Z0)>0
        Z0=T(1,2)/sinh(yeta);
    else
        Z0=-T(1,2)/sinh(yeta);
    end
    n0=1i*yeta/D/k(n);
    if imag(n0)>0
        n0=1i*yeta/D/k(n);
    else
        n0=-1i*yeta/D/k(n);
    end
    epsilon(n)=n0/Z0;
    mu(n)=n0*Z0;
end
figure
plot(w,epsilon,'k');
hold on 
plot(w,mu);
%axis([0.1602 0.2316 -0.25 0.25]);