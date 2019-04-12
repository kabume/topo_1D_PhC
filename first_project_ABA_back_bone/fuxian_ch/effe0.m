%半无限长的有效epsilon和mu,Z就是有效阻抗，ne必须是复数才能画出论文中的图
%This program could draw the figure 3. h==>l, D==>Lambda
clc,clear
L=100;D=200;h0=L/2;
Delta=-0.5;
h=(1+Delta)*L;h1=L-h/2;h2=h;ne=2;c=3e8;
w=0:0.001:0.6;
k=w*1e7*ne/c;
Tc1=1+1j*tan(k*h0)/2;Tc2=1j*tan(k*h0)/2;Tc3=-1j*tan(k*h0)/2;Tc4=1-1j*tan(k*h0)/2;
Z0=1/ne;
P=[1,1;1/Z0,-1/Z0];
num=length(w);
t=zeros(1,num);r=zeros(1,num);epsilon=zeros(1,num);mu=zeros(1,num);
for n=1:num
    Tc=[Tc1(n),Tc2(n);Tc3(n),Tc4(n)];
    Th1=[exp(1j*k(n)*h1),0;0,exp(-1j*k(n)*h1)];
    Th2=[exp(1j*k(n)*h2),0;0,exp(-1j*k(n)*h2)];
    T=Th1*Tc*Th2*Tc*Th1;
    M=P*T*P^-1;
    %t(n)=abs(1./T(2,2));r(n)=abs(-T(2,1)./T(2,2));
    Z=(M(1,1)-M(2,2)-sqrt((M(1,1)-M(2,2))^2+4*M(2,1)*M(1,2)))/(2*M(2,1));
    if real(Z)>0
         Z=(M(1,1)-M(2,2)-sqrt((M(1,1)-M(2,2))^2+4*M(2,1)*M(1,2)))/(2*M(2,1));
    else
        Z=(M(1,1)-M(2,2)+sqrt((M(1,1)-M(2,2))^2+4*M(2,1)*M(1,2)))/(2*M(2,1));
    end
%     r=abs((Z-Z0)/(Z+Z0));
%     t=1-r;
%     t0=exp(1i*k(n)*D)*t;
%     n0=acos(((1-(r^2-t0^2))/(2*t0))/(k(n)*D));
%     if imag(n0)>0
%         n0=acos(abs((1-(r^2-t0^2))/(2*t0))/(k(n)*D));
%     else
%         n0=-acos(abs((1-(r^2-t0^2))/(2*t0))/(k(n)*D));
%     end
%     z0=sqrt(((1+r)^2-t0^2)/((1-r)^2-t0^2));
%     if real(z0)>0
%         z0=sqrt(((1+r)^2-t0^2)/((1-r)^2-t0^2));
%     else
%         z0=-sqrt(((1+r)^2-t0^2)/((1-r)^2-t0^2));
%     end
    epsilon(n)=(2+2i)/Z;
    mu(n)=(2+2i)*Z;
end
figure
plot(w,epsilon,'k');
hold on 
plot(w,mu);
axis([0 0.6 -10 10])
grid on