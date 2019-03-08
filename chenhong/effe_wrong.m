%单层介质的有效mu和epsilon
%This program could draw the figure 3. h==>l, D==>Lambda
clc,clear
L=100;D=2*L;h0=L/2;
Delta=-0.5;
h=(1+Delta)*L;h1=L-h/2;h2=h;ne=2;c=3e8;
w=0:0.001:0.6;
k=w*1e7*ne/c;
Tc1=1+1j*tan(k*h0)/2;Tc2=1j*tan(k*h0)/2;Tc3=-1j*tan(k*h0)/2;Tc4=1-1j*tan(k*h0)/2;
num=length(w);
t=zeros(1,num);r=zeros(1,num);epsilon=zeros(1,num);mu=zeros(1,num);
for n=1:num
    Tc=[Tc1(n),Tc2(n);Tc3(n),Tc4(n)];
    Th1=[exp(1j*k(n)*h1),0;0,exp(-1j*k(n)*h1)];
    Th2=[exp(1j*k(n)*h2),0;0,exp(-1j*k(n)*h2)];
    T=Th1*Tc*Th2*Tc*Th1;
    t=abs(1/T(2,2));r=abs(-T(2,1)/T(2,2));
    t0=exp(1i*k(n)*D)*t;
    n0=acos(abs((1-(r^2-t0^2))/(2*t0))/(k(n)*D));
    if imag(n0)>0
        n0=acos((1-(r^2-t0^2))/(2*t0))/(k(n)*D);
    else
        n0=-acos((1-(r^2-t0^2))/(2*t0))/(k(n)*D);
    end
    z0=sqrt(((1+r)^2-t0^2)/((1-r)^2-t0^2));
    if real(z0)>0
        z0=sqrt(((1+r)^2-t0^2)/((1-r)^2-t0^2));
    else
        z0=-sqrt(((1+r)^2-t0^2)/((1-r)^2-t0^2));
    end
    epsilon(n)=n0/z0;
    mu(n)=n0*z0;
end
%figure
plot(w,epsilon,'k');
hold on 
plot(w,mu,'g');
axis([0 0.6 -4 4])