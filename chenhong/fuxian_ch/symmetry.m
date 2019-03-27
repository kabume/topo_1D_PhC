%随频率变化的透反射图
%This program could draw the figure 3. h==>l, D==>Lambda
L=100;D=2*L;h0=L/2;
Delta=-0.8;
h=(1+Delta)*L;h1=L-h/2;h2=h;ne=2;c=3e8;Z0=1/ne;
w=0:0.001:0.6;
k=w*1e7*ne/c;
P=[1,1;1/Z0,-1/Z0];
Tc1=1+1j*tan(k*h0)/2;Tc2=1j*tan(k*h0)/2;Tc3=-1j*tan(k*h0)/2;Tc4=1-1j*tan(k*h0)/2;
num=length(w);
t=zeros(1,num);r=zeros(1,num);
for n=1:num
    Tc=[Tc1(n),Tc2(n);Tc3(n),Tc4(n)];
    Th1=[exp(1j*k(n)*h1),0;0,exp(-1j*k(n)*h1)];
    Th2=[exp(1j*k(n)*h2),0;0,exp(-1j*k(n)*h2)];
    T=Th1*Tc*Th2*Tc*Th1;
    M=P*T*P^-1;
    Z=(M(1,1)-M(2,2)-sqrt((M(1,1)-M(2,2))^2+4*M(2,1)*M(1,2)))/(2*M(2,1));
    if real(Z)>0
         Z=(M(1,1)-M(2,2)-sqrt((M(1,1)-M(2,2))^2+4*M(2,1)*M(1,2)))/(2*M(2,1));
    else
        Z=(M(1,1)-M(2,2)+sqrt((M(1,1)-M(2,2))^2+4*M(2,1)*M(1,2)))/(2*M(2,1));
    end
    %t(n)=abs(1./T(2,2));r(n)=abs(-T(2,1)./T(2,2));
    r(n)=abs((Z-Z0)/(Z+Z0));
end
figure
plot(w,r,'k');