%奇对称和偶对称的场图，参考CT2014PRX中的附录2，可通过场分布图判断zak phase
%This program could draw the figure 3. h==>l, D==>Lambda
clc,clear
L=100;D=2*L;h0=L/2;
Delta=-0.8; %此时边界态的频率为0.1602，0.2316，0.3322，0.3768，0.4052
h=(1+Delta)*L;h1=L-h/2;h2=h;ne=2;c=3e8;
w=input('输个数');
k=w*1e7*ne/c;
Tc1=1+1j*tan(k*h0)/2;Tc2=1j*tan(k*h0)/2;Tc3=-1j*tan(k*h0)/2;Tc4=1-1j*tan(k*h0)/2;
Tc=[Tc1,Tc2;Tc3,Tc4];
Th1=[exp(1j*k*h1),0;0,exp(-1j*k*h1)];
Th2=[exp(1j*k*h2),0;0,exp(-1j*k*h2)];
T=Th1*Tc*Th2*Tc*Th1;
A=cos(k*D)-(tan(k*h0).^2).*cos(k*D)/4-tan(k*h0).*sin(k*D)+(tan(k*h0).^2).*cos(k*Delta*D)/4;%cosqD
B=real(acos(A));
C1=T(1,2);
C2=exp(1i*B)-T(1,1);
z1=0:0.01:h1;
W1=C1*exp(1i*k*z1)+C2*exp(-1i*k*z1);%最终结果
figure
plot(z1,abs(W1));
grid on