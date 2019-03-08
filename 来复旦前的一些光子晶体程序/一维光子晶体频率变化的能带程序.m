clc,clear
n1=2;n2=1;d1=40;d2=60;N=10;c=3e8;
C2=[1;0];%右端的透射振幅与反射振幅
omega=0:0.01:80;
omega=40;
%z1=d1;z2=d1;
k1=omega*1e6*n1/c;
e1=k1*d1;
k2=omega*1e6*n2/c;
e2=k2*d2;
c1=0.5*1i*(k1/k2+k2/k1);
c2=0.5*1i*(k1/k2-k2/k1);
num=length(omega);
M11=exp(1i*e1)*(cos(e2)+c1*sin(e2));M12=exp(-1i*e1)*(-c2*sin(e2));
M21=exp(1i*e1)*c2*sin(e2);M22=exp(-1i*e1)*(cos(e2)-c1*sin(e2));
A0=acos(cos(e1).*cos(e2)-c1.*sin(e1).*sin(e2)/1i);%Ka
%A=exp(1i*k2*d1);B=exp(-1i*k2*d1);C=exp(1i*k2*d1);D=-exp(-1i*k2*d1);
%E=exp(1i*k1*d1);F=exp(-1i*k1*d1);G=n1/n2*exp(1i*k1*d1);H=-n1/n2*exp(-1i*k1*d1);
%I=M12;J=exp(1i*A0)-M11;
%S11=((E*I+F*J)*D-(G*I+H*J)*B)/(A*D-B*C);S12=((E*I+F*J)*C-(G*I+H*J)*A)/(B*C-A*D);
%T=zeros(1,num);%为T预分配内存
% for j=1:num
%     M1=[exp(1i*e1(j))*(cos(e2(j))+c1*sin(e2(j))),exp(-1i*(e1(j)))*(-c2*sin(e2(j)));exp(1i*e1(j))*c2*sin(e2(j)),exp(-1i*e1(j))*(cos(e2(j))-c1*sin(e2(j)))];
%     M=M1^N;
%     C1=C2.*inv(M);%左端的透射振幅和反射振幅
%     T(j)=abs(C2(1)/C1(1))^2;
% end

%figure
%plot(omega,T,'k')
%E1=M12*exp(1i*k1*(z1+d1/2))+(exp(1i*A0)-M11)*exp(-1i*k1*(z1+d1/2));%Slab 1电场强度
%E2=S11*exp(1i*k2*(z2+d1/2))+S12*exp(-1i*k2*(z2+d1/2));%Slab 2电场强度
plot(A0,omega,'k')
hold on
plot(-A0,omega,'k')
