%利用传输矩阵计算单层介质的透射率，郭硕鸿第三版p136


[LL] = DL(2,5,100e-9,100e-9)
function [LL]=DL(n1,n2,h1,h2)
a=h1+h2;
c = 3e8;
C2=[1;0];%右端的透射振幅与反射振幅
%lambda=100:1:900;
omega=0:0.001e15:4e15;
%f=3e8./lambda;
%k1=2*pi*n1./(lambda*1e-9);
k1=omega*n1/c;
e1=k1*h1;
%k2=2*pi*n2./(lambda*1e-9);
k2=omega*n2/c;
e2=k2*h2;
c1=0.5*1i*(k1/k2+k2/k1);
c2=0.5*1i*(k1/k2-k2/k1);%公式部分分解
num=length(omega);
T=zeros(1,num);%为T预分配内存
for j=1:num
    M1=[exp(1i*e1(j))*(cos(e2(j))+c1*sin(e2(j))),exp(-1i*(e1(j)))*(-c2*sin(e2(j)));exp(1i*e1(j))*c2*sin(e2(j)),exp(-1i*e1(j))*(cos(e2(j))-c1*sin(e2(j)))];%特征矩阵
    N=10;%周期数
    M=M1^N;%多个周期光子晶体的特征矩阵
    C1=C2.*inv(M);%左端的透射振幅和反射振幅
    T(j)=abs(C2(1)/C1(1));
    R(j) = abs(-M(2,1)./M(2,2));
    K(j)=acos(0.5*(M1(1,1)+M1(2,2)))/(h1+h2);
    LL(j)=cos(K(j)*a);
    y(j)=1;
    %T(j)=abs((M(1,1)*M(2,2)-M(2,1)*M(1,2))./((A2*M(2,2)-B2*M(1,2))))^2;%计算出的透射率表达式
end

%subplot(2,2,3)
figure
plot(omega/1e7,T)
hold on
plot(omega/1e7,R)
%plot(K*a/2/pi,2*pi./lambda)
%G=diff(T);
%plot(G)
end