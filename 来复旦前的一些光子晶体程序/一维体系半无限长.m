%一维体系半无限长
clc,clear
n1=2.35;
n2=1.38;
h1=63.8e-9;
h2=108.6449e-9;
C2=[1;0];%右端的透射振幅与反射振幅
lambda=100:1:900;
k1=2*pi*n1./(lambda*1e-9);
e1=k1*h1;
k2=2*pi*n2./(lambda*1e-9);
e2=k2*h2;
c1=0.5*1i*(k1/k2+k2/k1);
c2=0.5*1i*(k1/k2-k2/k1);
num=length(lambda);
R=zeros(1,num);%为T预分配内存
for j=1:num
    M1=[exp(1i*e1(j))*(cos(e2(j))+c1*sin(e2(j))),exp(-1i*(e1(j)))*(-c2*sin(e2(j)));exp(1i*e1(j))*c2*sin(e2(j)),exp(-1i*e1(j))*(cos(e2(j))-c1*sin(e2(j)))];
    A=(M1(1,1)+M1(2,2)+sqrt((M1(1,1)+M1(2,2))^2-4))/2;
    R(j)=abs((A-M1(1,1))/M1(1,2));
    if R(j)>1
        A=(M1(1,1)+M1(2,2)-sqrt((M1(1,1)+M1(2,2))^2-4))/2;
        R(j)=abs((A-M1(1,1))/M1(1,2));
    else
        A=(M1(1,1)+M1(2,2)+sqrt((M1(1,1)+M1(2,2))^2-4))/2;
        R(j)=abs((A-M1(1,1))/M1(1,2));
    end
    B(j)=real(acos((M1(1,1)+M1(2,2))/2));
%     N=10;
%     M=M1^N;
%     C1=C2.*inv(M);%左端的透射振幅和反射振幅
%     T(j)=abs(C2(1)/C1(1))^2;
end
plot(lambda,R,'k')
hold on
plot(lambda,B)