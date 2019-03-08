
%���ô��������㵥����ʵ�͸���ʣ���˶�������p136
function [LL]=DL(n1,n2,h1,h2)
a=h1+h2;
C2=[1;0];%�Ҷ˵�͸������뷴�����
%lambda=100:1:900;
omega=0:0.01e7:3e8;
%f=3e8./lambda;
%k1=2*pi*n1./(lambda*1e-9);
k1=omega/n1;
e1=k1*h1;
%k2=2*pi*n2./(lambda*1e-9);
k2=omega/n2;
e2=k2*h2;
c1=0.5*1i*(k1/k2+k2/k1);
c2=0.5*1i*(k1/k2-k2/k1);%��ʽ���ַֽ�
num=length(omega);
T=zeros(1,num);%ΪTԤ�����ڴ�
for j=1:num
    M1=[exp(1i*e1(j))*(cos(e2(j))+c1*sin(e2(j))),exp(-1i*(e1(j)))*(-c2*sin(e2(j)));exp(1i*e1(j))*c2*sin(e2(j)),exp(-1i*e1(j))*(cos(e2(j))-c1*sin(e2(j)))];%��������
    N=10;%������
    M=M1^N;%������ڹ��Ӿ������������
    C1=C2.*inv(M);%��˵�͸������ͷ������
    T(j)=abs(C2(1)/C1(1))^2;
    K(j)=acos(0.5*(M1(1,1)+M1(2,2)))/(h1+h2);
    LL(j)=cos(K(j)*a);
    y(j)=1;
    %T(j)=abs((M(1,1)*M(2,2)-M(2,1)*M(1,2))./((A2*M(2,2)-B2*M(1,2))))^2;%�������͸���ʱ��ʽ
end

%subplot(2,2,3)
figure
plot(omega/1e7,LL)
hold on
plot(omega/1e7,y)
%plot(K*a/2/pi,2*pi./lambda)
%G=diff(T);
%plot(G)
end
