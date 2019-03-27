clear
close all
clc;

%����
epsilon_0=8.854187817*10^(-12);
c=299792458;

%��ȡÿ��Ŀ�Ⱥͽ�糣��
data=xlsread('data.xlsx');
m=length(data(1,:));       % m:�ܲ���
d=data(2,:)*10^(-9);       % d:ÿ����
epsilon=data(1,:);         % epsilon:ÿ�����Խ�糣��
mur=ones(1,m);             % mur:��Դŵ���

%ָ��Ƶ��
omega=2*pi*(10^14)*input('������ָ��Ƶ�ʣ�*10^14��');

%����ÿһ���������
refractivity=zeros(1,m);
for m1=1:m
  refractivity=(epsilon.*mur).^0.5;
end
%����ÿһ���k
k=zeros(1,m+1);
k(1,1:m)=refractivity*omega/c;  k(1,m+1)=omega/c;

ElectronicField=zeros(2,m);  ElectronicField1=zeros(2,m);
W=[1;0];

for m1=1:m
  current_number=m-m1+1;
  d_temp=sum(d(1,1:current_number));    %d_temp:��ǰ��������Ӧ��λ������
  M_temp1=[exp(1i*k(1,current_number)*d(1,current_number)),exp(-1i*k(1,current_number)*d(1,current_number));k(1,current_number)*exp(1i*k(1,current_number)*d(1,current_number)),-k(1,current_number)*exp(-1i*k(1,current_number)*d(1,current_number))];
  M_temp2=[1,1;k(1,current_number+1),-k(1,current_number+1)];
  ElectronicField(:,current_number)=(M_temp1\M_temp2)*W;
%  TI=[exp(-1i*k(1,current_number)*d(1,current_number)),0;0,exp(1i*k(1,current_number)*d(1,current_number))];
  W=ElectronicField(:,current_number);
%  ElectronicField1(:,current_number)=[exp(1i*k(1,current_number)*d(1,current_number)),0;0,exp(-1i*k(1,current_number)*d(1,current_number))]*ElectronicField(:,current_number);
end

figure(1)
hold on
XPointNumber=100;
d_modify=zeros(1,m+1);
x=zeros(1,m*XPointNumber-m+1);
ElectronicMagnitude=zeros(1,m*XPointNumber-m+1);
ElectronicMagnitude_imag=zeros(1,m*XPointNumber-m+1);
ElectronicMagnitude_real=zeros(1,m*XPointNumber-m+1);
y=zeros(1,m*XPointNumber-m+1);
for m1=2:m+1
  d_modify(1,m1)=sum(d(1,1:(m1-1)));
end

for m1=1:m
  x_temp=linspace(d_modify(1,m1),d_modify(1,m1+1),XPointNumber);
  x(1,XPointNumber*(m1-1)+1-m1+1:XPointNumber*m1-m1+1)=x_temp;
  delta_x=d(1,m1)/(XPointNumber-1);
  for m2=1:XPointNumber
%    current_number=XPointNumber-m2+1;
%    x_current=delta_x*(current_number-1);
   %  if m1==1
       d_temp=delta_x*(m2-1);
    % else
     % d_temp=delta_x*(m2-1);%+sum(d(1,1:m1-1));
     %end
     TI=[exp(1i*k(1,m1)*d_temp),0;0,exp(-1i*k(1,m1)*d_temp)];
     W=TI*ElectronicField(:,m1);
     ElectronicMagnitude_real(1,XPointNumber*(m1-1)+m2-m1+1)=(real(W(1,1)+W(2,1)));
     ElectronicMagnitude_imag(1,XPointNumber*(m1-1)+m2-m1+1)=(imag(W(1,1)+W(2,1)));
     ElectronicMagnitude(1,XPointNumber*(m1-1)+m2-m1+1)=(abs(W(1,1)+W(2,1)));
     y(1,XPointNumber*(m1-1)+m2-m1+1)=epsilon(1,m1);
  end
% plot([(x_temp(1,1)-sum(d)/2)*10^9,(x_temp(1,end)-sum(d)/2)*10^9],[ElectronicMagnitude(1,XPointNumber*(m1-1)+1-m1+1),ElectronicMagnitude(1,XPointNumber*m1-m1+1)],'red*')
end
figure(1);
plot(x*10^9,ElectronicMagnitude_real,'-r',x*10^9,ElectronicMagnitude_imag,'-b','Linewidth',1);
hold on
plot(x*10^9,y*10,'-black');
legend('ʵ�����','�鲿���','epsilon*10')
xlabel('x��������/nm');  ylabel('���');
grid on
figure(2);
plot(x*10^9,ElectronicMagnitude,'-','Linewidth',1);
hold on
plot(x*10^9,y*500,'-black');
legend('��ֵ','epsilon*500')
xlabel('x��������/nm');  ylabel('��ֵ');