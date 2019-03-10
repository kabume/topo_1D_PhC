```matlab
%测试数据
clc,clear
%paramaters
L=100;D=2*L;h0=L/2;
Delta=-0.8:0.001:0.8;
h=(1+Delta)*L;h1=L-h/2;h2=h;ne=2;c=3e8;
w=0:0.001:0.6;
k=w*1e7*ne/c;
tic
semiinfinite(ne,k,h0,w,Delta,h1,h2)
toc
%bloch(w,k,D,h0,Delta)
%onecell(k,h0,w,Delta,h1,h2)
```



```matlab
%随频率变化的透反射图
%This program could draw the figure 3. h==>l, D==>Lambda
L=100;D=2*L;h0=L/2;
Delta=-0.8;
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
    t(n)=abs(1./T(2,2));r(n)=abs(-T(2,1)./T(2,2));
end
figure
plot(w,r,'k');
hold on 
plot(w,t);

```

```Matlab
%随频率和Δ变化的反射率和辐角的二维图,单个cell
function onecell(k,h0,w,Delta,h1,h2)
    Tc1=1+1j*tan(k*h0)/2;Tc2=1j*tan(k*h0)/2;Tc3=-1j*tan(k*h0)/2;Tc4=1-1j*tan(k*h0)/2;
    num1=length(w);
    num2=length(Delta);
    r=zeros(num2,num1);Arg=zeros(num2,num1);
    for m=1:num2
        for n=1:num1
            Tc=[Tc1(n),Tc2(n);Tc3(n),Tc4(n)];
            Th1=[exp(1j*k(n)*h1(m)),0;0,exp(-1j*k(n)*h1(m))];
            Th2=[exp(1j*k(n)*h2(m)),0;0,exp(-1j*k(n)*h2(m))];
            T=Th1*Tc*Th2*Tc*Th1;
            r(m,n)=abs(-T(2,1)./T(2,2));
            Arg(m,n)=angle(-T(2,1)./T(2,2));
        end
    end
    figure
    imagesc(Delta,w,r);
    figure
    imagesc(Delta,w,Arg);
end
```

```matlab
%随频率和Δ变化的反射率和辐角的二维图,单个cell和半无限长
function semiinfinite(ne,k,h0,w,Delta,h1,h2)
    %transfer matrix of Tc
    Z0=1/ne;
    Tc1=1+1j*tan(k*h0)/2;Tc2=1j*tan(k*h0)/2;Tc3=-1j*tan(k*h0)/2;Tc4=1-1j*tan(k*h0)/2;
    P=[1,1;1/Z0,-1/Z0];
    num1=length(w);
    num2=length(Delta);
    r=zeros(num2,num1);Arg=zeros(num2,num1);
    r1=zeros(num2,num1);Arg1=zeros(num2,num1);
    parfor n=1:num1
        for m=1:num2
            Tc=[Tc1(n),Tc2(n);Tc3(n),Tc4(n)];
            Th1=[exp(1j*k(n)*h1(m)),0;0,exp(-1j*k(n)*h1(m))];
            Th2=[exp(1j*k(n)*h2(m)),0;0,exp(-1j*k(n)*h2(m))];
            T=Th1*Tc*Th2*Tc*Th1;
            M=P*T*P^-1;
            Z=(M(1,1)-M(2,2)-sqrt((M(1,1)-M(2,2))^2+4*M(2,1)*M(1,2)))/(2*M(2,1));
            if real(Z)>0
                 Z=(M(1,1)-M(2,2)-sqrt((M(1,1)-M(2,2))^2+4*M(2,1)*M(1,2)))/(2*M(2,1));
            else
                Z=(M(1,1)-M(2,2)+sqrt((M(1,1)-M(2,2))^2+4*M(2,1)*M(1,2)))/(2*M(2,1));
            end
            r1(m,n)=abs((Z-Z0)/(Z+Z0));
            Arg1(m,n)=angle((Z-Z0)/(Z+Z0));
            r(m,n)=abs(-T(2,1)./T(2,2));
            Arg(m,n)=angle(-T(2,1)./T(2,2));
        end
    end
    subplot(2,2,1)
    imagesc(w,Delta,r);
    title("|r| of one cell")
    view([-90 90]);
    box('on');
    axis('ij');
    xlabel("\omega");ylabel("\Delta");
    subplot(2,2,2)
    imagesc(w,Delta,Arg);
    title("Arg of one cell")
    view([-90 90]);
    box('on');
    axis('ij');
    xlabel("\omega");ylabel("\Delta");
    subplot(2,2,3)
    imagesc(w,Delta,r1);
    title("|r| of semi-infinite systems")
    view([-90 90]);
    box('on');
    axis('ij');
    xlabel("\omega");ylabel("\Delta");
    subplot(2,2,4)
    imagesc(w,Delta,Arg1);
    title("Arg of semi-infinite systems")
    view([-90 90]);
    box('on');
    axis('ij');
    xlabel("\omega");ylabel("\Delta");
%     figure
%     imagesc(w,Delta,Arg1);
%     title("Arg of semi-infinite systems")
%     
end
```

```matlab
%Calculate the band diagram 
function bloch(w,k,D,h0,Delta)
    A=cos(k*D)-(tan(k*h0).^2).*cos(k*D)/4-tan(k*h0).*sin(k*D)+(tan(k*h0).^2).*cos(k*Delta*D)/4;%cosqD
    B=real(acos(A));
    plot(B,w,'k');
    hold on
    plot(-B,w,'k');
    hold off
end
```

```matlab
%奇对称和偶对称的场图，参考CT2014PRX中的附录2，可通过场分布图判断zak phase
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
```

```matlab
%半无限长的有效epsilon和mu,Z就是有效阻抗，ne必须是复数才能画出论文中的图，图不太精确，不过有效参数的正负判断和原图中的一致，文中定义的ne也是有效折射率
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
    epsilon(n)=(2+2i)/Z;
    mu(n)=(2+2i)*Z;
end
figure
plot(w,epsilon,'k');
hold on 
plot(w,mu);
axis([0 0.6 -10 10])
grid on
```

```matlab
%AB中加个杆。化简结果和仿真表明，r和Delta无关，Delta只和Tc的形式以及B层中的k和l有关
% This program could simulate the paper of Li and Zww
clc,clear
%paramaters
L=100;D=2*L;h0=L/2;
Delta=0:001:D;%The length of l_0 to the left end of A
ne=1;%temp
ne1=5;%slab A
ne2=1.5;%slab B
h1=L;h2=L;
c=3e8;
w=0:0.001:0.6;
k=w*1e7*ne/c;
k1=w*1e7*ne1/c;
k2=w*1e7*ne2/c;
Z0=1/ne;% todo
Tc1=1+1j*tan(k*h0)/2;Tc2=1j*tan(k*h0)/2;Tc3=-1j*tan(k*h0)/2;Tc4=1-1j*tan(k*h0)/2;
P=[1,1;1/Z0,-1/Z0];
num1=length(w);
num2=length(Delta);
r=zeros(num2,num1);Arg=zeros(num2,num1);
r1=zeros(num2,num1);Arg1=zeros(num2,num1);
parfor n=1:num1%omega
    for m=1:num2%Delta
        Tc=[Tc1(n),Tc2(n);Tc3(n),Tc4(n)];
        if m<num2/2
            Th10=[exp(1j*k1(n)*Delta(m)),0;0,exp(-1j*k1(n)*Delta(m))];%todo:Th1 and Th2 will change when h0 change
            Th11=[exp(1j*k1(n)*(h1-Delta(m))),0;0,exp(-1j*k1(n)*(h1-Delta(m)))];
            Th2=[exp(1j*k2(n)*h2),0;0,exp(-1j*k2(n)*h2)];
            T=Th10*Tc*Th11*Th2;
        else
            Th1=[exp(1j*k1(n)*h1),0;0,exp(-1j*k1(n)*h1)];%todo:Th1 and Th2 will change when h0 change
            Th20=[exp(1j*k2(n)*(Delta(m)-h1)),0;0,exp(-1j*k2(n)*(Delta(m)-h1))];
            Th21=[exp(1j*k2(n)*(D-Delta(m))),0;0,exp(-1j*k2(n)*(D-Delta(m)))];
            T=Th1*Th20*Tc*Th21;
        end
        M=P*T*P^-1;
        Z=(M(1,1)-M(2,2)-sqrt((M(1,1)-M(2,2))^2+4*M(2,1)*M(1,2)))/(2*M(2,1));
        if real(Z)>0
            Z=(M(1,1)-M(2,2)-sqrt((M(1,1)-M(2,2))^2+4*M(2,1)*M(1,2)))/(2*M(2,1));
        else
            Z=(M(1,1)-M(2,2)+sqrt((M(1,1)-M(2,2))^2+4*M(2,1)*M(1,2)))/(2*M(2,1));
        end
        r1(m,n)=abs((Z-Z0)/(Z+Z0));
        Arg1(m,n)=angle((Z-Z0)/(Z+Z0));
%         if Arg1(m,n)<-pi+0.001&&Arg1(m,n)>-pi-0.001
%             Arg1(m,n)=pi;
%         end
        r(m,n)=abs(-T(2,1)./T(2,2));
        Arg(m,n)=angle(-T(2,1)./T(2,2));
    end
end
subplot(2,2,1)
imagesc(w,Delta,r);
title("|r| of one cell")
view([-90 90]);
box('on');
axis('ij');
xlabel("\omega");ylabel("\Delta");
subplot(2,2,2)
imagesc(w,Delta,Arg);
title("Arg of one cell")
view([-90 90]);
box('on');
axis('ij');
xlabel("\omega");ylabel("\Delta");
subplot(2,2,3)
imagesc(w,Delta,r1);
title("|r| of semi-infinite systems")
view([-90 90]);
box('on');
axis('ij');
xlabel("\omega");ylabel("\Delta");
subplot(2,2,4)
imagesc(w,Delta,Arg1);
title("Arg of semi-infinite systems")
view([-90 90]);
box('on');
axis('ij');
xlabel("\omega");ylabel("\Delta");
```



#### Symmetry

* $\Delta=-0.8$ 

第一个对称点(0,0)

第二个对称点(3.142,0.1602)

第三个对称点(3.142,0.2316)

第四个对称点(0,0.3322)

第五个对称点(0,0.3768)

第六个对称点(3.142,0.4052)

#### Questions

* **fig2**
  * 图中半无限长的相位图，我画出来会有很多点相位恰好相反，abs后与原图差不多，其他三个图和原图差不多。==答案：对于$cos$函数$-\pi=\pi​$。==
* **fig3**：
  * acos一般有正数也有负数，怎么确定正负呢？比如参考文献**39**中的公式**(3)**利用$cos(nkd)$求$n$。==对于band gap图，本身就是对称的，其他的根据文献来确定正负。==
  * 如何通过传输矩阵法画场分布图？以此来确认奇对称和偶对称。==通过传输矩阵往前推，要画出无限长体系,参考CT2014PRX中的附录2画出无限长体系的第一个单元的场分布图就可以判断对称性和zak phase==
  * 如何利用公式(5)和公式(6)算出Zak phase。