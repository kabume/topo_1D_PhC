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
%画出AB结构+共振柱的反射谱
%提供阻抗法和布洛赫法算Arg和r
% This program could simulate the paper of Li and Zww
clc,clear
%paramaters
da=100;db=100;D=da+db;l0=da/2;
dx=0:0.1:D;%The length of l_0 to the left end of A
%ne=1;%temp
na=3.2;%slab A
nb=1;%slab B  
nc=2;
c=3e8;
% Z0=377;
% P=[1,1;1/Z0,-1/Z0];
w=0:0.001:1.6;
ka=w*1e7*na/c;
kb=w*1e7*nb/c;
kc=w*1e7*nc/c;
num1=length(w);
num2=length(dx);
r0=zeros(num2,num1);
Arg=zeros(num2,num1);
parfor nn=1:num1%omega
    for mm=1:num2%Delta
        if mm<db/D*num2
%             m=[exp(1i*ka(nn)*da),exp(-1i*ka(nn)*da);ka(nn)*exp(1i*ka(nn)*da),-ka(nn)*exp(-1i*ka(nn)*da)];
%             n=[1,1;kb(nn),-kb(nn)];%[exp(1i*kb(nn)*da),exp(-1i*kb(nn)*da);kb(nn)*exp(1i*kb(nn)*da),-kb(nn)*exp(-1i*kb(nn)*da)];
%             o=[exp(1i*kb(nn)*(db-dx(mm))),exp(-1i*kb(nn)*(db-dx(mm)));kb(nn)*exp(1i*kb(nn)*(db-dx(mm))),-kb(nn)*exp(-1i*kb(nn)*(db-dx(mm)))];
%             p=[1,1;kc(nn),-kc(nn)];
%             q=[1+1i*tan(kc(nn)*l0)/2,1i*tan(kc(nn)*l0)/2;-1i*tan(kc(nn)*l0)/2,1-1i*tan(kc(nn)*l0)/2];
%             r=[1,1;kc(nn),-kc(nn)];
%             s=[1,1;kb(nn),-kb(nn)];
%             t=[exp(1i*kb(nn)*dx(mm)),exp(-1i*kb(nn)*dx(mm));kb(nn)*exp(1i*kb(nn)*dx(mm)),-kb(nn)*exp(-1i*kb(nn)*dx(mm))];
%             u=[1,1;ka(nn),-ka(nn)];
%             T=inv(u)*t*inv(s)*r*q*inv(p)*o*inv(n)*m;
            A=((.125*1i*kb(nn)^2*kc(nn)-.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)*exp(2*1i*kb(nn)*(db-dx(mm)))+((.25*kb(nn)*1i*ka(nn)*kc(nn)+.125*1i*kb(nn)^2*kc(nn)+.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)+.5*ka(nn)*kb(nn)^2+.25*kb(nn)^3+.25*ka(nn)^2*kb(nn))*exp(2*1i*kb(nn)*db)+(.125*1i*kb(nn)^2*kc(nn)-.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)*exp(2*1i*kb(nn)*dx(mm))+(-.25*kb(nn)*1i*ka(nn)*kc(nn)+.125*1i*kb(nn)^2*kc(nn)+.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)+.5*ka(nn)*kb(nn)^2-.25*kb(nn)^3-.25*ka(nn)^2*kb(nn))*exp(1i*(da*ka(nn)-db*kb(nn)))/(ka(nn)*kb(nn)^2);
            B=((-.25*kb(nn)*1i*ka(nn)*kc(nn)+.125*1i*kb(nn)^2*kc(nn)+.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)*exp(2*1i*kb(nn)*(db-dx(mm)))+((.125*1i*kb(nn)^2*kc(nn)-.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)+.25*kb(nn)^3-.25*ka(nn)^2*kb(nn))*exp(2*1i*kb(nn)*db)+(.25*kb(nn)*1i*ka(nn)*kc(nn)+.125*1i*kb(nn)^2*kc(nn)+.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)*exp(2*1i*kb(nn)*dx(mm))+(.125*1i*kb(nn)^2*kc(nn)-.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)-.25*kb(nn)^3+.25*ka(nn)^2*kb(nn))*exp(-1i*(da*ka(nn)+db*kb(nn)))/(ka(nn)*kb(nn)^2); 
            C=((-.25*kb(nn)*1i*ka(nn)*kc(nn)-.125*1i*kb(nn)^2*kc(nn)-.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)*exp(2*1i*kb(nn)*(db-dx(mm)))+((-.125*1i*kb(nn)^2*kc(nn)+.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)-.25*kb(nn)^3+.25*ka(nn)^2*kb(nn))*exp(2*1i*kb(nn)*db)+(.25*kb(nn)*1i*ka(nn)*kc(nn)-.125*1i*kb(nn)^2*kc(nn)-.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)*exp(2*1i*kb(nn)*dx(mm))+(-.125*1i*kb(nn)^2*kc(nn)+.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)+.25*kb(nn)^3-.25*ka(nn)^2*kb(nn))*exp(1i*(da*ka(nn)-db*kb(nn)))/(ka(nn)*kb(nn)^2); 
            E=((-.125*1i*kb(nn)^2*kc(nn)+.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)*exp(2*1i*kb(nn)*(db-dx(mm)))+((.25*kb(nn)*1i*ka(nn)*kc(nn)-.125*1i*kb(nn)^2*kc(nn)-.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)+.5*ka(nn)*kb(nn)^2-.25*kb(nn)^3-.25*ka(nn)^2*kb(nn))*exp(2*1i*kb(nn)*db)+(-.125*1i*kb(nn)^2*kc(nn)+.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)*exp(2*1i*kb(nn)*dx(mm))+(-.25*kb(nn)*1i*ka(nn)*kc(nn)-.125*1i*kb(nn)^2*kc(nn)-.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)+.5*ka(nn)*kb(nn)^2+.25*kb(nn)^3+.25*ka(nn)^2*kb(nn))*exp(-1i*(da*ka(nn)+db*kb(nn)))/(ka(nn)*kb(nn)^2);
            T=[A,B;C,E];
        else    
            m=[exp(1i*ka(nn)*(da+db-dx(mm))),exp(-1i*ka(nn)*(da+db-dx(mm)));ka(nn)*exp(1i*ka(nn)*(da+db-dx(mm))),-ka(nn)*exp(-1i*ka(nn)*(da+db-dx(mm)))];
            n=[1,1;kc(nn),-kc(nn)];
            o=[1+1i*tan(kc(nn)*l0)/2,1i*tan(kc(nn)*l0)/2;-1i*tan(kc(nn)*l0)/2,1-1i*tan(kc(nn)*l0)/2];
            p=[1,1;kc(nn),-kc(nn)];
            q=[1,1;ka(nn),-ka(nn)];
            r=[exp(1i*ka(nn)*(dx(mm)-db)),exp(-1i*ka(nn)*(dx(mm)-db));ka(nn)*exp(1i*ka(nn)*(dx(mm)-db)),-ka(nn)*exp(-1i*ka(nn)*(dx(mm)-db))];
            s=[1,1;kb(nn),-kb(nn)];
            t=[exp(1i*kb(nn)*db),exp(-1i*kb(nn)*db);kb(nn)*exp(1i*kb(nn)*db),-kb(nn)*exp(-1i*kb(nn)*db)];
            u=[1,1;ka(nn),-ka(nn)];
            T=inv(u)*t*inv(s)*r*inv(q)*p*o*inv(n)*m;
        end
         AAA=T(1,1)+T(2,2);
         Temp0=AAA^2-4;
         Temp1=(AAA-sqrt(Temp0))/2;
         r0(mm,nn)=abs((Temp1-T(1,1))/T(1,2));
         if r0(mm,nn)<=1
            Temp1=(AAA-sqrt(Temp0))/2;
            r0(mm,nn)=abs((Temp1-T(1,1))/T(1,2));
            Arg(mm,nn)=angle((Temp1-T(1,1))/T(1,2));
         else
            Temp1=(AAA+sqrt(Temp0))/2;
            r0(mm,nn)=abs((Temp1-T(1,1))/T(1,2));
            Arg(mm,nn)=angle((Temp1-T(1,1))/T(1,2));
        end
%        M=P*T*P^-1;
%        Z=(M(1,1)-M(2,2)-sqrt((M(1,1)-M(2,2))^2+4*M(2,1)*M(1,2)))/(2*M(2,1));
%        if real(Z)>=0
%             Z=(M(1,1)-M(2,2)-sqrt((M(1,1)-M(2,2))^2+4*M(2,1)*M(1,2)))/(2*M(2,1));
%        else
%            Z=(M(1,1)-M(2,2)+sqrt((M(1,1)-M(2,2))^2+4*M(2,1)*M(1,2)))/(2*M(2,1));
%        end
%         Temp(mm,nn)=(Z-Z0)/(Z+Z0);
%         r0(mm,nn)=abs(Temp(mm,nn));
%         Arg(mm,nn)=angle(Temp(mm,nn));
    end
end
figure
imagesc(w,dx,r0);
title("|r| of semiinfinite")
view([-90 90]);
box('on');
axis('ij');
xlabel("\omega");ylabel("\Delta");
figure
imagesc(w,dx,Arg);
title("Arg of semiinfinite")
view([-90 90]);
box('on');
axis('ij');
xlabel("\omega");ylabel("\Delta");
figure
plot(w,r0(round(50/D*num2),:))%奇异点处的反射谱[0.275,50]
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

#### Debug

$T_{li}=\left(
\begin{array}{cccc}
e^{jkl_i} &  0\\
0 & e^{-jkl_i}  \\
\end{array}
\right) $,$i=1,2,3$分别代表A，B，C

$T_{li}$的斜对角等于0，且$e^{jkl_{B_1}}e^{jkl_{B_2}}=e^{jkl_{B}}$,对于传输矩阵$T_AT_BT_C$来说，当把B层分为两段但是总长度不变时，$Z_{T_AT_BT_C}=Z_{T_{B_1}T_AT_{B_2}T_C}$。