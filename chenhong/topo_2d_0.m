%两个柱子的结构
% This program could simulate the paper of Li and Zww
clc,clear
%paramaters
L=100;D=2*L;h0=L/2;
Delta=0:001:D;%The length of l_0 to the left end of A
delta=30;
ne=2;%temp
ne1=3;%slab A
ne2=1;%slab B
h1=50;h2=150;
c=3e8;
w=0:0.001:0.6;
k=w*1e7*ne/c;
k1=w*1e7*ne1/c;
k2=w*1e7*ne2/c;
Z0=1;% todo
Tc1=1+1j*tan(k*h0)/2;Tc2=1j*tan(k*h0)/2;Tc3=-1j*tan(k*h0)/2;Tc4=1-1j*tan(k*h0)/2;
P=[1,1;1/Z0,-1/Z0];
num1=length(w);
num2=length(Delta);
r=zeros(num2,num1);Arg=zeros(num2,num1);
r1=zeros(num2,num1);Arg1=zeros(num2,num1);
parfor n=1:num1%omega
    for m=1:num2%Delta
        Tc=[Tc1(n),Tc2(n);Tc3(n),Tc4(n)];
        if m<(h1-delta)/D*num2
            T_A1=[exp(1j*k1(n)*Delta(m)),0;0,exp(-1j*k1(n)*Delta(m))];
            T_A2=[exp(1j*k1(n)*delta),0;0,exp(-1j*k1(n)*delta)];
            T_A3=[exp(1j*k1(n)*(h1-Delta(m)-delta)),0;0,exp(-1j*k1(n)*(h1-Delta(m)-delta))];
            T_B1=[exp(1j*k2(n)*h2),0;0,exp(-1j*k2(n)*h2)];
            T=T_A1*Tc*T_A2*Tc*T_A3*T_B1;
        elseif m>=(h1-delta)/D*num2 && m<h1/D*num2
            T_A1=[exp(1j*k1(n)*Delta(m)),0;0,exp(-1j*k1(n)*Delta(m))];
            T_A2=[exp(1j*k1(n)*(h1-Delta(m))),0;0,exp(-1j*k1(n)*(h1-Delta(m)))];
            T_B1=[exp(1j*k2(n)*(delta-h1+Delta(m))),0;0,exp(-1j*k2(n)*(delta-h1+Delta(m)))];
            T_B2=[exp(1j*k2(n)*(D-delta-Delta(m))),0;0,exp(-1j*k2(n)*(D-delta-Delta(m)))];
            T=T_A1*Tc*T_A2*T_B1*Tc*T_B2;
        elseif m>=h1/D*num2 && m<(D-delta)/D*num2
            T_A1=[exp(1j*k1(n)*h1),0;0,exp(-1j*k1(n)*h1)];
            T_B1=[exp(1j*k2(n)*(Delta(m)-h1)),0;0,exp(-1j*k2(n)*(Delta(m)-h1))];
            T_B2=[exp(1j*k2(n)*delta),0;0,exp(-1j*k2(n)*delta)];
            T_B3=[exp(1j*k2(n)*(D-delta-Delta(m))),0;0,exp(-1j*k2(n)*(D-delta-Delta(m)))];
            T=T_A1*T_B1*Tc*T_B2*Tc*T_B3;
        else
            T_A1=[exp(1j*k1(n)*(delta-D+Delta(m))),0;0,exp(-1j*k1(n)*(delta-D+Delta(m)))];
            T_A2=[exp(1j*k1(n)*(h1-delta+D-Delta(m))),0;0,exp(-1j*k1(n)*(h1-delta+D-Delta(m)))];
            T_B1=[exp(1j*k2(n)*(h2-D+Delta(m))),0;0,exp(-1j*k2(n)*(h2-D+Delta(m)))];
            T_B2=[exp(1j*k2(n)*(D-Delta(m))),0;0,exp(-1j*k2(n)*(D-Delta(m)))];
            %T=T_A1*Tc*T_A2*T_B1*Tc*T_B2;
            T=T_A2*T_B1*Tc*T_B2*T_A1*Tc;
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