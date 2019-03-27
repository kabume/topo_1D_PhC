%model D,Dr.Li 's paper + 2 column
% This program could simulate the paper of Li and Zww
clc,clear
%paramaters
da=100e-9;db=100e-9;D=da+db;l0=da;
dx=0:1:200;%The length of l_0 to the left end of A
dx=dx*1e-9;
%ne=1;%temp
na=3.2;%slab A
nb=1;%slab B  
nc=2;
c=3e8;
Z0=377;
P=[1,1;1/Z0,-1/Z0];
w=0:0.01:6;
ka=w*1e15*na/c;
kb=w*1e15*nb/c;
kc=w*1e15*nc/c;
num1=length(w);
num2=length(dx);
r0=zeros(num2,num1);
Arg=zeros(num2,num1);
parfor nn=1:num1%omega
    for mm=1:num2%Delta
        if mm<db/D*num2
            m=[exp(1i*ka(nn)*da),exp(-1i*ka(nn)*da);ka(nn)*exp(1i*ka(nn)*da),-ka(nn)*exp(-1i*ka(nn)*da)];
            n=[1,1;kb(nn),-kb(nn)];%[exp(1i*kb(nn)*da),exp(-1i*kb(nn)*da);kb(nn)*exp(1i*kb(nn)*da),-kb(nn)*exp(-1i*kb(nn)*da)];
            o=[exp(1i*kb(nn)*(db-dx(mm))),exp(-1i*kb(nn)*(db-dx(mm)));kb(nn)*exp(1i*kb(nn)*(db-dx(mm))),-kb(nn)*exp(-1i*kb(nn)*(db-dx(mm)))];
            p=[1,1;kc(nn),-kc(nn)];
            q=[1+1i*tan(kc(nn)*l0)/2,1i*tan(kc(nn)*l0)/2;-1i*tan(kc(nn)*l0)/2,1-1i*tan(kc(nn)*l0)/2];
            r=[1,1;kc(nn),-kc(nn)];
            s=[1,1;kb(nn),-kb(nn)];
            t=[exp(1i*kb(nn)*dx(mm)),exp(-1i*kb(nn)*dx(mm));kb(nn)*exp(1i*kb(nn)*dx(mm)),-kb(nn)*exp(-1i*kb(nn)*dx(mm))];
            u=[1,1;ka(nn),-ka(nn)];
            %T=inv(u)*t*inv(s)*r*q*inv(p)*o*inv(n)*m;
            T=inv(u)*t*inv(s)*r*q*inv(p)*o*inv(n)*m;
%             A=((.125*1i*kb(nn)^2*kc(nn)-.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)*exp(2*1i*kb(nn)*(db-dx(mm)))+((.25*kb(nn)*1i*ka(nn)*kc(nn)+.125*1i*kb(nn)^2*kc(nn)+.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)+.5*ka(nn)*kb(nn)^2+.25*kb(nn)^3+.25*ka(nn)^2*kb(nn))*exp(2*1i*kb(nn)*db)+(.125*1i*kb(nn)^2*kc(nn)-.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)*exp(2*1i*kb(nn)*dx(mm))+(-.25*kb(nn)*1i*ka(nn)*kc(nn)+.125*1i*kb(nn)^2*kc(nn)+.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)+.5*ka(nn)*kb(nn)^2-.25*kb(nn)^3-.25*ka(nn)^2*kb(nn))*exp(1i*(da*ka(nn)-db*kb(nn)))/(ka(nn)*kb(nn)^2);
%             B=((-.25*kb(nn)*1i*ka(nn)*kc(nn)+.125*1i*kb(nn)^2*kc(nn)+.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)*exp(2*1i*kb(nn)*(db-dx(mm)))+((.125*1i*kb(nn)^2*kc(nn)-.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)+.25*kb(nn)^3-.25*ka(nn)^2*kb(nn))*exp(2*1i*kb(nn)*db)+(.25*kb(nn)*1i*ka(nn)*kc(nn)+.125*1i*kb(nn)^2*kc(nn)+.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)*exp(2*1i*kb(nn)*dx(mm))+(.125*1i*kb(nn)^2*kc(nn)-.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)-.25*kb(nn)^3+.25*ka(nn)^2*kb(nn))*exp(-1i*(da*ka(nn)+db*kb(nn)))/(ka(nn)*kb(nn)^2); 
%             C=((-.25*kb(nn)*1i*ka(nn)*kc(nn)-.125*1i*kb(nn)^2*kc(nn)-.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)*exp(2*1i*kb(nn)*(db-dx(mm)))+((-.125*1i*kb(nn)^2*kc(nn)+.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)-.25*kb(nn)^3+.25*ka(nn)^2*kb(nn))*exp(2*1i*kb(nn)*db)+(.25*kb(nn)*1i*ka(nn)*kc(nn)-.125*1i*kb(nn)^2*kc(nn)-.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)*exp(2*1i*kb(nn)*dx(mm))+(-.125*1i*kb(nn)^2*kc(nn)+.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)+.25*kb(nn)^3-.25*ka(nn)^2*kb(nn))*exp(1i*(da*ka(nn)-db*kb(nn)))/(ka(nn)*kb(nn)^2); 
%             E=((-.125*1i*kb(nn)^2*kc(nn)+.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)*exp(2*1i*kb(nn)*(db-dx(mm)))+((.25*kb(nn)*1i*ka(nn)*kc(nn)-.125*1i*kb(nn)^2*kc(nn)-.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)+.5*ka(nn)*kb(nn)^2-.25*kb(nn)^3-.25*ka(nn)^2*kb(nn))*exp(2*1i*kb(nn)*db)+(-.125*1i*kb(nn)^2*kc(nn)+.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)*exp(2*1i*kb(nn)*dx(mm))+(-.25*kb(nn)*1i*ka(nn)*kc(nn)-.125*1i*kb(nn)^2*kc(nn)-.125*1i*ka(nn)^2*kc(nn))*tan(kc(nn)*l0)+.5*ka(nn)*kb(nn)^2+.25*kb(nn)^3+.25*ka(nn)^2*kb(nn))*exp(-1i*(da*ka(nn)+db*kb(nn)))/(ka(nn)*kb(nn)^2);
%             T=[A,B;C,E];
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
         if Temp1<=1
            Temp1=(AAA-sqrt(Temp0))/2;
            r0(mm,nn)=abs((Temp1-T(1,1))/T(1,2));
            Arg(mm,nn)=phase((Temp1-T(1,1))/T(1,2));
         else
            Temp1=(AAA+sqrt(Temp0))/2;
            r0(mm,nn)=abs((Temp1-T(1,1))/T(1,2));
            Arg(mm,nn)=phase((Temp1-T(1,1))/T(1,2));
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
% figure
% plot(w,r0(round(50/D*num2),:))