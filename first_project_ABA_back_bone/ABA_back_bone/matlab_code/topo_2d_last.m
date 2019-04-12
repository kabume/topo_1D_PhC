%model D,Dr.Li 's paper + 2 column
% This program could simulate the paper of Li and Zww
clc,clear
%paramaters
da=100e-9;db=100e-9;D=da+db;l0=da;
Dt=D*1e9;
dx0=0:0.1:Dt;
dx=dx0*1e-9;
na=3.2;
nb=1;
nc=2;
c=3e8;
w=0.1:0.01:15;
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
            T=inv(u)*t*inv(s)*r*q*inv(p)*o*inv(n)*m;
            
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
            if r0(mm,nn)<=1.0001&&r0(mm,nn)>=0.9999
                Temp1=(AAA-sqrt(Temp0))/2;
                r0(mm,nn)=abs((Temp1-T(1,1))/T(1,2));
                Arg(mm,nn)=angle((Temp1-T(1,1))/T(1,2));
            else
                Temp1=(AAA+sqrt(Temp0))/2;
                r0(mm,nn)=abs((Temp1-T(1,1))/T(1,2));
                Arg(mm,nn)=angle((Temp1-T(1,1))/T(1,2));
            end
        end
     end
end
figure
imagesc(w,dx0,r0);
ylim([0,D*1e9])
title("|r| of semiinfinite")
view([-90 90]);
box('on');
axis('ij');
xlabel("\omega");ylabel("\Delta");
figure
imagesc(w,dx0,Arg);
ylim([0,D*1e9])
view([-90 90]);
box('on');
axis('ij');
title("Arg of semiinfinite")
xlabel("\omega");ylabel("\Delta");