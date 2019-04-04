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
w=0.1:0.01:5;
ka=w*1e15*na/c;
kb=w*1e15*nb/c;
kc=w*1e15*nc/c;
num1=length(w);
num2=length(dx);
r0=zeros(num2,num1);
Arg=zeros(num2,num1);
for nn=1:num1%omega
    for mm=1:num2%Delta
        if mm<da/2/D*num2
            m=[exp(1i*ka(nn)*da),exp(-1i*ka(nn)*da);ka(nn)*exp(1i*ka(nn)*da),-ka(nn)*exp(-1i*ka(nn)*da)];
            n=[1,1;kb(nn),-kb(nn)];%[exp(1i*kb(nn)*da),exp(-1i*kb(nn)*da);kb(nn)*exp(1i*kb(nn)*da),-kb(nn)*exp(-1i*kb(nn)*da)];
            o=[exp(1i*kb(nn)*(db-dx(mm))),exp(-1i*kb(nn)*(db-dx(mm)));kb(nn)*exp(1i*kb(nn)*(db-dx(mm))),-kb(nn)*exp(-1i*kb(nn)*(db-dx(mm)))];
            p=[1,1;kc(nn),-kc(nn)];
            q=[1+1i*tan(kc(nn)*l0)/2,1i*tan(kc(nn)*l0)/2;-1i*tan(kc(nn)*l0)/2,1-1i*tan(kc(nn)*l0)/2];
            r=[1,1;kc(nn),-kc(nn)];
            s=[1,1;kb(nn),-kb(nn)];
            t=[exp(1i*kb(nn)*dx(mm)),exp(-1i*kb(nn)*dx(mm));kb(nn)*exp(1i*kb(nn)*dx(mm)),-kb(nn)*exp(-1i*kb(nn)*dx(mm))];
            u=[1,1;ka(nn),-ka(nn)];
            %A(mm,nn)=-exp(-1i*ka(nn)*da)*((exp(-1i*kb(nn)*dx(mm))*(ka(nn)-kb(nn))*(1i*kc(nn)*tan(kc(nn)*l0)-2*kb(nn))-1i*kc(nn)*(ka(nn)+kb(nn))*tan(kc(nn)*l0)*exp(1i*kb(nn)*dx(mm)))*(ka(nn)+kb(nn))*exp(-1i*kb(nn)*(-dx(mm)+db))-(1i*kc(nn)*(ka(nn)-kb(nn))*tan(kc(nn)*l0)*exp(-1i*kb(nn)*dx(mm))-(1i*kc(nn)*tan(kc(nn)*l0)+2*kb(nn))*(ka(nn)+kb(nn))*exp(1i*kb(nn)*dx(mm)))*exp(1i*kb(nn)*(-dx(mm)+db))*(ka(nn)-kb(nn)))/(8*ka(nn)*kb(nn)^2);
            %T=inv(u)*t*inv(s)*r*q*inv(p)*o*inv(n)*m;
            T=[(-i * kc(nn) * (ka(nn) - kb(nn)) * tan(kc(nn) * l0) * (ka(nn) + kb(nn)) * exp(-i * kb(nn) * db + 2*i * ka(nn) * dx(mm)) - (ka(nn) - kb(nn)) ^ 2 * (i * kc(nn) * tan(kc(nn) * l0) + 2 * ka(nn)) * exp(i * ka(nn) * da + -i * kb(nn) * db) + (i * kc(nn) * (ka(nn) - kb(nn)) * tan(kc(nn) * l0) * exp(i * kb(nn) * db + 2*i * ka(nn) * dx(mm)) + (i * kc(nn) * tan(kc(nn) * l0) + 2 * ka(nn)) * (ka(nn) + kb(nn)) * exp(i * ka(nn) * da + i * kb(nn) * db)) * (ka(nn) + kb(nn))) / ka(nn) ^ 2 / kb(nn) / 8 (-i * sin(kc(nn) * l0) * kc(nn) * (ka(nn) - kb(nn)) ^ 2 * exp(i * ka(nn) * (da - 2 * dx(mm)) + -i * kb(nn) * db) - 4 * (ka(nn) + kb(nn)) * (-0.1e1 / 0.4e1*i * sin(kc(nn) * l0) * kc(nn) * (ka(nn) + kb(nn)) * exp(i * ka(nn) * (da - 2 * dx(mm)) + i * kb(nn) * db) + (i * cos(kc(nn) * l0) * ka(nn) + (kc(nn) * sin(kc(nn) * l0) / 0.2e1)) * (ka(nn) - kb(nn)) * sin(kb(nn) * db))) / ka(nn) ^ 2 / kb(nn) / cos(kc(nn) * l0) / 8; (-i * sin(kc(nn) * l0) * kc(nn) * (ka(nn) + kb(nn)) ^ 2 * exp(-i * (da - 2 * dx(mm)) * ka(nn) + -i * kb(nn) * db) + 4 * (ka(nn) - kb(nn)) * (0.1e1 / 0.4e1*i * sin(kc(nn) * l0) * kc(nn) * (ka(nn) - kb(nn)) * exp(-i * (da - 2 * dx(mm)) * ka(nn) + i * kb(nn) * db) + (i * cos(kc(nn) * l0) * ka(nn) - (kc(nn) * sin(kc(nn) * l0) / 0.2e1)) * sin(kb(nn) * db) * (ka(nn) + kb(nn)))) / ka(nn) ^ 2 / kb(nn) / cos(kc(nn) * l0) / 8 (-(i * kc(nn) * tan(kc(nn) * l0) - 2 * ka(nn)) * (ka(nn) + kb(nn)) ^ 2 * exp(-i * ka(nn) * da + -i * kb(nn) * db) - (i * kc(nn) * tan(kc(nn) * l0) * (ka(nn) + kb(nn)) * exp(-i * kb(nn) * db + -2*i * ka(nn) * dx(mm)) - (i * kc(nn) * tan(kc(nn) * l0) - 2 * ka(nn)) * (ka(nn) - kb(nn)) * exp(-i * ka(nn) * da + i * kb(nn) * db) + -i * kc(nn) * exp(i * kb(nn) * db + -2*i * ka(nn) * dx(mm)) * tan(kc(nn) * l0) * (ka(nn) + kb(nn))) * (ka(nn) - kb(nn))) / ka(nn) ^ 2 / kb(nn) / 8;];

            T12(mm,nn)=T(1,2);
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
%             A(mm,nn)=(exp(i*ka(nn)*(dx(mm)-db))*cos(kb(nn)*db)/2+i*sin(kb(nn)*db)*(kb(nn)/ka(nn)*cos(ka(nn)*(dx(mm)-db))+i*ka(nn)/kb(nn)*sin(ka(nn)*(dx(mm)-db)))/2)*exp(-i*ka(nn)*(da+db-dx(mm)))+(exp(i*ka(nn)*(dx(mm)-db))*cos(kb(nn)*db)/2/ka(nn)+i*sin(kb(nn)*db)*(ka(nn)/kb(nn)*cos(ka(nn)*(dx(mm)-db))+i*kb(nn)/ka(nn)*sin(ka(nn)*(dx(mm)-db)))/2/ka(nn))*exp(-i*ka(nn)*(da+db-dx(mm)))*(i*kc(nn)*tan(kc(nn)*l0)-ka(nn));
            
            T=inv(u)*t*inv(s)*r*inv(q)*p*o*inv(n)*m;
%             T12(mm,nn)=T(1,2);
%             if abs(A(mm,nn))==abs(T12(mm,nn))
%                 F(mm,nn)=1;
%             else
%                 F(mm,nn)=0;
%             end
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
k1=((pi./kb+db)/2)*1e9;
k1=shaixuank(k1);
k_1=((-pi./kb+db)/2)*1e9;
k_1=shaixuank(k_1);
k0=ones(num1,1)*50;

p0=-(-db.*ka+atan(ka.*kc.*tan(kc.*l0).*tan(kb.*db)./(ka.^2.*tan(kb.*db)-kb.^2.*tan(kb.*db)+kb.*kc.*tan(kc.*l0))))./ka*1e9;
p0=shaixuanp(p0);
p1=(-db.*ka+atan(ka.*kc.*tan(kc.*l0).*tan(kb.*db)./(ka.^2.*tan(kb.*db)-kb.^2.*tan(kb.*db)+kb.*kc.*tan(kc.*l0)))+pi)./ka*1e9;
p1=shaixuanp(p1);
p_1=-(-db.*ka+atan(ka.*kc.*tan(kc.*l0).*tan(kb.*db)./(ka.^2.*tan(kb.*db)-kb.^2.*tan(kb.*db)+kb.*kc.*tan(kc.*l0)))-pi)./ka*1e9;
p_1=shaixuanp(p_1);
p_2=-(-db.*ka+atan(ka.*kc.*tan(kc.*l0).*tan(kb.*db)./(ka.^2.*tan(kb.*db)-kb.^2.*tan(kb.*db)+kb.*kc.*tan(kc.*l0)))-2*pi)./ka*1e9;
p_2=shaixuanp(p_2);
p_3=-(-db.*ka+atan(ka.*kc.*tan(kc.*l0).*tan(kb.*db)./(ka.^2.*tan(kb.*db)-kb.^2.*tan(kb.*db)+kb.*kc.*tan(kc.*l0)))-3*pi)./ka*1e9;
p_3=shaixuanp(p_3);
figure
imagesc(w,dx0,r0);
ylim([0,D*1e9])
hold on
%plot(w,p0,'k',w,p1,'k',w,p_1,'k',w,p_2,'k',w,p_3,'k')
%plot(w,k0,'k',w,k1,'k',w,k_1,'k')
title("|r| of semiinfinite")
view([-90 90]);
box('on');
axis('ij');
xlabel("\omega");ylabel("\Delta");
hold off
figure
imagesc(w,dx0,Arg);
hold on
%plot(w,p0,'k',w,p1,'k',w,p_1,'k',w,p_2,'k',w,p_3,'k')
%plot(w,k0,'k',w,k1,'k',w,k_1,'k')
ylim([0,D*1e9])
view([-90 90]);
box('on');
axis('ij');
title("Arg of semiinfinite")
xlabel("\omega");ylabel("\Delta");

function k=shaixuank(k1)
    for ii=1:length(k1)
        if k1(ii)<=100&&k1(ii)>=0
            k(ii)=k1(ii);
        else
            k(ii)=NaN;
        end
    end
end
function p=shaixuanp(p1)
     for ii=1:length(p1)
        if p1(ii)<=200&&p1(ii)>=100
            p(ii)=p1(ii);
        else
            p(ii)=NaN;
        end
    end
end
