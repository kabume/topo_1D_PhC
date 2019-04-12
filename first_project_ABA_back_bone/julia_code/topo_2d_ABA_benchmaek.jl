#并行版本
#ABA结构
using SharedArrays
using Distributed
using BenchmarkTools
addprocs(36)#除了第一次运行，下次运行的时候要注释掉
function testbm()
    da=100e-9;db=100e-9;D=da+db;l0=da/3.2;#单位m
    Dt=D*1e9;
    dx0=0:0.1:Dt;
    dx=dx0*1e-9
    na=3.2;#slab A
    nb=1;#slab B
    nc=2;
    c=3e8;
    w=0.1:0.01:15;
    ka=w*1e15*na/c;
    kb=w*1e15*nb/c;
    kc=w*1e15*nc/c;
    num1=length(w);
    num2=length(dx);
    r0 = zeros(num2,num1);
    Arg = zeros(num2,num1);
    r0 = SharedArray{Float64}(num2, num1);
    Arg = SharedArray{Float64}(num2, num1);
    rr = CartesianIndices(size(r0));
    @sync begin#同步
        @distributed for k in rr
            mm,nn=k.I;
            if mm<da/2D*num2
                m=[exp(im*ka[nn]*dx[mm]) exp(-im*ka[nn]*dx[mm]);ka[nn]*exp(im*ka[nn]*dx[mm]) -ka[nn]*exp(-im*ka[nn]*dx[mm])];
                n=[1 1;kc[nn] -kc[nn]];
                o=[1+im*tan(kc[nn]*l0)/2 im*tan(kc[nn]*l0)/2;-im*tan(kc[nn]*l0)/2 1-im*tan(kc[nn]*l0)/2];
                p=[1 1;kc[nn] -kc[nn]];
                q=[1 1;ka[nn] -ka[nn]];
                r=[exp(im*ka[nn]*(da/2-dx[mm])) exp(-im*ka[nn]*(da/2-dx[mm]));ka[nn]*exp(im*ka[nn]*(da/2-dx[mm])) -ka[nn]*exp(-im*ka[nn]*(da/2-dx[mm]))];
                s=[1 1;kb[nn] -kb[nn]];
                t=[exp(im*kb[nn]*db) exp(-im*kb[nn]*db);kb[nn]*exp(im*kb[nn]*db) -kb[nn]*exp(-im*kb[nn]*db)];
                u=[1 1;ka[nn] -ka[nn]];
                v=[exp(im*ka[nn]*da/2) 0;0 exp(-im*ka[nn]*da/2)];
                T=v*inv(u)*t*inv(s)*r*inv(q)*p*o*inv(n)*m;
            elseif mm>=da/2D*num2 && mm<(da/2+db)/D*num2
                m=[exp(im*ka[nn]*da/2) exp(-im*ka[nn]*da/2);ka[nn]*exp(im*ka[nn]*da/2) -ka[nn]*exp(-im*ka[nn]*da/2)];
                n=[1 1;kb[nn] -kb[nn]];
                o=[exp(im*kb[nn]*(dx[mm]-da/2)) exp(-im*kb[nn]*(dx[mm]-da/2));kb[nn]*exp(im*kb[nn]*(dx[mm]-da/2)) -kb[nn]*exp(-im*kb[nn]*(dx[mm]-da/2))];
                p=[1 1;kc[nn] -kc[nn]];
                q=[1+im*tan(kc[nn]*l0)/2 im*tan(kc[nn]*l0)/2;-im*tan(kc[nn]*l0)/2 1-im*tan(kc[nn]*l0)/2];
                r=[1 1;kc[nn] -kc[nn]];
                s=[1 1;kb[nn] -kb[nn]];
                t=[exp(im*kb[nn]*(D-dx[mm]-da/2)) exp(-im*kb[nn]*(D-dx[mm]-da/2));kb[nn]*exp(im*kb[nn]*(D-dx[mm]-da/2)) -kb[nn]*exp(-im*kb[nn]*(D-dx[mm]-da/2))];
                u=[1 1;ka[nn] -ka[nn]];
                v=[exp(im*ka[nn]*da/2) 0;0 exp(-im*ka[nn]*da/2)];
                T=v*inv(u)*t*inv(s)*r*q*inv(p)*o*inv(n)*m;
            else
                m=[exp(im*ka[nn]*da/2) exp(-im*ka[nn]*da/2);ka[nn]*exp(im*ka[nn]*da/2) -ka[nn]*exp(-im*ka[nn]*da/2)];
                n=[1 1;kb[nn] -kb[nn]];
                o=[exp(im*kb[nn]*db) exp(-im*kb[nn]*db);kb[nn]*exp(im*kb[nn]*db) -kb[nn]*exp(-im*kb[nn]*db)];
                p=[1 1;ka[nn] -ka[nn]];
                q=[exp(im*ka[nn]*(dx[mm]-da/2-db)) exp(-im*ka[nn]*(dx[mm]-da/2-db));ka[nn]*exp(im*ka[nn]*(dx[mm]-da/2-db)) -ka[nn]*exp(-im*ka[nn]*(dx[mm]-da/2-db))];
                r=[1 1;kc[nn] -kc[nn]];
                s=[1+im*tan(kc[nn]*l0)/2 im*tan(kc[nn]*l0)/2;-im*tan(kc[nn]*l0)/2 1-im*tan(kc[nn]*l0)/2];
                t=[1 1;kc[nn] -kc[nn]];
                u=[1 1;ka[nn] -ka[nn]];
                v=[exp(im*ka[nn]*(D-dx[mm])) 0;0 exp(-im*ka[nn]*(D-dx[mm]))];
                T=v*inv(u)*t*s*inv(r)*q*inv(p)*o*inv(n)*m;
            end
            AAA=T[1,1]+T[2,2];
            Temp0=AAA^2-4;
            Temp1=(AAA-sqrt(Temp0))/2;
            r0[mm,nn]=abs((Temp1-T[1,1])/T[1,2]);
            if r0[mm,nn]<=1#导带取一种情况
                Temp1=(AAA-sqrt(Temp0))/2;
                r0[mm,nn]=abs((Temp1-T[1,1])/T[1,2]);
                Arg[mm,nn]=angle((Temp1-T[1,1])/T[1,2]);
            else#此时并不都是禁带，r0仍可能为1，所以取出r0等于1的值与上面导带的情况简并
                Temp1=(AAA+sqrt(Temp0))/2;
                r0[mm,nn]=abs((Temp1-T[1,1])/T[1,2]);
                if r0[mm,nn]<=1.0001&&r0[mm,nn]>=0.9999#1为int所以加一些小数位
                    Temp1=(AAA-sqrt(Temp0))/2;
                    r0[mm,nn]=abs((Temp1-T[1,1])/T[1,2]);
                    Arg[mm,nn]=angle((Temp1-T[1,1])/T[1,2]);
                else
                    Temp1=(AAA+sqrt(Temp0))/2;
                    r0[mm,nn]=abs((Temp1-T[1,1])/T[1,2]);
                    Arg[mm,nn]=angle((Temp1-T[1,1])/T[1,2]);
                end
            end
        end
    end
end
@benchmark testbm()
