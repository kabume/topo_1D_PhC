function [r0, Arg] = ABA(da, db, lc, na, nb, nc, w, dx)
    c = 3e8;
    D = da + db;
    ka = w * 1e15 * na / c;
    kb = w * 1e15 * nb / c;
    kc = w * 1e15 * nc / c;
    num1 = length(w);
    num2 = length(dx);
    r0=zeros(num2,num1);
    Arg=zeros(num2,num1);
    parfor nn = 1:num1
        for mm = 1:num2
            if mm < da / 2 / D * num2
                m = [exp(1i*ka(nn)*dx(mm)) exp(-1i*ka(nn)*dx(mm)); ka(nn)*exp(1i*ka(nn)*dx(mm)) -ka(nn)*exp(-1i*ka(nn)*dx(mm))];
                n = [1 1; kc(nn) -kc(nn)];
                o = [1+1i*tan(kc(nn)*lc)/2 1i*tan(kc(nn)*lc)/2; -1i*tan(kc(nn)*lc)/2 1-1i*tan(kc(nn)*lc)/2];
                p = [1 1; kc(nn) -kc(nn)];
                q = [1 1; ka(nn) -ka(nn)];
                r = [exp(1i*ka(nn)*(da/2-dx(mm))) exp(-1i*ka(nn)*(da/2-dx(mm))); ka(nn)*exp(1i*ka(nn)*(da/2-dx(mm))) -ka(nn)*exp(-1i*ka(nn)*(da/2-dx(mm)))];
                s = [1 1; kb(nn) -kb(nn)];
                t = [exp(1i*kb(nn)*db) exp(-1i*kb(nn)*db); kb(nn)*exp(1i*kb(nn)*db) -kb(nn)*exp(-1i*kb(nn)*db)];
                u = [1 1; ka(nn) -ka(nn)];
                v = [exp(1i*ka(nn)*da/2) 0; 0 exp(-1i*ka(nn)*da/2)];
                T = v * inv(u) * t * inv(s) * r * inv(q) * p * o * inv(n) * m;
            elseif mm >= da / 2 / D * num2 && mm < (da/2+db) / D * num2
                m = [exp(1i*ka(nn)*da/2) exp(-1i*ka(nn)*da/2); ka(nn)*exp(1i*ka(nn)*da/2) -ka(nn)*exp(-1i*ka(nn)*da/2)];
                n = [1 1; kb(nn) -kb(nn)];
                o = [exp(1i*kb(nn)*(dx(mm)-da/2)) exp(-1i*kb(nn)*(dx(mm)-da/2)); kb(nn)*exp(1i*kb(nn)*(dx(mm)-da/2)) -kb(nn)*exp(-1i*kb(nn)*(dx(mm)-da/2))];
                p = [1 1;kc(nn) -kc(nn)];
                q = [1+1i*tan(kc(nn)*lc)/2 1i*tan(kc(nn)*lc)/2; -1i*tan(kc(nn)*lc)/2 1-1i*tan(kc(nn)*lc)/2];
                r = [1 1; kc(nn) -kc(nn)];
                s = [1 1; kb(nn) -kb(nn)];
                t = [exp(1i*kb(nn)*(D-dx(mm)-da/2)) exp(-1i*kb(nn)*(D-dx(mm)-da/2)); kb(nn)*exp(1i*kb(nn)*(D-dx(mm)-da/2)) -kb(nn)*exp(-1i*kb(nn)*(D-dx(mm)-da/2))];
                u = [1 1; ka(nn) -ka(nn)];
                v = [exp(1i*ka(nn)*da/2) 0; 0 exp(-1i*ka(nn)*da/2)];
                T = v * inv(u) * t * inv(s) * r * q * inv(p) * o * inv(n) * m;
            else
                m = [exp(1i*ka(nn)*da/2) exp(-1i*ka(nn)*da/2); ka(nn)*exp(1i*ka(nn)*da/2) -ka(nn)*exp(-1i*ka(nn)*da/2)];
                n = [1 1; kb(nn) -kb(nn)];
                o = [exp(1i*kb(nn)*db) exp(-1i*kb(nn)*db); kb(nn)*exp(1i*kb(nn)*db) -kb(nn)*exp(-1i*kb(nn)*db)];
                p = [1 1; ka(nn) -ka(nn)];
                q = [exp(1i*ka(nn)*(dx(mm)-da/2-db)) exp(-1i*ka(nn)*(dx(mm)-da/2-db)); ka(nn)*exp(1i*ka(nn)*(dx(mm)-da/2-db)) -ka(nn)*exp(-1i*ka(nn)*(dx(mm)-da/2-db))];
                r = [1 1; kc(nn) -kc(nn)];
                s = [1+1i*tan(kc(nn)*lc)/2 1i*tan(kc(nn)*lc)/2; -1i*tan(kc(nn)*lc)/2 1-1i*tan(kc(nn)*lc)/2];
                t = [1 1; kc(nn) -kc(nn)];
                u = [1 1; ka(nn) -ka(nn)];
                v = [exp(1i*ka(nn)*(D-dx(mm))) 0; 0 exp(-1i*ka(nn)*(D-dx(mm)))];
                T = v * inv(u) * t * s * inv(r) * q * inv(p) * o * inv(n) * m;
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
end