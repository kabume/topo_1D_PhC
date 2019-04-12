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
    k1=((pi./kb+db)/2+da/2)*1e9;
    k1=shaixuank(k1, da, db);
    k_1=((-pi./kb+db)/2+da/2)*1e9;
    k_1=shaixuank(k_1, da, db);
    k0=ones(num1,1)*100;

    c1=atan2(-(tan(kb .* db) .* ka .^ 2 + tan(kb .* db) .* kb .^ 2) ./ ka ./ kb .* (ka .^ 4 .* tan(kb .* db) .^ 2 - kb .^ 4 .* tan(kb .* db) .^ 2 - 0.4e1 .* sqrt(tan(kb .* db) .^ 2 .* ka .^ 4 .* kb .^ 4 + ka .^ 4 .* kb .^ 4)) ./ (ka .^ 4 .* tan(kb .* db) .^ 2 + 0.2e1 .* ka .^ 2 .* kb .^ 2 .* tan(kb .* db) .^ 2 + kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* ka .^ 2 .* kb .^ 2) ./ 0.2e1 + (tan(kb .* db) .* ka .^ 2 - tan(kb .* db) .* kb .^ 2) ./ ka ./ kb ./ 0.2e1, -(ka .^ 4 .* tan(kb .* db) .^ 2 - kb .^ 4 .* tan(kb .* db) .^ 2 - 0.4e1 .* sqrt(tan(kb .* db) .^ 2 .* ka .^ 4 .* kb .^ 4 + ka .^ 4 .* kb .^ 4)) ./ (ka .^ 4 .* tan(kb .* db) .^ 2 + 0.2e1 .* ka .^ 2 .* kb .^ 2 .* tan(kb .* db) .^ 2 + kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* ka .^ 2 .* kb .^ 2)); 
    c2=atan2(-(tan(kb .* db) .* ka .^ 2 + tan(kb .* db) .* kb .^ 2) ./ ka ./ kb .* (ka .^ 4 .* tan(kb .* db) .^ 2 - kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* sqrt(tan(kb .* db) .^ 2 .* ka .^ 4 .* kb .^ 4 + ka .^ 4 .* kb .^ 4)) ./ (ka .^ 4 .* tan(kb .* db) .^ 2 + 0.2e1 .* ka .^ 2 .* kb .^ 2 .* tan(kb .* db) .^ 2 + kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* ka .^ 2 .* kb .^ 2) ./ 0.2e1 + (tan(kb .* db) .* ka .^ 2 - tan(kb .* db) .* kb .^ 2) ./ ka ./ kb ./ 0.2e1, -(ka .^ 4 .* tan(kb .* db) .^ 2 - kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* sqrt(tan(kb .* db) .^ 2 .* ka .^ 4 .* kb .^ 4 + ka .^ 4 .* kb .^ 4)) ./ (ka .^ 4 .* tan(kb .* db) .^ 2 + 0.2e1 .* ka .^ 2 .* kb .^ 2 .* tan(kb .* db) .^ 2 + kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* ka .^ 2 .* kb .^ 2));
    c3=atan2(-(-ka .^ 2 .* tan(kb .* db) - kb .^ 2 .* tan(kb .* db)) ./ ka ./ kb .* (ka .^ 4 .* tan(kb .* db) .^ 2 - kb .^ 4 .* tan(kb .* db) .^ 2 - 0.4e1 .* sqrt(tan(kb .* db) .^ 2 .* ka .^ 4 .* kb .^ 4 + ka .^ 4 .* kb .^ 4)) ./ (ka .^ 4 .* tan(kb .* db) .^ 2 + 0.2e1 .* ka .^ 2 .* kb .^ 2 .* tan(kb .* db) .^ 2 + kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* ka .^ 2 .* kb .^ 2) ./ 0.2e1 + (-ka .^ 2 .* tan(kb .* db) + kb .^ 2 .* tan(kb .* db)) ./ ka ./ kb ./ 0.2e1, -(ka .^ 4 .* tan(kb .* db) .^ 2 - kb .^ 4 .* tan(kb .* db) .^ 2 - 0.4e1 .* sqrt(tan(kb .* db) .^ 2 .* ka .^ 4 .* kb .^ 4 + ka .^ 4 .* kb .^ 4)) ./ (ka .^ 4 .* tan(kb .* db) .^ 2 + 0.2e1 .* ka .^ 2 .* kb .^ 2 .* tan(kb .* db) .^ 2 + kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* ka .^ 2 .* kb .^ 2));
    c4=atan2(-(-ka .^ 2 .* tan(kb .* db) - kb .^ 2 .* tan(kb .* db)) ./ ka ./ kb .* (ka .^ 4 .* tan(kb .* db) .^ 2 - kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* sqrt(tan(kb .* db) .^ 2 .* ka .^ 4 .* kb .^ 4 + ka .^ 4 .* kb .^ 4)) ./ (ka .^ 4 .* tan(kb .* db) .^ 2 + 0.2e1 .* ka .^ 2 .* kb .^ 2 .* tan(kb .* db) .^ 2 + kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* ka .^ 2 .* kb .^ 2) ./ 0.2e1 + (-ka .^ 2 .* tan(kb .* db) + kb .^ 2 .* tan(kb .* db)) ./ ka ./ kb ./ 0.2e1, -(ka .^ 4 .* tan(kb .* db) .^ 2 - kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* sqrt(tan(kb .* db) .^ 2 .* ka .^ 4 .* kb .^ 4 + ka .^ 4 .* kb .^ 4)) ./ (ka .^ 4 .* tan(kb .* db) .^ 2 + 0.2e1 .* ka .^ 2 .* kb .^ 2 .* tan(kb .* db) .^ 2 + kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* ka .^ 2 .* kb .^ 2));

    p110=(ka*da-c3)/2./ka; p111=(ka*da-c3-pi)/2./ka; p112=(ka*da-c3-2*pi)/2./ka; p113=(ka*da-c3-3*pi)/2./ka; p114=(ka*da-c3-4*pi)/2./ka;
    p120=(ka*da-c4)/2./ka; p121=(ka*da-c4-pi)/2./ka; p122=(ka*da-c4-2*pi)/2./ka; p123=(ka*da-c4-3*pi)/2./ka; p124=(ka*da-c4-4*pi)/2./ka;
    
    p110 = shaixuanp1(p110, da, db); p111 = shaixuanp1(p111, da, db); p112 = shaixuanp1(p112, da, db); p113 = shaixuanp1(p113, da, db); p114 = shaixuanp1(p114, da, db);
    p120 = shaixuanp1(p120, da, db); p121 = shaixuanp1(p121, da, db); p122 = shaixuanp1(p122, da, db); p123 = shaixuanp1(p123, da, db); p124 = shaixuanp1(p124, da, db);
    
    p210=(ka.*da+2*ka.*db-c1)./ka/2; p211=(ka.*da+2*ka.*db-c1+pi)./ka/2; p212=(ka.*da+2*ka.*db-c1+2*pi)./ka/2; p213=(ka.*da+2*ka.*db-c1+3*pi)./ka/2; p214=(ka.*da+2*ka.*db-c1+4*pi)./ka/2;
    p220=(ka.*da+2*ka.*db-c2)./ka/2; p221=(ka.*da+2*ka.*db-c2+pi)./ka/2; p222=(ka.*da+2*ka.*db-c2+2*pi)./ka/2; p223=(ka.*da+2*ka.*db-c2+3*pi)./ka/2; p224=(ka.*da+2*ka.*db-c2+4*pi)./ka/2;
    
    p210 = shaixuanp2(p210, da, db); p211 = shaixuanp2(p211, da, db); p212 = shaixuanp2(p212, da, db); p213 = shaixuanp2(p213, da, db); p214 = shaixuanp2(p214, da, db);
    p220 = shaixuanp2(p220, da, db); p221 = shaixuanp2(p221, da, db); p222 = shaixuanp2(p222, da, db); p223 = shaixuanp2(p223, da, db); p224 = shaixuanp2(p224, da, db);

    plot(p110, w, p111, w, p120, w, p121, w)
    hold on
    plot(p210, w, p211, w, p220, w, p221, w)
    plot(k0,w,k1,w,'k',k_1,w,'k')
    xlim([0,200])
    grid on

    function k=shaixuank(k1, da, db)
        for ii=1:length(k1)
            if k1(ii) <= (da/2+db)*1e9 && k1(ii)>=(da/2)*1e9
                k(ii)=k1(ii);
            else
                k(ii)=NaN;
            end
        end
    end

    function p1 = shaixuanp1(p, da, db)
        for ii=1:length(p)
            if p(ii)<=da/2 && p(ii)>=0
                p1(ii)=p(ii) * 1e9;
            else
                p1(ii)=NaN;
            end
        end    
    end

    function p2 = shaixuanp2(p, da, db)
        for ii=1:length(p)
            if p(ii)<=da+db && p(ii)>=db + da/2
                p2(ii)=p(ii) * 1e9;
            else
                p2(ii)=NaN;
            end
        end    
    end