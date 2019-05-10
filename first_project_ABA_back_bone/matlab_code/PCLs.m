%p123代表3阶，第2个解，A层左半部分
function PCLs(da, db, na, nb, w)
    c=3e8;
    ka=w*1e15*na/c;
    kb=w*1e15*nb/c;
    num1=length(w);
    k1=((pi./kb+db)/2+da/2)*1e9;
    k1=shaixuank(k1, da, db);
    k_1=((-pi./kb+db)/2+da/2)*1e9;
    k_1=shaixuank(k_1, da, db);
    k2=((2*pi./kb+db)/2+da/2)*1e9;
    k2=shaixuank(k2, da, db);
    k_2=((-2*pi./kb+db)/2+da/2)*1e9;
    k_2=shaixuank(k_2, da, db);
    k0=ones(num1,1)*(da + db) * 1e9 / 2 ;

    c1=atan2(-(tan(kb .* db) .* ka .^ 2 + tan(kb .* db) .* kb .^ 2) ./ ka ./ kb .* (ka .^ 4 .* tan(kb .* db) .^ 2 - kb .^ 4 .* tan(kb .* db) .^ 2 - 0.4e1 .* sqrt(tan(kb .* db) .^ 2 .* ka .^ 4 .* kb .^ 4 + ka .^ 4 .* kb .^ 4)) ./ (ka .^ 4 .* tan(kb .* db) .^ 2 + 0.2e1 .* ka .^ 2 .* kb .^ 2 .* tan(kb .* db) .^ 2 + kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* ka .^ 2 .* kb .^ 2) ./ 0.2e1 + (tan(kb .* db) .* ka .^ 2 - tan(kb .* db) .* kb .^ 2) ./ ka ./ kb ./ 0.2e1, -(ka .^ 4 .* tan(kb .* db) .^ 2 - kb .^ 4 .* tan(kb .* db) .^ 2 - 0.4e1 .* sqrt(tan(kb .* db) .^ 2 .* ka .^ 4 .* kb .^ 4 + ka .^ 4 .* kb .^ 4)) ./ (ka .^ 4 .* tan(kb .* db) .^ 2 + 0.2e1 .* ka .^ 2 .* kb .^ 2 .* tan(kb .* db) .^ 2 + kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* ka .^ 2 .* kb .^ 2)); 
    c2=atan2(-(tan(kb .* db) .* ka .^ 2 + tan(kb .* db) .* kb .^ 2) ./ ka ./ kb .* (ka .^ 4 .* tan(kb .* db) .^ 2 - kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* sqrt(tan(kb .* db) .^ 2 .* ka .^ 4 .* kb .^ 4 + ka .^ 4 .* kb .^ 4)) ./ (ka .^ 4 .* tan(kb .* db) .^ 2 + 0.2e1 .* ka .^ 2 .* kb .^ 2 .* tan(kb .* db) .^ 2 + kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* ka .^ 2 .* kb .^ 2) ./ 0.2e1 + (tan(kb .* db) .* ka .^ 2 - tan(kb .* db) .* kb .^ 2) ./ ka ./ kb ./ 0.2e1, -(ka .^ 4 .* tan(kb .* db) .^ 2 - kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* sqrt(tan(kb .* db) .^ 2 .* ka .^ 4 .* kb .^ 4 + ka .^ 4 .* kb .^ 4)) ./ (ka .^ 4 .* tan(kb .* db) .^ 2 + 0.2e1 .* ka .^ 2 .* kb .^ 2 .* tan(kb .* db) .^ 2 + kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* ka .^ 2 .* kb .^ 2));
    c3=atan2(-(-ka .^ 2 .* tan(kb .* db) - kb .^ 2 .* tan(kb .* db)) ./ ka ./ kb .* (ka .^ 4 .* tan(kb .* db) .^ 2 - kb .^ 4 .* tan(kb .* db) .^ 2 - 0.4e1 .* sqrt(tan(kb .* db) .^ 2 .* ka .^ 4 .* kb .^ 4 + ka .^ 4 .* kb .^ 4)) ./ (ka .^ 4 .* tan(kb .* db) .^ 2 + 0.2e1 .* ka .^ 2 .* kb .^ 2 .* tan(kb .* db) .^ 2 + kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* ka .^ 2 .* kb .^ 2) ./ 0.2e1 + (-ka .^ 2 .* tan(kb .* db) + kb .^ 2 .* tan(kb .* db)) ./ ka ./ kb ./ 0.2e1, -(ka .^ 4 .* tan(kb .* db) .^ 2 - kb .^ 4 .* tan(kb .* db) .^ 2 - 0.4e1 .* sqrt(tan(kb .* db) .^ 2 .* ka .^ 4 .* kb .^ 4 + ka .^ 4 .* kb .^ 4)) ./ (ka .^ 4 .* tan(kb .* db) .^ 2 + 0.2e1 .* ka .^ 2 .* kb .^ 2 .* tan(kb .* db) .^ 2 + kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* ka .^ 2 .* kb .^ 2));
    c4=atan2(-(-ka .^ 2 .* tan(kb .* db) - kb .^ 2 .* tan(kb .* db)) ./ ka ./ kb .* (ka .^ 4 .* tan(kb .* db) .^ 2 - kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* sqrt(tan(kb .* db) .^ 2 .* ka .^ 4 .* kb .^ 4 + ka .^ 4 .* kb .^ 4)) ./ (ka .^ 4 .* tan(kb .* db) .^ 2 + 0.2e1 .* ka .^ 2 .* kb .^ 2 .* tan(kb .* db) .^ 2 + kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* ka .^ 2 .* kb .^ 2) ./ 0.2e1 + (-ka .^ 2 .* tan(kb .* db) + kb .^ 2 .* tan(kb .* db)) ./ ka ./ kb ./ 0.2e1, -(ka .^ 4 .* tan(kb .* db) .^ 2 - kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* sqrt(tan(kb .* db) .^ 2 .* ka .^ 4 .* kb .^ 4 + ka .^ 4 .* kb .^ 4)) ./ (ka .^ 4 .* tan(kb .* db) .^ 2 + 0.2e1 .* ka .^ 2 .* kb .^ 2 .* tan(kb .* db) .^ 2 + kb .^ 4 .* tan(kb .* db) .^ 2 + 0.4e1 .* ka .^ 2 .* kb .^ 2));

    p110=(ka*da-c3)/2./ka; p112=(ka*da-c3-2*pi)/2./ka; p114=(ka*da-c3-4*pi)/2./ka;
    p120=(ka*da-c4)/2./ka; p122=(ka*da-c4-2*pi)/2./ka; p124=(ka*da-c4-4*pi)/2./ka;
    
    p110 = shaixuanp1(p110, da, db); p112 = shaixuanp1(p112, da, db); p114 = shaixuanp1(p114, da, db);
    p120 = shaixuanp1(p120, da, db); p122 = shaixuanp1(p122, da, db); p124 = shaixuanp1(p124, da, db);
    
    p210=(ka.*da+2*ka.*db-c1)./ka/2; p212=(ka.*da+2*ka.*db-c1+2*pi)./ka/2; p214=(ka.*da+2*ka.*db-c1+4*pi)./ka/2;
    p220=(ka.*da+2*ka.*db-c2)./ka/2; p222=(ka.*da+2*ka.*db-c2+2*pi)./ka/2; p224=(ka.*da+2*ka.*db-c2+4*pi)./ka/2;
    
    p210 = shaixuanp2(p210, da, db); p212 = shaixuanp2(p212, da, db); p214 = shaixuanp2(p214, da, db);
    p220 = shaixuanp2(p220, da, db); p222 = shaixuanp2(p222, da, db); p224 = shaixuanp2(p224, da, db);

    
    
    figure
    plot(p110, w, 'r',p112, w, 'g--', p114, w, 'b:', p120, w, 'r', p122, w,'g--', p124, w, 'b:', 'LineWidth',2);
    legend("0", "pi", "2pi")
    hold on
    plot(p210, w, 'r',p212, w, 'g--', p214, w, 'b:', p220, w, 'r', p222, w,'g--', p224, w, 'b:', 'LineWidth',2);
    plot(k0, w, 'k', k1, w, 'g--', k_1, w, 'g--', k2, w,'b:', k_2, w, 'b:', 'LineWidth',2);
    function k=shaixuank(k1, da, db)
        %k = zeros(length(k1));
        for ii=1:length(k1)
            if k1(ii) <= (da/2+db)*1e9 && k1(ii)>=(da/2)*1e9
                k(ii)=k1(ii);
            else
                k(ii)=NaN;
            end
        end
    end

    function p1 = shaixuanp1(p, da, ~)
        %p1 = zeros(length(p));
        for ii=1:length(p)
            if p(ii)<=da/2 && p(ii)>=0
                p1(ii)=p(ii) * 1e9;
            else
                p1(ii)=NaN;
            end
        end    
    end

    function p2 = shaixuanp2(p, da, db)
        %p2 = zeros(length(p));
        for ii=1:length(p)
            if p(ii)<=da+db && p(ii)>=db + da/2
                p2(ii)=p(ii) * 1e9;
            else
                p2(ii)=NaN;
            end
        end    
    end
end