function [uw_coupling] = si_calc_uw_RWA(t,rf_freq,phimw)

uw_coupling_x = zeros(16);
uw_coupling_y = zeros(16);

    for m = 2:-1:-3
        mw_sigma_x = zeros(16,16);
        mw_sigma_y = zeros(16,16);
        jjj = 13-m;
        kkk = 4-m;
        mw_sigma_x(jjj,kkk) = 1;
        mw_sigma_x(kkk,jjj) = 1;
%         mw_sigma_y(jjj,kkk) = i;
%         mw_sigma_y(kkk,jjj) = -i;
        mw_sigma_y(jjj,kkk) = -i;
        mw_sigma_y(kkk,jjj) = i;
        
        uw_coupling_x = uw_coupling_x+(mw_sigma_x)*(bgrape_ClebschGordan(3,1,4,m,1,m+1)) * ...
            cos( 2*(m-3)*rf_freq*t + phimw); 
        uw_coupling_y = uw_coupling_y+(mw_sigma_y)*(bgrape_ClebschGordan(3,1,4,m,1,m+1)) * ...
            sin( 2*(m-3)*rf_freq*t + phimw);
    end
    
    uw_coupling = uw_coupling_x + uw_coupling_y;
end