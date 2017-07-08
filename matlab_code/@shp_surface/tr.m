function X_o = tr(clks, L_max)

    %clks = S.X_o;% check whether S has a dimension which is one. i.e. we have only one shape
    nc = round(length(clks)/3);
    xclks = clks(1:nc);yclks = clks(nc+1:2*nc); zclks = clks(2*nc+1:3*nc);
    trunc = (L_max+1)^2;
    lmax_in = sqrt(length(clks)/3)-1;
    if lmax_in<L_max,
        xclks(trunc) = 0;
        yclks(trunc) = 0;
        zclks(trunc) = 0;
    else
        xclks = xclks(1:trunc);
        yclks = yclks(1:trunc);
        zclks = zclks(1:trunc);
    end
    X_o = [xclks(:)' yclks(:)' zclks(:)']';
