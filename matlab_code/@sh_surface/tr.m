function S = tr(S, L_max)


    clks = S.xc;% check whether S has a dimension which is one. i.e. we have only one shape
    nc = round(length(clks));
    xclks = clks(1:nc);
    trunc = (L_max+1)^2;
    lmax_in = sqrt(length(clks)/3)-1;
    if lmax_in<L_max,
        xclks(trunc) = 0;
    else
        xclks = xclks(1:trunc);
    end
    S.xc = xclks(:);
    S.L_max = L_max;
    S.basis = sh_basis(L_max,S.basis.gdim);
