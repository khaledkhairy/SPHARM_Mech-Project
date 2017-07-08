function xc = tr_xc(xc, L_max)


trunc = (L_max+1)^2;
lmax_in = sqrt(length(xc))-1;
if lmax_in<L_max,
    xc(trunc) = 0;
else
    xc = xc(1:trunc);
end
