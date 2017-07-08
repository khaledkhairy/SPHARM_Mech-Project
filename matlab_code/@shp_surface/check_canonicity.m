function pass = check_canonicity(X_o)
verbose = 0;pass = 0;
[xc yc zc] = shp_surface.get_xyz_clks(X_o);
if find(abs(xc(2:4))==max(abs(xc(2:4))))==3 && find(abs(yc(2:4))==max(abs(yc(2:4))))==1 &&...
        find(abs(zc(2:4))==max(abs(zc(2:4))))==2
    if verbose, disp('Validating -- Major elements found on diagonal');end
    if sign(xc(4))== -1 && sign(yc(2))==1 && sign(zc(3))==-1,
        if verbose, disp('Validating -- Major elements have correct signs');end
        pass = 1;
    end
end