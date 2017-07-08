function [xc, yc, zc, nc] = get_xyz_clks(clks)
% just cut the clks vector into three and return the three vectors
nc = round(length(clks)/3);
xc = clks(1:nc);
yc = clks(nc+1:2*nc);
zc = clks(2*nc+1:3*nc);

xc = xc(:);
yc = yc(:); 
zc = zc(:);
