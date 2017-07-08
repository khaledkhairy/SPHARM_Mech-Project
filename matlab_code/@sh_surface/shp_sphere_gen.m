function [xclks yclks zclks] = shp_sphere_gen(L_max,R,hfn)
% generates the set of coefficients that correspond to a sphere
% which means filling in the L = 1 coefficients
xclks = zeros(4,1);
yclks = zeros(4,1);
zclks = zeros(4,1);

xclks(4) = R/hfn(1,1);
yclks(2) = R/hfn(1,-1);
zclks(3) = R/hfn(1,0);
