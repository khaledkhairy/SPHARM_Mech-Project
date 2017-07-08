function X_rot = rotate_shp(X_o,ang)
%% rotate the shape described by the spherical harmonic coefficiencts X_o
%% by the Euler angles a b g (radians). Returns the rotated coefficients

a = ang(1);b = ang(2);g = ang(3);
[xc, yc, zc] = shp_surface.get_xyz_clks(X_o);
%% rotation conventions are y-z-y
x = xc(1);xc(1) = 0;
y = yc(1);yc(1) = 0;
z = zc(1);zc(1) = 0;

Rg = rot_mx(g,2);
Rb = rot_mx(b,1);
Ra = rot_mx(a,2);
C = [xc(:)';yc(:)';zc(:)'];
cr = [Ra*Rb*Rg*C];
cr = cr';
tx = [cr(:,1)];tx(1) = x;
ty = [cr(:,2)];ty(1) = y;
tz = [cr(:,3)]; tz(1) = z;
X_rot = [tx; ty; tz ]; 



