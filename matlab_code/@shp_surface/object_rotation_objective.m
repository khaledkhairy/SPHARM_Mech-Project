function res = object_rotation_objective(ang,X_o,D,nY1, nY2, plot_flag)
ang = real(ang);
ang(1) = mod(ang(1), 2*pi);
ang(2) = mod(ang(2), 2*pi);
ang(3) = mod(ang(3), 2*pi);

if plot_flag
    disp(ang);
    X_o1 = rotate_shp(X_o,ang);
    dfig(1);clf;plot_shps(cs2nocs((X_o1)),3,'red');view(3);lighting gouraud;camlight;
end
%% rotate the spherical harmonics according to ang
X_r = shp_surface.rotate_shp(X_o,ang);
[xc yc zc] = shp_surface.get_xyz_clks(X_r);
%% % evaluate the surface vector at theta = 0 (assumed to be aligned with the longest axis of the object's first order ellipsoid)
% t = 0;p = 0;
% X1 = xc(2) * ylk_cos_sin(1,-1,p, t)/N_LK(1,-1) + xc(3) * ylk_cos_sin(1,0,p, t)/N_LK(1,0) + xc(4) * ylk_cos_sin(1,1,p, t)/N_LK(1,1);
% Y1 = yc(2) * ylk_cos_sin(1,-1,p, t)/N_LK(1,-1) + yc(3) * ylk_cos_sin(1,0,p, t)/N_LK(1,0) + yc(4) * ylk_cos_sin(1,1,p, t)/N_LK(1,1);
% Z1 = zc(2) * ylk_cos_sin(1,-1,p, t)/N_LK(1,-1) + zc(3) * ylk_cos_sin(1,0,p, t)/N_LK(1,0) + zc(4) * ylk_cos_sin(1,1,p, t)/N_LK(1,1);
%% evaluate the length of the surface vector at phi = 0 and theta = pi/2 assumed to be aligned with the objects first order ellipsoid shortest axis
% t = pi/2;p = 0;
% X2 = xc(2) * ylk_cos_sin(1,-1,p, t)/N_LK(1,-1) + xc(3) * ylk_cos_sin(1,0,p, t)/N_LK(1,0) + xc(4) * ylk_cos_sin(1,1,p, t)/N_LK(1,1);
% Y2 = yc(2) * ylk_cos_sin(1,-1,p, t)/N_LK(1,-1) + yc(3) * ylk_cos_sin(1,0,p, t)/N_LK(1,0) + yc(4) * ylk_cos_sin(1,1,p, t)/N_LK(1,1);
% Z2 = zc(2) * ylk_cos_sin(1,-1,p, t)/N_LK(1,-1) + zc(3) * ylk_cos_sin(1,0,p, t)/N_LK(1,0) + zc(4) * ylk_cos_sin(1,1,p, t)/N_LK(1,1);

%%% a minimum is achieved if the Northpole (i.e. the point at theta = 0
%%% (and phi = 0)) is along the longest axis of the first-order ellipsoid
%%% is along the z axis
%%% and simultaneously the point where the Greenwich meridian meets the
%%% equator (i.e. phi = 0 and theta = pi/2), is along the shortest axis, is
%%% a long the x axis
Z1 = nY1*zc(2:4);
X2 = nY2*xc(2:4);
res = (abs(Z1)-max(D))^2 + (abs(X2)-min(D))^2;
if plot_flag, disp([abs(Z1) max(D) abs(X2) min(D) res]);end
