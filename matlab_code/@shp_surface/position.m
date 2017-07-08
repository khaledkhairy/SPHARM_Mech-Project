function obj = position(obj, pos)
%% set new translation position according to pos
pos = pos;
Xtmp = obj.X_o;
[xc2 yc2 zc2] = shp_surface.get_xyz_clks(Xtmp);
xc2(1)= pos(1);
yc2(1) = pos(2);
zc2(1) = pos(3);
X_o = [xc2(:)' yc2(:)' zc2(:)'];
obj.X_o = X_o;
obj = update(obj);