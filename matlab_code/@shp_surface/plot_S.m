function plot_S(obj)
%% plot on top of the shape the scalar field curvedness as defined in Duncan and Olson
dfig;
obj.needs_updating = 1;
obj = obj.update_full;
C = reshape(obj.S,size(obj.x));
[~, x y z] = update(obj);
surf(double(x),double(y),double(z), C,  'EdgeColor', 'none');
axis equal;axis off;
lighting gouraud;
camlight;
if obj.use_camorbit, for i=1:36,camorbit(10,0,'camera');drawnow;end;end
axis vis3d;


