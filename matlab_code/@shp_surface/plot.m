function plot(obj)
[~, x y z] = update(obj);
%surf(double(x),double(y),double(z), obj.basis.p-pi, 'EdgeColor', obj.edge_color);
surf(double(x),double(y),double(z), 'EdgeColor', obj.edge_color);
axis equal;
view(3);
lighting gouraud;
camlight;
if obj.use_camorbit, for i=1:36,camorbit(10,0,'camera');drawnow;end;end
axis off;
