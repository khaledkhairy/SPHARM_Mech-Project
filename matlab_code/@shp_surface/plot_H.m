function plot_H(obj)
%% plot on top of the shape the scalar field H
% dfig;
obj.needs_updating = 1;
obj = update_full(obj);
C = reshape(obj.H,size(obj.x));
% %     minH = -10 ;  C(C<minH) = minH;
% %     maxH = 10 ;  C(C>maxH) = maxH;
[~, x y z] = update(obj);
lo = -0.1;hi = 0.1;
C(C<lo) = lo;
C(C>hi) = hi;
surf(double(x),double(y),double(z), C,  'EdgeColor', 'none');
axis equal;axis off;
% lighting gouraud;camlight;
if obj.use_camorbit, for i=1:36,camorbit(10,0,'camera');drawnow;end;end
% axis vis3d;


