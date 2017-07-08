function plot_K(obj)
%% plot on top of the shape the scalar field K (Gaussian curvature)
dfig;
obj.needs_updating = 1;
    obj = obj.update_full;
    C = reshape(obj.KG,size(obj.x));
    [~, x y z] = update(obj);
    surf(double(x),double(y),double(z), C,  'EdgeColor', 'none');
    axis equal;
    view(3);
    lighting phong;
    camlight;
    axis off;


    