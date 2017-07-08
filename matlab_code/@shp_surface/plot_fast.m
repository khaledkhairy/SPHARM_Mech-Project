function plot_fast(obj)
[~, x y z] = update(obj);
surf(double(x),double(y),double(z),ones(size(z)).*obj.Eb, 'EdgeColor', obj.edge_color);

