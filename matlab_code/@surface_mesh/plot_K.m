function plot_K(obj)
obj = obj.props;
[A, V, v, F_areas_o, h, H, wb, da, N, K, k_g, dA] = triangulated_props(obj.X, obj.F, 0);
patch('Vertices', obj.X, 'Faces', obj.F,'FaceColor', 'interp', 'FaceVertexCData', K, 'FaceAlpha',1);
axis equal;
lighting gouraud;camlight;axis vis3d;
if obj.use_camorbit, for i=1:36,camorbit(10,0,'camera');drawnow;end;end