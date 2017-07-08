function plot_H(obj)
obj = obj.props;
patch('Vertices', obj.X, 'Faces', obj.F,...
    'FaceColor', 'interp', 'FaceVertexCData', obj.H, 'FaceAlpha',1,...
    'EdgeColor','none');
axis equal;
%lighting gouraud;camlight;axis vis3d;
if obj.use_camorbit, for i=1:36,camorbit(10,0,'camera');drawnow;end;end