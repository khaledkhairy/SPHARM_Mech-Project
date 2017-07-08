function plot_field(obj, nsf)
if nargin == 1, nsf = 1;end
if isempty(obj.sf),
    plot(obj);
else
    s = obj.sf{nsf};
    s = s{2};
    patch('Vertices', obj.X, 'Faces', obj.F,'FaceVertexCData',s(:),'FaceColor', 'flat','EdgeColor','none','FaceAlpha', 1);
    axis on;axis equal;xlabel('x');ylabel('y');zlabel('z');
end