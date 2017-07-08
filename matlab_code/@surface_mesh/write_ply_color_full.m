function write_ply_color(obj,Path)

X = obj.X;
F  = obj.F - 1;

[fid,Msg] = fopen(Path,'wt');
if fid == -1, error(Msg); end
fprintf(fid,'ply\nformat ascii 1.0\ncomment created by write_ply_color function in matlab\n');
fprintf(fid,'element vertex %d\n',size(X,1));
fprintf(fid,'property float x\nproperty float y\nproperty float z\n');
% fprintf(fid,'property uchar red\nproperty uchar green\nproperty uchar blue\n');
for i = 1:length(obj.sf)
    fprintf(fid,'property uchar %s\n', obj.sf{i}{1});
end

fprintf(fid,'element face %d\n',size(F,1));
% fprintf(fid,'element face %d\n',300);
fprintf(fid,'property list uchar int vertex_indices\nend_header\n');
for i = 1:size(X,1)
    fprintf(fid,'%f %f %f ', -X(i,1), X(i,2), X(i,3));
    for k = 1:length(obj.sf)
       fprintf(fid,'%d ',obj.sf{k}{2}(i));
    end
    fprintf(fid,'\n');
end
for i = 1:size(F,1)
   fprintf(fid, '3 %d %d %d\n', F(i,1), F(i,2), F(i,3)); 
end
fclose(fid);