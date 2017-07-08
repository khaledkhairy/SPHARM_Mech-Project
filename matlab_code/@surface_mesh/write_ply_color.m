function write_ply_color(obj,Path, indx)

X = obj.X;
F  = obj.F - 1;
if numel(indx)==1
c = obj.sf{indx}{2};
c = round(mat2gray(c)*255);
c = [c c c];
elseif numel(indx)==2
    c = [round(mat2gray(obj.sf{indx(1)}{2})*255) ...
        round(mat2gray(obj.sf{indx(2)}{2}) * 255) ...
        zeros(size(obj.sf{indx(1)}{2}))];
elseif numel(indx)==3
        c = [round(mat2gray(obj.sf{indx(1)}{2}) * 255) ...
             round(mat2gray(obj.sf{indx(2)}{2}) * 255) ...
             round(mat2gray(obj.sf{indx(3)}{2}) * 255)];
else
    warning('no color information included in the file');
    c = ones(size(X)) * 255;
end
    
[fid,Msg] = fopen(Path,'wt');
if fid == -1, error(Msg); end
fprintf(fid,'ply\nformat ascii 1.0\ncomment created by write_ply_color function in matlab\n');
fprintf(fid,'element vertex %d\n',size(X,1));
fprintf(fid,'property float x\nproperty float y\nproperty float z\n');
fprintf(fid,'property uchar red\nproperty uchar green\nproperty uchar blue\n');
fprintf(fid,'element face %d\n',size(F,1));
% fprintf(fid,'element face %d\n',300);
fprintf(fid,'property list uchar int vertex_indices\nend_header\n');
for i = 1:size(X,1)
    fprintf(fid,'%f %f %f %d %d %d\n', X(i,2), X(i,1), X(i,3), c(i,1), c(i,2), c(i,3));
end
for i = 1:size(F,1)
   fprintf(fid, '3 %d %d %d\n', F(i,1), F(i,2), F(i,3)); 
end
fclose(fid);