function writeVTK(obj, filename,u)
% vtk export
t = obj.X;
p = obj.F;
if nargin>2,
u = obj.sf{u}{2};
else
    u = ones(length(t),1);
end
[np,dim]=size(p);
[nt]=size(t,1);
celltype=[3,5,10];

FID = fopen(strcat(filename,'.vtk'),'w+');
fprintf(FID,'# vtk DataFile Version 2.0\nUnstructured Grid Example\nASCII\n');
fprintf(FID,'DATASET UNSTRUCTURED_GRID\n');

fprintf(FID,'POINTS %d float\n',np);
s='%f %f %f \n';
P=[p zeros(np,3-dim)];
fprintf(FID,s,P');

fprintf(FID,'CELLS %d %d\n',nt,nt*(dim+2));
s='%d ';
for k=1:dim+1
    s=horzcat(s,{' %d'});
end
s=cell2mat(horzcat(s,{' \n'}));
fprintf(FID,s,[(dim+1)*ones(nt,1) t-1]');

fprintf(FID,'CELL_TYPES %d\n',nt);
s='%d\n';
fprintf(FID,s,celltype(dim)*ones(nt,1));

fprintf(FID,'POINT_DATA %s\nSCALARS u float 1\nLOOKUP_TABLE default\n',num2str(np));
s='%f\n';
fprintf(FID,s,u);

fclose(FID);