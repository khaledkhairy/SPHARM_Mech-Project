function export_stl(s,fn)
%%% export the shp shape as stl file
if nargin<2, fn = s.name;end
m = get_mesh(s);
write_stl(m, fn);

