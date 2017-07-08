function write_stl(m, fn)
%%% export the shp shape as stl file
if nargin<2, fn = 'untitled';end
str1 = sprintf('%s.raw',fn);
str2 = sprintf('%s.stl',fn);
m.mesh2raw(double(m.X),double(m.F), str1);
m.raw2stl(str1, str2);