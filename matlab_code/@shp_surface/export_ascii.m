function export_ascii(s, fn)
% example: s = shp_surface('bowling_pin');fn = 'bowling_pin.shp3';export_ascii(s, fn);
%
if nargin<2, fn = 'untitled.shp3';end
n_components = 3 + numel(s.sf);
nc = length(s.X_o)/3;
fid = fopen(fn,'w');
%% print the number of shape vectors given in this file
fprintf(fid,'n_shapes = %d\n',1);
%% print L_max and the number of components of the shape vector 

fprintf(fid, 'L_max = %d\n', s.L_max);
fprintf(fid, 'n_components = %d\n',n_components);

%% print the names of the fields
fprintf(fid,'%s\t%s\t%s','x', 'y', 'z');
for (fix = 1:numel(s.sf))
        field = s.sf{fix};
        fprintf(fid,'\t%s',field{1});
        if length(field{2}.xc)~=nc,
            disp(['Adjusting L_max for field: ' field{1}]);
            %s.sf{fix}{2}= s.sf{fix}{2}.xc(1:nc);
            s.sf{fix}{2}=tr(s.sf{fix}{2},s.L_max);
        end
end
fprintf(fid,'\n');
 
%% print the values
for(ix = 1:nc)
    fprintf(fid,'%e\t%e\t%e', s.xc(ix), s.yc(ix), s.zc(ix));
    for (fix = 1:numel(s.sf))
        field = s.sf{fix};
        fprintf(fid,'\t%e',field{2}.xc(ix));
    end
    fprintf(fid,'\n');
end

    
fclose(fid);