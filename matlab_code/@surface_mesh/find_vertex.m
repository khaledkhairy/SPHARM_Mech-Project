function [ind] = find_vertex(obj, p)
%% returns the index of a vertex closes to this position
thresh = 1e-2;
indx = find(obj.X(:,1)>(p(1)-thresh) & obj.X(:,1)<(p(1)+thresh));
indy = find(obj.X(:,2)>(p(2)-thresh) & obj.X(:,2)<(p(2)+thresh));
indz = find(obj.X(:,3)>(p(3)-thresh) & obj.X(:,3)<(p(3)+thresh));
ind = intersect(indx, indy);
ind = intersect(ind, indz);
