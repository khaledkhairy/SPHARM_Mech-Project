function pv = f2p(obj,fv)
% converts the fv values, i.e. values defined on the faces to pv, i.e. the
% average values of the faces at the points
% pv has dimension size(obj.X,1), fv has dimension size(obj.F,1);
obj = edge_info(obj);
for ix = 1:numel(obj.face_memb) % loop over the number of points
    fvec = obj.face_memb{ix};   % get the face indices that each point is member of
    pv(ix) = 0;
    for fix = 1:numel(fvec)
        pv(ix) = pv(ix) + fv(fvec(fix));    % add the contribution of the member face to the point
    end
    pv(ix) = pv(ix)/numel(fvec);            % take the average
end
pv = pv(:);
        