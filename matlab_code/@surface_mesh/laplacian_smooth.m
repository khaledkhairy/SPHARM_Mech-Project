function [X,F, E, L] = laplacian_smooth(X,F, L, maxiter)
%% Smooth a given surface triangulation by locally adjusting vertex
%% positions to the average of neighbors. This function is used as a helper
%% function for preparing meshes for curvature-flow based inflation.

beta = .001;
if nargin<4, maxiter = 10;end
if nargin == 2,
    maxiter = 10;
    %%%%%%%%%%%% generate the cell array of links
    for ix = 1:length(X),   % loop over the vertices
        fmemb = ismember(F, ix);% find the faces that vertex ix belongs to% find the rows in faces where ix occurs
        [ig] = ind2sub(size(fmemb), find(fmemb==1)); % ig is a column vector that indicates in which row of F we can find ix
        links = [];
        for ik = 1:length(ig),
            links = [links F(ig(ik),:)]; %#ok<AGROW>
        end
        L{ix} = unique(links(links~=ix));   % only record the links that are not ix
    end;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

for iter = 1:maxiter,
    for ix = 1:length(X),
        vec = X(ix,:);
        vix = L{ix};
        for ld = 1:30,
            X(ix,:) = X(ix,:) + beta* sum(X(vix,:)-vec(ones(length(vix),1),:))/length(vix);
        end
    end
end