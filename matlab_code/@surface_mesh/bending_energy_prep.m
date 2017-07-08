function [Eb, H, M, ddA, A, V, n, F_areas] = bending_energy_prep(obj, Co)
%% The routine calculates the Area
%% (A), the volume (V) and the local mean and total mean curvatures (H and h) of the shape.
%%
M = [];
dA_vec = [];
if size(obj.F,2)==4,
    disp('Quad mesh restriction. Properties calculation not implemented yet.');
else
    %if obj.needs_updating
        X = obj.X;
        C = obj.F;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        A = 0; V = 0; H = 0; h = 0; v = 0; F_areas = 0;
        %% quick calculation of the area and volume
        u = X(:,1); v = X(:,2); w = X(:,3);
        x1 = u(C(:,1)); y1 = v(C(:,1));z1 =  w(C(:,1));x2 = u(C(:,2)); y2 = v(C(:,2));z2 =  w(C(:,2));x3 = u(C(:,3)); y3 = v(C(:,3));z3 =  w(C(:,3));
        q = [x2-x1 y2-y1 z2-z1]; r = [x3-x1 y3-y1 z3-z1];
        crossqpr = cross(q,r,2);
        twoA = sqrt(sum(crossqpr.^2,2));   % take the norm
        A = sum(twoA)/2;                    % this is the total area
        F_areas = twoA/2;                   % this is the vector of face areas
        n = crossqpr./twoA(:,ones(3,1));
        V = abs(sum(1/3*(dot(n,[x1 y1 z1], 2).*twoA./2)));
        Vo = 4/3*pi*(A/4/pi)^(3/2);
        v = V/Vo;
        
        %%% for triangle quality determination
        p = [x3-x2 y3-y2 z3-z2];
        d1 = sqrt(q(:,1).^2 + q(:,2).^2 + q(:,3).^2);
        d2 = sqrt(r(:,1).^2 + r(:,2).^2 + r(:,3).^2);
        d3 = sqrt(p(:,1).^2 + p(:,2).^2 + p(:,3).^2);
        quality = 4*F_areas*sqrt(3)./(d1.^2 + d2.^2 + d3.^2);
        %%%% Now calculate the local mean curvature H at a vertex
        %% Remember: Gauss curvature is the angle defect at a vector
        %%          Mean curvature is edge length x dihedral angle at edge
        H = zeros(size(u));
        M = zeros(size(u));
        ddA = zeros(size(u));
        for ix = 1:length(X),   %loop over the vertices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            V_a = X(ix,:);
            %%% To calculate the curvature at some vertex V_a
            %%  we need to find the triangles that it is part of.
            [r c] = ind2sub(size(C),find(C==ix)); % the rows define the triangles indexed into C
            %% Calculate the dihedral angles of all unique permutaions of the
            %% triangle pairs (planes)
            [I2 I1] = ind2sub([length(r) length(r)],find( tril(ones(length(r)),-1)~=0)); % all possible combinations of triangles
            %% I2 and I1 are to be understood as pairs of indices into r. A value
            %% in r is an index into a row in C, which in turn contains indices
            %% into the points (rows) in X.
            % so let us loop over the permutations and calculate the curvature H
            H(ix) = 0;
            theta = [];
            for I = 1:length(I2)
                % each permutation selects a triangle pair
                r1 = r(I2(I));r2 = r(I1(I));  % the two rows(triangles) in C that are considered
                tr1 = C(r1,:);tr2 = C(r2,:); % the two triangles (they contain indices into rows of X)
                % two of these indices are identical (this is the edge, the length
                % of which we need to calculate
                % we know that one of the vertices that are on the edge is V_a
                % now we need to find the other vertex V_e
                tr1r = tr1((tr1~=ix));tr2r = tr2((tr2~=ix));
                rvrs = (length(tr1)-1):-1:1;
                indx = max(tr1r.*(tr1r==tr2r) + tr1r(rvrs).*(tr1r(rvrs)==tr2r));
                if indx >0,  %i.e. if they share an edge
                    V_e = X(indx,:);
                    V_far = X(tr2(tr2~=ix & tr2 ~=indx),:);
                    % the distance between the two vertices of the common edge
                    Lij   = sqrt((V_a(1)-V_e(1))^2 +(V_a(2)-V_e(2))^2+(V_a(3)-V_e(3))^2);
                    %% the surface normals are
                    n1 = n(r1,:);n2 = n(r2,:);
                    theta(I) = acos(dot(n1,n2));
                    %%% we need to determine if they are convex or nonconvex
                    %%% first let us complete the Hessian normal form of the planes
                    %%% of the two intersecting triangles
                    P1 = -(n1(1)*V_a(1) + n1(2)*V_a(2) + n1(3)*V_a(3));
                    P2 = -(n2(1)*V_a(1) + n2(2)*V_a(2) + n2(3)*V_a(3));
                    %%% Now we need to calculate whether the far point of tr2 lies
                    %%% in the half-space of the normal direction (local H is -ve)
                    %%% or on the anti-normal direction (local curvature is +ve)
                    s = sign(dot(n1,V_far)+P1);
                    H(ix) = H(ix) + Lij * real(theta(I))/4 * (s); % as in F.J.thesis
                end
            end
            %     (theta*180/pi)';
            % The average area of triangles around the vertex V_a is given by
            ddA(ix) = sum(F_areas(r))/3;
            if ddA(ix)==0,M(ix) = 0;warning('error in area calculation II');
            else
                H(ix) = H(ix)/ddA(ix);
                M(ix) = H(ix).*ddA(ix);  % this is correct (as in F.J.thesis)
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
%         h = sum(real(-M))./A;
         Eb  = 2* sum((M-Co.*ddA).^2./ddA)/8/pi;   %bending energy
%         dA_vec = ddA;
%         r = sqrt(A/4/pi);
%         dA  = 2 *sum(-M);
%         dAo = 8 * pi * r;
%         da  = dA/dAo;
        
%         
%         obj.A = A;
%         obj.V = V;
%         obj.v = v;
%         obj.Eb = Eb;
%         obj.dA = dA;
%         obj.da = da;
%         obj.F_areas = F_areas;
%         obj.h = h;
%         obj.H = H;
%         obj.quality = quality;
%         obj.needs_updating = 0;
    %end
end