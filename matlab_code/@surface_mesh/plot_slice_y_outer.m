function [X_in X_out F] = plot_slice_y_outer(obj, d, yval, cflag)
%% plot the sliced object
if nargin ==1, d = 1;yval = 0;cflag = 1;end
X = obj.X;
F = obj.F;
centers = [];
if nargin <4, cflag = 0;end
Fs = F;
if abs(cflag) == 1,
    colordef black;
    h = figure('Color', [0 0 0]);
    set(h,'InvertHardcopy','off');
else
    colordef white;
    h = figure('Color', [1 1 1]);
    set(h,'InvertHardcopy','on');
end
%% %% translate to center of mass (optional)
% cm = sum(X,1)./size(X,1);
% cm = cm(ones(1,size(X,1)),:);
% X(:,1) = X(:,1)-cm(:,1);
% X(:,2) = X(:,2)-cm(:,2);
% X(:,3) = X(:,3)-cm(:,3);
% if ~isempty(centers),
%     cm = sum(X,1)./size(X,1);
%     cm = cm(ones(1,size(centers,1)),:);
%     centers(:,1) = centers(:,1)-cm(:,1);
%     centers(:,3) = centers(:,3)-cm(:, 3);
% end
%% delete unwanted centers (optional) and plot the centers if they are defined
% if ~isempty(centers)
%     
%     %indx = find(centers(:,2)>yval-3); centers(indx,:) = [];
%     
%     s = shp_surface(1);s.use_camorbit = 0;
%     r = 1;
%     for ix = 1:size(centers,1)
%         s = scale(s,r);
%         s = set_center(s,[centers(ix,:)]);
%         plot_shps(s,2,'greynone');hold on;
%     end
%     
% end
%% calculate normals inner and outer shell
[A, V, v, F_areas, h, H, wb, da, N, K, k_g, dA] = triangulated_props(X, F, 0);
%Xs = X - N .* d/2;
%% find the triangles that intersect that surface
large = 1e3;
X2 = [-large yval -large; 0 yval large; large yval -large];F2 = [1 2 3];
[int_indx, TP] = tri_plane_intersect(X,F, X2, F2, [], 0); % determine the triangle indices that intersect the plane
% patch('vertices',X,'faces',C,'FaceColor','red', 'FaceAlpha', 0.8);  % the middle surface
%% find the surface after cutting
[X F] = get_cut_surface(X, F, int_indx, yval);
[E,L,face_memb] = edge_info(X,F);
%% %% delete orphan points --- faces need to be changed as well
% % % % del_ixX = [];
% % % % for ix = 1:numel(L),
% % % %     t = L{ix};
% % % %     if isempty(t),
% % % %         del_ixX = [del_ixX; ix];
% % % %     end;
% % % % end;
% % % % % X(del_ixX,:) = [];
% % % % N(del_ixX,:) = [];
% % % % L(del_ixX) = [];
% % % % [X, F] = removeisolatednode(X,F);   % also takes care of the face indexing
% % % % [E,L,face_memb] = edge_info(X,F);   % repeat, now that we have new face and vertex indices
% % % % %% plot the points on the cutting plane (optional)
c_ix = find(X(:,2)==yval);  % here are the indices of points along the contour to y = yval
%plot(X(c_ix,1), X(c_ix,3),'.'); 

%% clean the surface of very small triangles (optional);
% [A, V,v, F_areas] = triangulated_props(X, F, 0);
% thresh = 1/1000;
% indx = find(F_areas/norm(F_areas)<thresh);
% disp(['Deleting ' num2str(length(indx)) ' triangles.']);
% F(indx,:) = [];

%% % now that we have the normals at each point we calculate the surface that is a distance d removed from it
X_in = X + N .* d;  %%% now generate the inner surface
X_out = X;% - N .* d/2;  %%% now generate the outer surface

% X_in = X;
% X_out = X-N.*d;

X_in(c_ix,2) = yval;
X_out(c_ix,2) = yval;
%% make the connections to get a cut through the shell thickness
% [X F] = connect_inner_and_outer_shells(X_in, X_out,F, c_ix)

%% define the cut-mesh surface (using the contour indices) and plot it 
view(-180,0);
for ix = 1:length(c_ix),  % loop through the connected pairs
    p1 = X_in(c_ix(ix),:);
    p2 = X_out(c_ix(ix),:);
    %plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)],'k-');hold on;
    linx = L{c_ix(ix)};
    for lix = 1:length(linx)
        v1 = X_in(linx(lix),:);
        v2 = X_out(linx(lix),:);
        if v1(2)==yval,  p3 = v1; p4  = v2;
            %p = [p1;p2;p3;p4];patch(p(:,1), p(:,2), p(:,3),'w','EdgeColor','none');
            p = [p1;p3;p4;p2];patch(p(:,1), p(:,2), p(:,3),'w','EdgeColor','k');
        end
    end
    
end
%%
if sign(cflag) == 1,
%    patch('vertices',X_out,'faces',F,'facecolor','green','edgecolor',[0.2 0.2 0.2], 'FaceAlpha', 1.0);
    patch('vertices',X_out,'faces',F,'facecolor','green','edgecolor','none', 'FaceAlpha', 1.0);
    patch('vertices',X_in,'faces',F,'facecolor','red','edgecolor','none', 'FaceAlpha', 1.0);
    %patch('vertices',Xs,'faces',Fs,'facecolor','green','edgecolor','none','FaceAlpha',0.5);
else
    patch('vertices',X_in,'faces',F,'facecolor','green','edgecolor','none', 'FaceAlpha', 1.0);
    patch('vertices',X_out,'faces',F,'facecolor','red','edgecolor','none', 'FaceAlpha', 1.0);
end

axis tight;axis off;axis equal;
view(136,-2);camlight
lighting phong;camlight
cameramenu
%plot_mesh(X_in,F);hold on;plot_mesh(X_out, F);
%%
function [X F c_ix] = get_cut_surface(X, F, int_indx, yval)
%% of the intersecting triangles, project the vertex that is closest to the plane yval onto yval
contour = [];
c_ix = [];
large = 1e13;
for ix = 1:numel(int_indx),         % loop over intersecting triangles
    %%% find the vertex closest to the plane yval
    tri = F(int_indx(ix),:);
    dmin = large;
    minix = 1;
    for dix = 1:3
        v = X(tri(dix),:);d = abs(v(2)-yval);
        if d<dmin, dmin = d;minix = dix;end
    end
    %%% set the yvalue of that vertex to yval
    v = X(tri(minix),:);
    X(tri(minix),:) = [v(1) yval v(3)];
    %%% store this contour point as a 2-vector
    contour = [contour;v(1) v(3)];   % we only need to store the x and z values
    c_ix = [c_ix;tri(minix)];
end
% %plot(contour(:,1)', contour(:,2)','*');
%%  delete all triangles with vertices of values larger than yval, except those that only have one value
del_tri = [];
ylix_vec = [];
for ix = 1:size(F,1),   % loop over the triangles
    tri = F(ix,:);
    cum = 0;
    ylix = 0;
    for dix = 1:3
        v = X(tri(dix),:);
        if v(2)>yval,
            cum = cum +1;
            ylix = dix;
        end
    end
    if cum>1,
        del_tri = [del_tri;ix];
    elseif cum==1
        %%% the vertex with yval larger is stored in ylix
        ylix_vec = [ylix_vec;tri(ylix)]; % store the index of the point whose y value will be projected
    end
end
F(del_tri,:) = [];
%% now project the one vertex in the triangles identified in the last step
for ix = 1:numel(ylix_vec)
    v1 = X(ylix_vec(ix),:);
    X(ylix_vec(ix),:) = [v1(1) yval v1(3)];
end

%% clean the surface of triangles that are coplanar with the y plane y = yval
del_tri = [];
for ix = 1:size(F,1),   % loop over the triangles
    tri = F(ix,:);
    cum = 0;
    for dix = 1:3
        v = X(tri(dix),:);
        if v(2)== yval, cum = cum +1;end
    end
    if cum == 3, del_tri = [del_tri;ix];end
end
F(del_tri,:) = [];
%%
function [int_indx, TP] = tri_plane_intersect(X1,C1, X2, C2, TP, flag)
%%% find which of the triangles defined in the meshes X1 C1 intersect
%%% with any triangles in X2 C2. Call once to obtain TP and then include TP
%%% in the function call for speedup.
%%%
%%% Output: res = 0 for no intersection, 1 for intersection of at least one triangle pair
%%%   Author: Khaled Khairy --- based on Moeller code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
res = 0;
int_indx = [];
if nargin<6, flag = 0;end
if nargin < 5 || isempty(TP)
    %% generate triangle pair list
    %%% although the generation of this is slow it needs to be done only once
    %%% for other calls include TP in the call.
    TP = [];
    disp('Building triangle pair database....')
    for ix = 1:size(C1,1)
        %disp([num2str(ix) ' of ' num2str(size(C1,1))]);
        TP = [TP;[ix 1]];
    end
    disp('Done!');
end
%% % generate the triangle lists T1 and T2 corresponding to those pairs
nTP = length(TP);
T1 = zeros(nTP,3, 'uint16');
T2 = zeros(nTP,3, 'uint16');

T1= [C1(TP(:,1),1) C1(TP(:,1),2) C1(TP(:,1),3)];
T2= [C2(TP(:,2),1) C2(TP(:,2),2) C2(TP(:,2),3)];

%% generate the vertices (starting at this point we have to have the current X calculated beforehand)
V21 = [X2(T2(:,1),1) X2(T2(:,1),2) X2(T2(:,1),3)]; % coordinates of the first vertex of the T2 triangles
V22 = [X2(T2(:,2),1) X2(T2(:,2),2) X2(T2(:,2),3)];
V23 = [X2(T2(:,3),1) X2(T2(:,3),2) X2(T2(:,3),3)];

V11 = [X1(T1(:,1),1) X1(T1(:,1),2) X1(T1(:,1),3)]; % coordinates of the first vertex of the T1 triangles
V12 = [X1(T1(:,2),1) X1(T1(:,2),2) X1(T1(:,2),3)];
V13 = [X1(T1(:,3),1) X1(T1(:,3),2) X1(T1(:,3),3)];

%% calculate surface normals and distances to construct the plane equations for T2
%N2 = cross( (V22-V21), (V23-V21) , 2);  % triangle normal
N2 = kk_cross((V22-V21),(V23-V21));
%d2 = -dot(N2,V21,2);                      % to complete the plane equation N2 . X + d2 = 0
d2 = -kk_dot(N2,V21);
%% Calculate distances of vertices of T1 triangles to T2 planes
dV11 = kk_dot(N2,V11) + d2;
dV12 = kk_dot(N2,V12) + d2;
dV13 = kk_dot(N2,V13) + d2;

%% for every pair test whether all vertices are on the plane of the corresponding triangle or not
%%% The point is that if the vertices are not on the plane and all
%%% distances have the same sign, then there is no intersection
Dzero1 = uint8(dV11==0) + uint8(dV12==0) + uint8(dV13==0);  % zero means no point is on the plane
Dsign1 = abs((sign(dV11)) + (sign(dV12)) + (sign(dV13)))~=3;% zero means all signs are the same (i.e.no intersection)
% one means possible intersection

%% calculate surface normals and distances to construct the plane equations for T1
N1 = kk_cross( (V12-V11), (V13-V11));  % triangle normal
d1 = -kk_dot(N1,V11);                      % to complete the plane equation N2 . X + d2 = 0
%% Calculate distances of vertices of T2 triangles to T1 planes
dV21 = kk_dot(N1,V21) + d1;
dV22 = kk_dot(N1,V22) + d1;
dV23 = kk_dot(N1,V23) + d1;
Dzero2 = uint8(dV21==0) + uint8(dV22==0) + uint8(dV23==0);  % zero means no point is on the plane
Dsign2 = abs((sign(dV21)) + (sign(dV22)) + (sign(dV23)))~=3;    % test whether they are all the same sign

% % %%% talk to me
% % if any(Dzero1), disp('Dzero1: found triangle(s) fully coplanar to corresponding triangle');end
% % if any(Dzero2), disp('Dzero2: found triangle(s) fully coplanar to corresponding triangle');end


indx = find(and((Dsign1==1),Dsign2==1)); % get the indices of the possibly intersecting pairs
%% % resolve possible triangle intersection
if ~isempty(indx)
    %disp('possible triangle intersection detected');
    N1I = N1(indx,:);N2I = N2(indx,:);
    DI = kk_cross(N1I,N2I);  % this is the line  along the plane of both triangles, i.e. it intersects two edges in each triangle of a pair
    
    V11I = V11(indx,:);V12I = V12(indx,:);V13I = V13(indx,:);
    V21I = V21(indx,:);V22I = V22(indx,:);V23I = V23(indx,:);
    
    dV11I = dV11(indx,:);dV12I = dV12(indx,:);dV13I = dV13(indx,:);
    dV21I = dV21(indx,:);dV22I = dV22(indx,:);dV23I = dV23(indx,:);
    
    %%% now we need to find the pair of points that is on one side and separate
    %%% it from the one point that is on the other side. We will rearrange the
    %%% vertices so that V1 and V2 are on one side and V3 is on the other. We
    %%% will then calculate the intersection
    
    %%% let us do this using for loops first
    for ix = 1:size(DI,1),       % loop over the triangle pairs that will be checked for intersection
        D = DI(ix,:);
        V0 = V11I(ix,:);    % coordinates of the first vertex of the first triangle in the ixth pair to be compared
        V1 = V12I(ix,:);    % coordinates of the second vertex of the first triangle in the ixth pair to be compared
        V2 = V13I(ix,:);    % coordinates of the third vertex of the first triangle in the ixth pair to be compared
        
        U0 = V21I(ix,:);
        U1 = V22I(ix,:);
        U2 = V23I(ix,:);
        
        dv0 = dV11I(ix,:); dv1 = dV12I(ix,:);dv2 = dV13I(ix,:);
        dv0dv1 = dv0*dv1;dv0dv2 = dv0*dv2;
        
        du0 = dV21I(ix,:); du1 = dV22I(ix,:);du2 = dV23I(ix,:);
        du0du1 = du0*du1;du0du2 = du0*du2;
        %/* compute and index to the largest component of D */
        index = find(abs(D) == max(abs(D)));
        index = index(1);
        %/* this is the simplified projection onto L*/
        vp0=V0(index);
        vp1=V1(index);
        vp2=V2(index);
        
        up0=U0(index);
        up1=U1(index);
        up2=U2(index);
        
        %%%%%%
        [isect10 isect11] = COMPUTE_INTERVALS(vp0,vp1,vp2,dv0,dv1,dv2,dv0dv1,dv0dv2);
        %/* compute interval for triangle 2 */
        [isect20 isect21] = COMPUTE_INTERVALS(up0,up1,up2,du0,du1,du2,du0du1,du0du2);
        
        [isect10, isect11] = kk_sort(isect10,isect11);
        [isect20, isect21] = kk_sort(isect20,isect21);
        if(isect11<isect20 || isect21<isect10),
            %res = 0;
            %disp('No intersection based on intervals');
        else
            res = res + 1;
            int_indx = [int_indx;indx(ix)];
            %disp('Intersection detected based on interval')
            if flag, break;end
        end
    end
else
    %res = 0;
    %disp('No intersection based on preliminary test');
end
function [a , b] = kk_sort(a,b)
if a>b, c = a; a = b;b = c;end
function [isect0, isect1] = COMPUTE_INTERVALS(VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2)
if(D0D1>0.0),
    
    %/* here we know that D0D2<=0.0 */                   \
    %/* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
    [isect0 isect1] = ISECT(VV2,VV0,VV1,D2,D0,D1);
    
elseif(D0D2>0.0),
    
    %/* here we know that d0d1<=0.0 */
    [isect0 isect1] = ISECT(VV1,VV0,VV2,D1,D0,D2);
    
elseif(D1*D2>0.0 || D0~=0.0),
    
    %/* here we know that d0d1<=0.0 or that D0!=0.0 */
    [isect0 isect1] =ISECT(VV0,VV1,VV2,D0,D1,D2);
    
elseif(D1~=0.0),
    
    [isect0 isect1] = ISECT(VV1,VV0,VV2,D1,D0,D2);
    
elseif(D2~=0.0),
    
    [isect0 isect1] = ISECT(VV2,VV0,VV1,D2,D0,D1);
    
else
    
    %/* triangles are coplanar */
    disp('triangles are coplanar. Test for intersection not implemented yet');
    isect0 = nan;
    isect1 = nan;
    %return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2);
    
end
function [isect0, isect1] = ISECT(VV0,VV1,VV2,D0,D1,D2)
isect0=VV0+(VV1-VV0)*D0/(D0-D1);
isect1=VV0+(VV2-VV0)*D0/(D0-D2);








