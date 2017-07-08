function [res TP] = self_intersection(obj, nico, TP)
%% detect self-intersection based on intersection of triangles (Moeller method)
%% run with two arguments to generate a database
%% run with three to use an existing TP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1, nico = 3;end

res = 1;    % default is self-intersection = 1

[X,C]=surface_mesh.sphere_mesh_gen(nico);

[t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
% [x y z] = kk_sph2cart(t,p,1);
%%% generate the new basis
[L, K] = obj.indices_gen(1:(obj.L_max + 1)^2); M = length(L);N = length(t);
Y_LK  = zeros(N, M, 'single');
for S = 1:length(L),
    Y_LK(:,S) = obj.basis.ylk_bosh(L(S),K(S),p',t')'; % uses bosh version
end; 
X = Y_LK(:,1:length(obj.xc))* [obj.xc(:) obj.yc(:) obj.zc(:)];

%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 2 || nargin ==1
    fn = sprintf('DATA_self_nico_%d.mat',nico);
    if exist(fn)==2, load(fn,'TP');
    else    % i.e. we are neither given TP nor can we load it then
        %% generate triangle pair list that excludes neighbors
        %% i.e. we assume that locally triangles don't intersect
        %% although the generation of this is slow it needs to be done only once
        TP = [];
        %disp('Building triangle pair database....')
        for ix = 1:size(C,1)
            disp(['Building triangle pair database....  ' num2str(ix) ' of ' num2str(size(C,1))]);
            for jx = ix:size(C,1)
                if ix ~=jx
                    %% test whether triangle ix and triangle jx have one or more common vertices
                    t1 = C(ix,:);
                    t2 = C(jx,:);
                    if ~sum(ismember(t1,t2))
                        TP = [TP;[ix jx]];
                    end
                end
            end
        end
        fn = sprintf('DATA_self_nico_%d.mat',nico);
        str = sprintf('save %s TP;',fn);eval(str);
        disp(['Saved database under filename: ' fn]);
    end
end

% % %% generate triangle pair list that includes neighbors
% % %% problematic as shared vertices are regarded as intersections
% % nT = size(C,1);
% % TP = zeros((nT^2-nT)/2,2, 'uint16');
% % pos = 1;
% % for ix = 1:length(C)
% %     TP(pos:pos+nT-ix-1,1) = ix;
% %     TP(pos:pos+nT-ix-1,2) = ix+1:nT;
% %     pos = pos + nT-ix;
% % end

%%% generate the triangle lists T1 and T2 corresponding to those pairs
nTP = length(TP);
T1 = zeros(nTP,3, 'uint16');
T2 = zeros(nTP,3, 'uint16');

T1= [C(TP(:,1),1) C(TP(:,1),2) C(TP(:,1),3)];
T2= [C(TP(:,2),1) C(TP(:,2),2) C(TP(:,2),3)];

%% generate the vertices (starting at this point we have to have the current X calculated beforehand)
V21 = [X(T2(:,1),1) X(T2(:,1),2) X(T2(:,1),3)]; % coordinates of the first vertex of the T2 triangles
V22 = [X(T2(:,2),1) X(T2(:,2),2) X(T2(:,2),3)];
V23 = [X(T2(:,3),1) X(T2(:,3),2) X(T2(:,3),3)];

V11 = [X(T1(:,1),1) X(T1(:,1),2) X(T1(:,1),3)]; % coordinates of the first vertex of the T1 triangles
V12 = [X(T1(:,2),1) X(T1(:,2),2) X(T1(:,2),3)];
V13 = [X(T1(:,3),1) X(T1(:,3),2) X(T1(:,3),3)];

%% calculate surface normals and distances to construct the plane equations for T2
N2 = cross( (V22-V21), (V23-V21) , 2);  % triangle normal
d2 = -dot(N2,V21,2);                      % to complete the plane equation N2 . X + d2 = 0
%% Calculate distances of vertices of T1 triangles to T2 planes
dV11 = dot(N2,V11,2) + d2;
dV12 = dot(N2,V12,2) + d2;
dV13 = dot(N2,V13,2) + d2;

%% for every pair test whether all vertices are on the plane of the corresponding triangle or not
%%% The point is that if the vertices are not on the plane and all
%%% distances have the same sign, then there is no intersection
Dzero1 = uint8(dV11==0) + uint8(dV12==0) + uint8(dV13==0);  % zero means no point is on the plane
Dsign1 = abs((sign(dV11)) + (sign(dV12)) + (sign(dV13)))~=3;% zero means all signs are the same (i.e.no intersection)
% one means possible intersection

%% calculate surface normals and distances to construct the plane equations for T1
N1 = cross( (V12-V11), (V13-V11) , 2);  % triangle normal
d1 = -dot(N1,V11,2);                      % to complete the plane equation N2 . X + d2 = 0
%% Calculate distances of vertices of T2 triangles to T1 planes
dV21 = dot(N1,V21,2) + d1;
dV22 = dot(N1,V22,2) + d1;
dV23 = dot(N1,V23,2) + d1;
Dzero2 = uint8(dV21==0) + uint8(dV22==0) + uint8(dV23==0);  % zero means no point is on the plane
Dsign2 = abs((sign(dV21)) + (sign(dV22)) + (sign(dV23)))~=3;    % test whether they are all the same sign

% % %%% talk to me
% % if any(Dzero1), disp('Dzero1: found triangle(s) fully coplanar to corresponding triangle');end
% % if any(Dzero2), disp('Dzero2: found triangle(s) fully coplanar to corresponding triangle');end


indx = find(and((Dsign1==1),Dsign2==1)); % get the indices of the possibly intersecting pairs
%%% resolve possible triangle intersection
if ~isempty(indx)
    %disp('possible triangle intersection detected');
    N1I = N1(indx,:);N2I = N2(indx,:);
    DI = cross(N1I,N2I,2);  % this is the line  along the plane of both triangles, i.e. it intersects two edges in each triangle of a pair
    
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
             res = 0; 
             %disp('No intersection based on intervals');
         else
             res = 1;
             %disp('Intersection detected based on interval')
             break;
         end
    end
else
            res = 0;
            disp('No intersection based on preliminary test');
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











