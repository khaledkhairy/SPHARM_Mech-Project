function [XF X C Y_LK_out t p] = get_mesh(obj, nico, Y_LK_in, C)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
persistent Y_LK
if nargin == 1, nico = 3;end
obj = obj.update;
%% using subdivisions of icosahedron
%if nico > 6, nico = 6; disp(['Icosahedron subdivision: ' num2str(nico)]);end;
if nargin<3 % then build the basis
    %%%%[X,C]=surface_mesh.sphere_mesh_gen(nico);
    [X,C] = mesh_gen(nico,0);
    [t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
    [x y z] = kk_sph2cart(t,p,1);
    %%% generate the new basis (at the vertices of the icosahedron)
    [L, K] = shp_surface.indices_gen(1:(obj.L_max + 1)^2); M = length(L);N = length(t);
    if size(Y_LK,1)~= N || size(Y_LK,2)~= M
        disp(['Required dimensions: ' num2str(N) ' x ' num2str(M)]);
        disp(['Current dimensions (in memory): ' num2str(size(Y_LK))]);
        pause(3);
        disp('Generating new basis ...');
        
        Y_LK  = zeros(N, M, 'single');
        for S = 1:length(L),
            Y_LK(:,S) = obj.basis.ylk_bosh(L(S),K(S),p',t')'; % uses bosh version
        end;
    end
    X = Y_LK(:,1:length(obj.xc))* [obj.xc(:) obj.yc(:) obj.zc(:)];
    XF = surface_mesh(X,C);
    if ~isempty(obj.sf)
        %%%% loop over the scalar fields and add them to the surface mesh
        %%%% object
        for ix = 1:numel(obj.sf)
            if length(obj.sf{ix}{2}.xc)~=length(obj.xc),
                disp(['Adjusting L_max for field: ' obj.sf{ix}{1}]);
                obj.sf{ix}{2}.xc= trunc_sh(obj.sf{ix}{2}.xc, obj.L_max);
            end
            s = obj.sf{ix};
            XF.sf{ix} = {s{1}, Y_LK(:,1:length(obj.xc))* [s{2}.xc]};
        end
    end
    XF.t = t;
    XF.p = p;
    Y_LK_out = Y_LK;
else
    
    %disp([size(Y_LK) size(C) nico]);
    
    X = Y_LK_in(:,1:length(obj.xc))* [obj.xc(:) obj.yc(:) obj.zc(:)];
    XF = [];
    Y_LK_out = Y_LK_in;
end

%%%%%%%%%%%%%%%%%%%%%%%%
function [X, F,x, y, z,  A, V, v, F_areas, h, H, Eb, da] = mesh_gen(dim, flag)
A = 12.5664;

if(flag),
    %%%% old approximate method
    P = partsphere(dim^2);
    x = reshape(P(1,:),dim,dim);
    y = reshape(P(2,:),dim,dim);
    z = reshape(P(3,:),dim,dim);
else
    % % %%%% using subdivisions of icosahedron
    if dim>6, dim = 6;end
    [X,F]=BuildSphere(dim);
    [t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
    [x y z] = kk_sph2cart(t,p,sqrt(A/4/pi));
end

X = double([x(:) y(:) z(:)]);


[C] = convhulln(X, {'Qt'});c = reshape(C,length(C)*3,1);c = unique(c);c = sort(c);
X = X(c,:);
whos X x
C = convhulln(X, {'Qt'});c = reshape(C,length(C)*3,1);c = unique(c);c = sort(c);
F = C;

[A, V, v, F_areas, h, H, Eb, da] = triangulated_props(X, F, 0);

x = X(:,1);
y = X(:,2);
z = X(:,3);

