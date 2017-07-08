function [XF X C Y_LK t p] = get_mesh(obj, nico, Y_LK, C)
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 1, nico = 3;end
obj = obj.update;
%% using subdivisions of icosahedron
%if nico > 6, nico = 6; disp(['Icosahedron subdivision: ' num2str(nico)]);end;
if nargin<3 % then build the basis
    [X,C]=surface_mesh.sphere_mesh_gen(nico);
    [t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
    [x y z] = kk_sph2cart(t,p,1);
    %%% generate the new basis (at the vertices of the icosahedron)
    [L, K] = sh_surface.indices_gen(1:(obj.L_max + 1)^2); M = length(L);N = length(t);
    Y_LK  = zeros(N, M, 'single');
    for S = 1:length(L),
        Y_LK(:,S) = obj.basis.ylk_bosh(L(S),K(S),p',t')'; % uses bosh version
    end;
    r = Y_LK(:,1:length(obj.xc))* [obj.xc(:)];
    [x y z] = kk_sph2cart(t,p,r);
    X = [x(:);y(:);z(:)];
    XF = surface_mesh(X,C);
else
    r = Y_LK(:,1:length(obj.xc))* [obj.xc(:)];
    [x y z] = kk_sph2cart(t,p,r);
    X = [x(:);y(:);z(:)];
    XF = surface_mesh(X,C);
end


