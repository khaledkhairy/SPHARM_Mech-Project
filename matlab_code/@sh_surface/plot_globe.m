function plot_globe(obj, nico)

persistent Y_LK_pg

if nargin<2,nico = 4;end;
[X,C]=surface_mesh.sphere_mesh_gen(nico);
[t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
[x y z] = kk_sph2cart(t,p,1);
%%% generate the new basis (at the vertices of the icosahedron)
[L, K] = sh_surface.indices_gen(1:(obj.L_max + 1)^2); M = length(L);N = length(t);
if size(Y_LK_pg,1)~= N || size(Y_LK_pg,2)~= M
    disp('Generating new basis ...');
    Y_LK_pg  = zeros(N, M, 'single');
    for S = 1:length(L),
        Y_LK_pg(:,S) = obj.basis.ylk_bosh(L(S),K(S),p',t')'; % uses bosh version
    end
end
r = Y_LK_pg(:,1:length(obj.xc))* [obj.xc(:)];

patch('Vertices', X, 'Faces', C,'FaceVertexCData',r(:),'FaceColor', 'flat','EdgeColor','none','FaceAlpha', 1);

axis equal;view(3);%lighting phong;camlight;
if obj.use_camorbit, for i=1:36,camorbit(10,0,'camera');drawnow;end;end
axis on;
ylabel('y');
xlabel('x');
zlabel('z');
