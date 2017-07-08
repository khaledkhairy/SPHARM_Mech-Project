function [X C t p] = plot_pretty(obj, nico, option)
%
%   USAGE: S is n x m where n is the number of shapes and m is dimension
%   along which the coefficients are stored (i.e. 3 x that for each
%   coordinate).
%   option: 'red' 'blue' 'volume' 'area' 'random'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
persistent Y_LK_pp
if nargin == 1, nico = 3; option = 'greynone';end
if nargin == 2, option = 'rednone';end

obj = obj.update;
%% using subdivisions of icosahedron
if nico > 6, nico = 6; end;

[X,C]=surface_mesh.sphere_mesh_gen(nico);
%[X, C] = mesh_gen_rand(2);       % looks better than the icosahedron subdivision

[t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
[x y z] = kk_sph2cart(t,p,1);
%%% generate the new basis
[L, K] = obj.indices_gen(1:(obj.L_max + 1)^2); M = length(L);N = length(t);
if size(Y_LK_pp,1)~= N || size(Y_LK_pp,2)~= M
    disp([N M]);
    disp(size(Y_LK_pp));
    pause(3);
    disp('Generating new basis ...');
    
    Y_LK_pp  = zeros(N, M, 'single');
end
for S = 1:length(L),
    Y_LK_pp(:,S) = obj.basis.ylk_bosh(L(S),K(S),p',t')'; % uses bosh version
end;

%%%
[X] = plot_handler(obj, obj.xc, obj.yc, obj.zc, Y_LK_pp, C, option);
view(-90,-90);camlight;
daspect([1 1 1]);axis on;lighting gouraud;view(3);camlight;
% cameramenu;
axis on;ylabel('y');xlabel('x');zlabel('z');
if obj.use_camorbit, for i=1:36,camorbit(10,0,'camera');drawnow;end;end
%axis vis3d;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X] = plot_handler(obj, xclks, yclks, zclks, Y_LK, C, option)

X = Y_LK(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];
%%%%%% look at shape
if isnumeric(option),
    patch('Vertices', X, 'Faces', C,'FaceVertexCData',option*ones(length(X),1),'FaceColor', 'flat','EdgeColor','none','FaceAlpha', 1);
end
if strcmpi(option,'red'),
    patch('Vertices', X, 'Faces', C,'FaceColor', 'r','EdgeColor','k','FaceAlpha',1);
end
if strcmpi(option,'rednone'),
    patch('Vertices', X, 'Faces', C,'FaceColor', 'r','EdgeColor','none','FaceAlpha',1);
end
if strcmpi(option,'rednone_transparent'),
    patch('Vertices', X, 'Faces', C,'FaceColor', 'r','EdgeColor','none','FaceAlpha',0.4);
end
if strcmpi(option,'greennone'),
    patch('Vertices', X, 'Faces', C,'FaceColor', 'g','EdgeColor','none','FaceAlpha',1);
end
if strcmpi(option,'greynone'),
    patch('Vertices', X, 'Faces', C,'FaceColor', [0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',1);
end

if strcmpi(option,'blue'),
    patch('Vertices', X, 'Faces', C,'FaceColor', 'b','EdgeColor','none','FaceAlpha',1);
end
if strcmpi(option,'bluenone'),
    patch('Vertices', X, 'Faces', C,'FaceColor', 'b','EdgeColor','none','FaceAlpha',1);
end
if strcmpi(option,'volume'),
    patch('Vertices', X, 'Faces', C,'FaceVertexCData',obj.V,'FaceColor', 'flat','EdgeColor','none','FaceAlpha', 1);
end
if strcmpi(option,'area'),
    patch('Vertices', X, 'Faces', C,'FaceVertexCData',obj.A,'FaceColor', 'flat','EdgeColor','none','FaceAlpha', 1);
end
if strcmpi(option,'invisible'),
        % do nothing
end
if strcmpi(option,'random')
     a = rand(1,3);      % random RGB color specification
     if sum(a)> 2.8; a = [1 0 1];end
    patch('Vertices', X, 'Faces', C,'FaceColor', a,'EdgeColor','none','FaceAlpha',1);
end
%lighting phong;


