function [X C sf Y1 Y2] = plot_field(obj, nico, nsf, subflag, Y1, Y2)
% subflag =1 means we use icosahedron subdivision
% subflag =0 means uniform random points (nico is used as dim in this case
%  nsf can be the index of the required existing field
% or a sh_surface object
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
persistent Y_LK_pf Y_LK_sf
if nargin ==6,
    Y_LK_pf = Y1;
    Y_LK_sf = Y2;
end
if nargin == 1, nico = 4;nsf = 1;subflag = 1;
elseif nargin ==2,
    nsf = nico;nico = 3;subflag = 1;
elseif nargin==3,subflag = 1;
end

obj = obj.update;
if ~strcmp('sh_surface',class(nsf)), % i.e. if nsf is not a sh_surface object (then it is assumed to be the field index/indices to be plotted)
    if length(nsf) == 1
        sf = obj.sf{nsf};disp(['Plotting scalar field:  ' sf{1}]);
        s = sf{2};%
    else
        for ix = 1:length(nsf)
            sf = obj.sf{nsf(ix)};disp(['Plotting scalar field:  ' sf{1}]);
            s{ix} = sf{2};%
        end
    end
    field_name = sf{1};
else   %nsf is in fact an sh_surface object
    s = nsf;
end
%% generate subdivisions of icosahedron
if subflag==1,
    if nico > 6, nico = 6; disp(['Icosahedron subdivision: ' num2str(nico)]);end;
    [X,C]=surface_mesh.sphere_mesh_gen(nico);
else
    [X, C] = mesh_gen_rand(nico);       % looks better than the icosahedron subdivision
end
[t p] = kk_cart2sph(X(:,1),X(:,2),X(:,3));
% [x y z] = kk_sph2cart(t,p,1);
%% % generate the new basis for shape_outline
[L, K] = obj.indices_gen(1:(obj.L_max + 1)^2); M = length(L);N = length(t);
if size(Y_LK_pf,1)~= N || size(Y_LK_pf,2)~= M
    disp([N M]);
    disp(size(Y_LK_pf));
    pause(3);
    disp('Generating new basis ...');
    
    Y_LK_pf  = zeros(N, M, 'single');
    for S = 1:length(L),
        disp(['Vector ' num2str(S) ' of ' num2str(length(L))]);
        Y_LK_pf(:,S) = obj.basis.ylk_bosh(L(S),K(S),p',t')'; % uses bosh version
    end
end
%% generate the shape outline
X = Y_LK_pf(:,1:length(obj.xc))* [obj.xc(:) obj.yc(:) obj.zc(:)];

%% % generate the new basis for scalar field
if length(nsf) == 1
[L, K] = obj.indices_gen(1:(s.L_max + 1)^2); M = length(L);N = length(t);
else
    [L, K] = obj.indices_gen(1:(s{1}.L_max + 1)^2); M = length(L);N = length(t);
end
if size(Y_LK_sf,1)~= N || size(Y_LK_sf,2)~= M
    disp([N M]);
    disp(size(Y_LK_sf));
    %     pause(3);
    disp('Generating new basis ...');
    Y_LK_sf  = zeros(N, M, 'single');
    for S = 1:length(L),
        disp(['Vector ' num2str(S) ' of ' num2str(length(L))]);
        Y_LK_sf(:,S) = sh_basis.ylk_bosh(L(S),K(S),p',t')'; % uses bosh version
    end;
end
if length(nsf)==1,
    
    sf = Y_LK_sf(:,1:length(s.xc))*s.xc(:);%% generate the  scalar field
else
    counter = 0;
    for ix = [nsf]%1:length(nsf)
        counter = counter + 1;
        %%%% uncomment when using scalar fields with varying L_max
% % %         [L, K] = obj.indices_gen(1:(s{counter}.L_max + 1)^2); M = length(L);N = length(t);
% % %         Y_LK_sf  = zeros(N, M, 'single');
% % %         for S = 1:length(L),
% % %             Y_LK_sf(:,S) = sh_basis.ylk_bosh(L(S),K(S),p',t')'; % uses bosh version
% % %         end;
        sf{counter} = Y_LK_sf(:,1:length(s{counter}.xc))*s{counter}.xc(:);
    end
end
%% do the actual plotting
if length(nsf)>1,
    %     sf = [sf{1}(:) sf{2}(:) sf{3}(:)];
    
    sf = [sf{(1)}(:) sf{(2)}(:) sf{(3)}(:)];%zeros(size(sf{1}(:)))]; % twist in green, ftz in red
    sf(:,1) = sf(:,1).* (sf(:,1)>50);
    sf = [ sf(:,2) 255-sf(:,1) 255-sf(:,1)]/max(sf(:));
    % elseif nargin == 4
    %     if color ==1,     sf = [sf(:) zeros(size(sf(:))) zeros(size(sf(:)))]/max(sf(:));
    %     elseif color ==2, sf = [zeros(size(sf(:))) sf(:) zeros(size(sf(:)))]/max(sf(:));
    %     elseif color ==3, sf = [zeros(size(sf(:))) zeros(size(sf(:))) sf(:)]/max(sf(:));
    %     end
end

%%% load
%%% C:\KK_share\Projects\hot\Project_modeling_embryogenesis\code\modeling_drosophila_gastrulation_02\matlab.mat
%%% to save time
warning off;
figure;patch('Vertices', X, 'Faces', C,...
    'FaceVertexCData',sf,'FaceColor', 'flat','EdgeColor','none','FaceAlpha', 1);
daspect([1 1 1]);axis off;%lighting phong;camlight
view(0,-90);
%title(field_name);
warning on;
% figure;hist(sf(:),100);
    Y1 = Y_LK_pf;
    Y2 = Y_LK_sf;
%%
function [X, C] = mesh_gen_rand(dim)
P = partsphere(dim^2);
x = reshape(P(1,:),dim,dim);
y = reshape(P(2,:),dim,dim);
z = reshape(P(3,:),dim,dim);
X = double([x(:) y(:) z(:)]);
[C] = convhulln(X, {'Qt'});c = reshape(C,length(C)*3,1);c = unique(c);c = sort(c);
X = X(c,:);


