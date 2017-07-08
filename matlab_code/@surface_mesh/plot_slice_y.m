function [X F X_in X_out P] = plot_slice_y(obj, d, yval, flag)
%% plot the sliced object
if nargin ==1, d = 1;yval = 0;end
if nargin<4, flag  = 1;end
X = obj.X;
F = obj.F;
centers = [];
cflag = 1;
Fs = F;
if cflag == 1,
    colordef black;
    %h = figure('Color', [0 0 0]);
    set(gcf,'InvertHardcopy','off');
else
    colordef white;
    %h = figure('Color', [1 1 1]);
    set(gcf,'InvertHardcopy','on');
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
if sign(flag) == -1,
    [A, V, v, F_areas, h, H, wb, da, N, K, k_g, dA] = triangulated_props([X(:,1) -X(:,2) X(:,3)], F, 0);
else
    [A, V, v, F_areas, h, H, wb, da, N, K, k_g, dA] = triangulated_props(X, F, 0);
end
Xs = X - N .* d/2;
% %% find the triangles that intersect that surface
% large = 1e3;
% X2 = [-large yval -large; 0 yval large; large yval -large];F2 = [1 2 3];
% [int_indx, TP] = tri_plane_intersect(X,F, X2, F2, [], 0); % determine the triangle indices that intersect the plane
% % patch('vertices',X,'faces',C,'FaceColor','red', 'FaceAlpha', 0.8);  % the middle surface

%% find the surface after cutting
[X F] = get_cut_surface(X, F, yval);
[E,L,face_memb] = edge_info(X,F);

if sign(flag) == -1, X(:,2) = -X(:,2);end
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
X_in = X + N .* d/2;  %%% now generate the inner surface
X_out = X - N .* d/2;  %%% now generate the outer surface

% X_in = X;
% X_out = X-N.*d;

X_in(c_ix,2) = yval;
X_out(c_ix,2) = yval;
%% make the connections to get a cut through the shell thickness
%[X F] = connect_inner_and_outer_shells(X_in, X_out,F, c_ix)

%% define the cut-mesh surface (using the contour indices) and plot it 
view(-180,0);
P1 = [];
P2 = [];
P3 = [];
for ix = 1:length(c_ix),  % loop through the connected pairs
    p1 = X_in(c_ix(ix),:);
    p2 = X_out(c_ix(ix),:);
    %plot3([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)],'w-');hold on;
    linx = L{c_ix(ix)};
    for lix = 1:length(linx)
        v1 = X_in(linx(lix),:);
        v2 = X_out(linx(lix),:);
        if v1(2)==yval,  p3 = v1; p4  = v2;
            %p = [p1;p2;p3;p4];patch(p(:,1), p(:,2), p(:,3),'w','EdgeColor','none');
            p = [p1;p3;p4;p2];
            patch(p(:,1), p(:,2), p(:,3),'w','EdgeColor','b', 'FaceAlpha',1.0);
            P1 = [P1 p(:,1)];
            P2 = [P2 p(:,2)];
            P3 = [P3 p(:,3)];
        end
    end
    
end
P.P1 = P1;
P.P2 = P2;
P.P3 = P3;
%%
% if cflag == 1,
    if abs(flag)==1
    patch('vertices',X_out,'faces',F,'facecolor','green','edgecolor','none', 'FaceAlpha', 1);
    patch('vertices',X_in,'faces',F,'facecolor','red','edgecolor',[0.2 0.2 0.2], 'FaceAlpha', 1);
    elseif abs(flag) ==2
        patch('vertices',X_out,'faces',F,'facecolor','green','edgecolor','none', 'FaceAlpha', 0.9);
        patch('vertices',X_in,'faces',F,'facecolor','red','edgecolor',[0.2 0.2 0.2], 'FaceAlpha', 0.9);
    end
    %patch('vertices',Xs,'faces',Fs,'facecolor','green','edgecolor','none','FaceAlpha',0.5);
% else
%     patch('vertices',X_in,'faces',F,'facecolor','green','edgecolor','none', 'FaceAlpha', 1);
%     patch('vertices',X_out,'faces',F,'facecolor','red','edgecolor','none', 'FaceAlpha', 1);
% end

% % axis tight;axis off;axis equal;
% % view(-155,0);camlight;
% % view(136,-2);
% % lighting phong;camlight
% % cameramenu
%plot_mesh(X_in,F);hold on;plot_mesh(X_out, F);
%%










