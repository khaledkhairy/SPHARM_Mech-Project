function plot_difference(m1,m2)
[A1, V1, v1, F_areas1, h1, H1, wb1, da1, N1, K1, k_g1, dA1, n1] = ...
    triangulated_props(m1.X, m1.F, 0);
[A2, V2, v2, F_areas2, h2, H2, wb2, da2, N2, K2, k_g2, dA2, n2] = ...
    triangulated_props(m2.X, m2.F, 0);

s = 1;
V = [s*(m1.X(:,1)-m2.X(:,1)) s*(m1.X(:,2)-m2.X(:,2)) s*( m1.X(:,3)-m2.X(:,3))];
dotVN = dot(V,N1,2);
dotVN = dotVN(:,ones(1,size(N1,2)));
D = V - dotVN.*N1;    % direction perpendicular to the surface normal



M = sqrt((m1.X(:,1)-m2.X(:,1)).^2 + (m1.X(:,2)-m2.X(:,2)).^2 + (m1.X(:,3)-m2.X(:,3)).^2); % magnitude
sf_diff = {};
sf_diff{1} = 'distance with C_o = egg';
sf_diff{2} = M;
m1.sf{12} = sf_diff;

%figure;plot_field(m1,12);
patch('vertices',m1.X,'faces',m1.F,'FaceColor','green','edgecolor', 'none');
hold on;
h = quiver3(m1.X(:,1), m1.X(:,2),...
    m1.X(:,3),D(:,1),D(:,2),D(:,3),'LineWidth',1.5);
%adjust_quiver_arrowhead_size(h, 2);
axis off;axis equal; camlight;
lighting phong;
view(-90,90);camlight;
view(-90,-45);camlight;
view(-90,-72);camlight;
view(158,-30);camlight;
view(90,0);
set(gcf,'Color','white');
cameramenu