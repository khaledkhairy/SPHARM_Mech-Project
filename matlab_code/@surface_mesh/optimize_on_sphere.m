function obj = optimize_on_sphere(obj)
% optimize the t,p coordinates for mapping a surface mesh to the sphere
% We use a large scale algorithm in matlab's optimization toolbox
global itercount
% define initial quantities
X = obj.X;
F = obj.F;
t = obj.t;
p = obj.p;
Xo = [t(:);p(:)];
nvert = length(obj.X);
% determine sparsity pattern for jacobian in order to be able to use the large-scale algorithm
disp('Calculating Jacobian Sparsity Pattern');
JacPat = sparse(length(F),nvert);
for ix = 1:length(F),
    verts = F(ix,:);
    for vert = 1:length(verts),
        JacPat(ix,verts(vert)) = 1;
    end
end
JacPat = sparse([JacPat; JacPat]);
%    JacPat = sparse([JacPat;zeros(4,size(JacPat,2))]);

%% % Calculate the relative triangle areas of the original object
u = X(:,1); v = X(:,2); w = X(:,3); % the vertex coordinates
x1 = u(F(:,1)); y1 = v(F(:,1));z1 =  w(F(:,1));
x2 = u(F(:,2)); y2 = v(F(:,2));z2 =  w(F(:,2));
x3 = u(F(:,3)); y3 = v(F(:,3));z3 =  w(F(:,3));%%% generate the three vertices (coordinates) of the triangles in same order as F
q = [x2-x1 y2-y1 z2-z1]; r = [x3-x1 y3-y1 z3-z1];
crossqpr = cross(q,r,2);
twoA = sqrt(sum(crossqpr.^2,2));   % take the norm
A = sum(twoA)/2;                    % this is the total area
F_areas = twoA/2;                   % this is the vector of face areas
Areas_o_bar = F_areas/A;            % vector of relative areas for each triangle of original object
%% Calculate the angles A, B, C at each vertex of each face of the original object
V1 = [x1(:) y1(:) z1(:)];   % list of coordinates for first vertices of all triangles F
V2 = [x2(:) y2(:) z2(:)];
V3 = [x3(:) y3(:) z3(:)];
nV1 = norm_list(V1);
nV2 = norm_list(V2);
nV3 = norm_list(V3);

a = acos(dot(V1,V3,2)./nV1./nV3);  % length (in angle) of arc from one vertex to the other
b = acos(dot(V1,V2,2)./nV1./nV2);  
c = acos(dot(V2,V3,2)./nV2./nV3);  
s = (a + b + c)/2;

% Calculate the angles at the vertices
A_o = 2 * atan(sqrt(sin(s-b).*sin(s-c)./sin(s)./sin(s-a)));% angle between V1 and V3 (for all faces)
B_o = 2 * atan(sqrt(sin(s-a).*sin(s-c)./sin(s)./sin(s-b)));% angle between V1 and V2 (for all faces)
C_o = 2 * atan(sqrt(sin(s-a).*sin(s-b)./sin(s)./sin(s-c)));% angle between V2 and V3 (for all faces)

%% configure optimization and execute
plot_flag = 0;
itercount = 1;
options = optimset('MaxFunEvals', 400000,'DiffMaxChange', 1e0,'DiffMinChange', 1e-8,...
    'DerivativeCheck','off','GradObj','off','GradConstr','off',...
    'MaxIter', 20,'Display', 'iter','Diagnostics', 'on','LevenbergMarquardt','on',...
    'TolCon', 1e-6,'TolFun', 1e-16,'TolX', 1e-16, 'MaxSQPIter',30,...
    'LargeScale','on', 'JacobPattern',JacPat,'LineSearchType','quadcubic',...
    'PrecondBandWidth',0,'TypicalX',0 * ones(length(Xo),1),...
    'MaxPCGIter',10);
[Xf] = lsqnonlin(@optimize_on_sphere_objective,Xo,[],[],options,F,Areas_o_bar, A_o, B_o, C_o, plot_flag);
obj.t = Xf(1:nvert);
obj.p = Xf(nvert+1:end);

plot_flag = 1;
res = optimize_on_sphere_objective(Xf, F, Areas_o_bar, A_o, B_o, C_o, plot_flag);


function res = optimize_on_sphere_objective(Xo, F, Areas_o_bar, A_o, B_o, C_o, plot_flag)
% the objective function that returns the sum of squared residuals of
% differences in relative area between the triangles on the original object
% and the spherical ones on the unit sphere
global itercount
nvert = length(Xo)/2;
t = Xo(1:nvert);
p = Xo(nvert+1:end);
[u, v, w] = kk_sph2cart(t,p,1);
%% Calculate the Euler (geodesic) area of the spherical triangles
%%% Area of spherical triangle = R^2*(A+B+C-pi);
%%% where A, B and C are the angles in radians
x1 = u(F(:,1)); y1 = v(F(:,1));z1 =  w(F(:,1));
x2 = u(F(:,2)); y2 = v(F(:,2));z2 =  w(F(:,2));
x3 = u(F(:,3)); y3 = v(F(:,3));z3 =  w(F(:,3));%%% generate the three vertices (coordinates) of the triangles in same order as F

%%% for testing
% q = [x2-x1 y2-y1 z2-z1]; r = [x3-x1 y3-y1 z3-z1];
% crossqpr = cross(q,r,2);
% twoA = sqrt(sum(crossqpr.^2,2));   % take the norm
% Ap = sum(twoA)/2;                    % this is the total area
%dfig(4);plot_state(u,v,w, F, 1, 1, 0, 0);%% to look at some triangles on the sphere
%%% END TEST

V1 = [x1(:) y1(:) z1(:)];   % list of coordinates for first vertices of all triangles F
V2 = [x2(:) y2(:) z2(:)];
V3 = [x3(:) y3(:) z3(:)];

% nV1 = norm_list(V1);
% nV2 = norm_list(V2);
% nV3 = norm_list(V3);

% calculate the three angles that constitute the arcs between the vertices
% (expressed as angles)
% % % % ca = (dot(V1,V3,2));  % angle between V1 and V3 (for all faces)
% % % % cb = (dot(V1,V2,2));  % angle between V1 and V2 (for all faces)
% % % % cc = (dot(V2,V3,2));  % angle between V2 and V3 (for all faces)

a = acos(dot(V1,V3,2));  % angle between V1 and V3 (for all faces)
b = acos(dot(V1,V2,2));  % angle between V1 and V2 (for all faces)
c = acos(dot(V2,V3,2));  % angle between V2 and V3 (for all faces)
s = (a + b + c)/2;

% Calculate the angles at the vertices
A = 2 * atan(sqrt(sin(s-b).*sin(s-c)./sin(s)./sin(s-a)));
B = 2 * atan(sqrt(sin(s-a).*sin(s-c)./sin(s)./sin(s-b)));
C = 2 * atan(sqrt(sin(s-a).*sin(s-b)./sin(s)./sin(s-c)));

% the areas are given in terms of the Excess as:
% Area of spherical triangle = R^2 * [A + B + C -pi] = R^2 * E
% A, B and C must be in radians
% since we have a unit sphere R = 1 the area is simply E

Excess = A+B+C-pi;
Areas_bar = Excess/sum(Excess); % normalize the areas so that their sum is equal to one (i.e. relative area)

% calculate the residual

 res1 = abs(Areas_bar-Areas_o_bar);
% res1 = Areas_bar;
% res1 = s;

% res2 = (abs(A-A_o) + abs(B-B_o) + abs(C-C_o));
res2 = (abs(A-pi/2) + abs(B-pi/2) + abs(C-pi/2));

gamma = 0;%mean(res1)/mean(res2);
res =  res1 + gamma *res2;


itercount = itercount + 1;
if plot_flag || mod(itercount,1000)==0
    % plot the triangles with residual coloring
    dfig(4);cla
    patch('Vertices',[u v w],'Faces',F, 'FaceVertexCData',res,'FaceColor','flat');
    graphlims = [-1.1 1.1];xlim(graphlims);ylim(graphlims); zlim(graphlims);view(3);axis square;
    colorbar;
    drawnow;
end







