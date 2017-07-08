function [E_stretch, E_shear] = calc_Estrain(obj)
global F_areas_o invDM DM 
global diagx1 diagy1 diagz1 SPI diag_vec_1 diag_vec diag_vec_m1 diag_vec_2 diag_vec_m2
global diag_vec_ds diag_vec_1_ds diag_vec_m1_ds SPI2

X = double(obj.X);
C = obj.F;
nfaces = size(C,1);
u = X(:,1);v = X(:,2);w = X(:,3);


%% calculation of the area and volume
u(:) = X(:,1); v(:) = X(:,2); w(:) = X(:,3);
crossqpr = cross([u(C(:,2))-u(C(:,1)) v(C(:,2))-v(C(:,1)) w(C(:,2))-w(C(:,1))],[u(C(:,3))-u(C(:,1)) v(C(:,3))-v(C(:,1)) w(C(:,3))-w(C(:,1))],2);
twoA = sqrt(sum(crossqpr.*crossqpr,2)); A = sum(twoA)/2; F_areas = twoA/2;n = crossqpr./twoA(:,ones(3,1));
V = -sum(1/3*(dot(n,[u(C(:,1)) v(C(:,1)) w(C(:,1))], 2).*twoA./2));

%% calculate E_stretch (accoring to Lim)
a3 = -2; a4 = 8;k_stretch = 1;
ai = (F_areas./F_areas_o -1);     % area_invariant
E_stretch = k_stretch/2 * sum((ai.*ai + a3 * ai .*ai.*ai +a4 * ai.*ai.*ai.*ai).*F_areas_o) ;

%% MS shear energy



%%%%% build the necessary transformations as in Computer Graphics p. 217.
V1 = [u(C(:,1)) v(C(:,1)) w(C(:,1))];V2 = [u(C(:,2)) v(C(:,2)) w(C(:,2))];V3 = [u(C(:,3)) v(C(:,3)) w(C(:,3))];
%%% construct the translation matrix that brings V1 to the origin
T  = speye(4*nfaces);
% diagx1 = zeros(1,size(T,1));diagy1 = zeros(1,size(T,1));diagz1 = zeros(1,size(T,1));
diagx1(:) = 0;diagy1(:) = 0;diagz1(:) = 0;
diagx1(4:4:end) = -V1(:,1);T = T + spdiags(diagx1',3,size(T,1), size(T,2));
diagy1(4:4:end) = -V1(:,2);T = T + spdiags(diagy1',2,size(T,1), size(T,2));
diagz1(4:4:end) = -V1(:,3);T = T + spdiags(diagz1',1,size(T,1), size(T,2));
%%% construct the necessary rotation matrix
R = speye(4*nfaces); p1p2 = V2-V1;p1p3 = V3-V1;
normp1p2 = sqrt(p1p2(:,1).^2 + p1p2(:,2).^2 + p1p2(:,3).^2);
rz = p1p2./normp1p2(:,ones(3,1));
crossp1p3p1p2 = cross(p1p3,p1p2, 2);normcrossp1p3p1p2 = sqrt(crossp1p3p1p2(:,1).^2 + crossp1p3p1p2(:,2).^2 + crossp1p3p1p2(:,3).^2);
rx = crossp1p3p1p2./normcrossp1p3p1p2(:,ones(3,1));ry = cross(rz, rx, 2);
diag_vec = ones(1,4 * nfaces);diag_vec(1:4:end) =rx(:,1);diag_vec(2:4:end) = ry(:,2);diag_vec(3:4:end) = rz(:,3);
R = R + spdiags(diag_vec',0,4*nfaces, 4*nfaces);
% diag_vec_1 = zeros(1,4 * nfaces);diag_vec_m1 = zeros(1,4 * nfaces);diag_vec_2 = zeros(1,4 * nfaces);diag_vec_m2 = zeros(1,4 * nfaces);
diag_vec_1(:) = 0;diag_vec_m1(:) = 0;diag_vec_2(:) = 0;diag_vec_m2(:) = 0;
diag_vec_1(2:4:end) = rx(:,2);diag_vec_1(3:4:end) = ry(:,3);R = R + spdiags(diag_vec_1',1,size(R,1), size(R,2));
diag_vec_m1(1:4:end) = ry(:,1);diag_vec_m1(2:4:end) = rz(:,2);R = R + spdiags(diag_vec_m1',-1,size(R,1), size(R,2));
diag_vec_2(3:4:end) = rx(:,3);R = R + spdiags(diag_vec_2',2,size(R,1), size(R,2));
diag_vec_m2(1:4:end) = rz(:,1);R = R + spdiags(diag_vec_m2',-2,size(R,1), size(R,2));
R = R - SPI; %speye(4 * nfaces);
V1pr = reshape([V1 ones(length(V1),1)]',length(V1)*4,1);
V2pr = reshape([V2 ones(length(V2),1)]',length(V2)*4,1);
V3pr = reshape([V3 ones(length(V3),1)]',length(V3)*4,1);
V1r = R*T*V1pr; V2r = R * T *V2pr; V3r = R * T * V3pr; % transform
V1r = [reshape(V1r,4,length(V1))]';V2r = [reshape(V2r,4,length(V2))]';V3r = [reshape(V3r,4,length(V3))]';
%%% All vertices are rotated now, and all x-coordinates are zero.now we omit the x-coordinate and look at the triangles in 2-D only. note that all V1r are translated to zero.
V2r = [V2r(:,2) V2r(:,3)]; V3r = [V3r(:,2) V3r(:,3)];
%%% use the rotated vertices to build the necessary matrices for the calculation of the deformation
DS = SPI2;%speye(2*nfaces);
% diag_vec = zeros(1,2*nfaces);diag_vec_1 = zeros(1,2*nfaces);diag_vec_m1 = zeros(1,2*nfaces);
diag_vec_ds(:) = 0;diag_vec_1_ds(:) = 0;diag_vec_m1_ds(:) = 0;
diag_vec_ds(1:2:end) =V2r(:,1);diag_vec_ds(2:2:end) = V3r(:,2);DS = DS + spdiags(diag_vec_ds',0,2*nfaces, 2*nfaces);
diag_vec_1_ds(2:2:end)=V3r(:,1);DS = DS + spdiags(diag_vec_1_ds',1,2*nfaces, 2*nfaces);
diag_vec_m1_ds(1:2:end)=V2r(:,2);DS = DS + spdiags(diag_vec_m1_ds',-1,2*nfaces, 2*nfaces);
DS = DS - SPI2;%speye(size(DS));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G = 1/2 * (invDM' * DS' * DS * invDM - SPI2);;
diagG = diag(G);eps11 = diagG(1:2:end);eps22 = diagG(2:2:end);diag1G = diag(G,1);eps12 = diag1G(1:2:end);
b = -(eps11+eps22); c = eps11.*eps22 -eps12.*eps12; % a is always 1
E1 = (-b + sqrt(b.*b - 4.* c))./2;E2 = (-b - sqrt(b.*b - 4.* c))./2;%disp([E1 E2]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta  = (E1 + E2 - ai)./(1 + ai);b1 = 0.7; b2 = 0.75;k_shear = 1;
E_shear = k_shear *sum((beta + b1 * ai .*beta + b2 * beta.*beta).*F_areas_o) ;
% E = E + shear_fac * E_shear;