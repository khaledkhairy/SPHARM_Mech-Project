function precalc_Estrain(obj)
global diagx1 diagy1 diagz1 SPI diag_vec_1 diag_vec diag_vec_m1 diag_vec_2 diag_vec_m2
global diag_vec_ds diag_vec_1_ds diag_vec_m1_ds SPI2
global F_areas_o invDM DM 

X = double(obj.X);
C = obj.F;
nfaces = size(C,1);
u = X(:,1);v = X(:,2);w = X(:,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
crossqpr = cross([u(C(:,2))-u(C(:,1)) v(C(:,2))-v(C(:,1)) w(C(:,2))-w(C(:,1))],...
    [u(C(:,3))-u(C(:,1)) v(C(:,3))-v(C(:,1)) w(C(:,3))-w(C(:,1))],2);
twoA = sqrt(sum(crossqpr.*crossqpr,2));  
F_areas_o = twoA/2;
n = crossqpr./twoA(:,ones(3,1));

%%%%% build the necessary transformations as in Computer Graphics p. 217.
%%% Get the vertices 
V1 = [u(C(:,1)) v(C(:,1)) w(C(:,1))];V2 = [u(C(:,2)) v(C(:,2)) w(C(:,2))];V3 = [u(C(:,3)) v(C(:,3)) w(C(:,3))];
%%% construct the translation matrix that brings V1 to the origin
T  = speye(4*nfaces);
diagx1 = zeros(1,size(T,1));diagx1(4:4:end) = -V1(:,1);T = T + spdiags(diagx1',3,size(T,1), size(T,2));
diagy1 = zeros(1,size(T,1));diagy1(4:4:end) = -V1(:,2);T = T + spdiags(diagy1',2,size(T,1), size(T,2));
diagz1 = zeros(1,size(T,1));diagz1(4:4:end) = -V1(:,3);T = T + spdiags(diagz1',1,size(T,1), size(T,2));
%%% construct the necessary rotation matrix
R = speye(4*nfaces); p1p2 = V2-V1;p1p3 = V3-V1;
normp1p2 = sqrt(p1p2(:,1).^2 + p1p2(:,2).^2 + p1p2(:,3).^2);
rz = p1p2./normp1p2(:,ones(3,1));
crossp1p3p1p2 = cross(p1p3,p1p2, 2);normcrossp1p3p1p2 = sqrt(crossp1p3p1p2(:,1).^2 + crossp1p3p1p2(:,2).^2 + crossp1p3p1p2(:,3).^2);
rx = crossp1p3p1p2./normcrossp1p3p1p2(:,ones(3,1));ry = cross(rz, rx, 2);
diag_vec = ones(1,4 * nfaces);diag_vec(1:4:end) =rx(:,1);diag_vec(2:4:end) = ry(:,2);diag_vec(3:4:end) = rz(:,3);
R = R + spdiags(diag_vec',0,4*nfaces, 4*nfaces);
diag_vec_1 = zeros(1,4 * nfaces);diag_vec_1(2:4:end) = rx(:,2);diag_vec_1(3:4:end) = ry(:,3);R = R + spdiags(diag_vec_1',1,size(R,1), size(R,2));
diag_vec_m1 = zeros(1,4 * nfaces);diag_vec_m1(1:4:end) = ry(:,1);diag_vec_m1(2:4:end) = rz(:,2);R = R + spdiags(diag_vec_m1',-1,size(R,1), size(R,2));
diag_vec_2 = zeros(1,4 * nfaces);diag_vec_2(3:4:end) = rx(:,3);R = R + spdiags(diag_vec_2',2,size(R,1), size(R,2));
diag_vec_m2 = zeros(1,4 * nfaces);diag_vec_m2(1:4:end) = rz(:,1);R = R + spdiags(diag_vec_m2',-2,size(R,1), size(R,2));
R = R - speye(4 * nfaces);
V1pr = reshape([V1 ones(length(V1),1)]',length(V1)*4,1);
V2pr = reshape([V2 ones(length(V2),1)]',length(V2)*4,1);
V3pr = reshape([V3 ones(length(V3),1)]',length(V3)*4,1);
V1r = R*T*V1pr; V2r = R * T *V2pr; V3r = R * T * V3pr; % transform
V1r = [reshape(V1r,4,length(V1))]';V2r = [reshape(V2r,4,length(V2))]';V3r = [reshape(V3r,4,length(V3))]';
%%% All vertices are rotated now, and all x-coordinates are zero.now we omit the x-coordinate and look at the triangles in 2-D only. note that all V1r are translated to zero.
V2r = [V2r(:,2) V2r(:,3)]; V3r = [V3r(:,2) V3r(:,3)];
%%% use the rotated vertices to build the necessary matrices for the calculation of the deformation
DM = speye(2*nfaces);%invDM = sparse(nfaces*2);
diag_vec = zeros(1,2*nfaces);diag_vec(1:2:end) =V2r(:,1);diag_vec(2:2:end) = V3r(:,2);DM = DM + spdiags(diag_vec',0,size(DM,1), size(DM,2));
diag_vec_1 = zeros(1,2*nfaces);diag_vec_1(2:2:end)=V3r(:,1);DM = DM + spdiags(diag_vec_1',1,size(DM,1), size(DM,2));
diag_vec_m1 = zeros(1,2*nfaces);diag_vec_m1(1:2:end)=V2r(:,2);DM = DM + spdiags(diag_vec_m1',-1,size(DM,1), size(DM,2));
DM = DM-speye(size(DM));

% invDM = inv(DM); 
%%% manually invert DM to make the inversion possible for fine meshes
indDM = 1:2:2*nfaces;
invDM = sparse(size(DM));
for ix = 1:nfaces,      % loop over the triangles (direction cosines rotation, and prep. DM)
    Dm = DM(indDM(ix):indDM(ix)+1,indDM(ix):indDM(ix)+1);
    invDM(indDM(ix):indDM(ix)+1,indDM(ix):indDM(ix)+1) = inv(Dm);
end


clear T R normp1p2 rz crossp1p3p1p2
clear  V1pr V2pr V3pr V1r V2r V3r
%%

diagx1 = zeros(1,4*nfaces);diagy1 = zeros(1,4*nfaces);diagz1 = zeros(1,4*nfaces);
SPI = speye(4*nfaces); SPI2 = speye(2*nfaces);
diag_vec_1 = zeros(1,4 * nfaces);diag_vec_m1 = zeros(1,4 * nfaces);diag_vec_2 = zeros(1,4 * nfaces);diag_vec_m2 = zeros(1,4 * nfaces);
diag_vec_ds = zeros(1,2*nfaces);diag_vec_1_ds = zeros(1,2*nfaces);diag_vec_m1_ds = zeros(1,2*nfaces);
end