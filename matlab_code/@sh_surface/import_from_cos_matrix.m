function [s cos_mx_basis] = import_from_cos_matrix(s,mx, cos_mx_basis)
%%   generate the field coefficients from an imported  matrix
% the matrix was a greyscale image (e.g. cos_view_sn_coarse.tif)
% (maybe modified by some photoeditor) and imported by imtool
%incr = 0.1;
persistent Y_LK_sf U_sf S_sf V_sf invS_sf
if nargin == 3, 
    Y_LK_sf = cos_mx_basis.Y_LK_sf;
    U_sf = cos_mx_basis.U_sf;
    S_sf = cos_mx_basis.S_sf;
    V_sf = cos_mx_basis.V_sf;
    invS_sf = cos_mx_basis.invS_sf;
end
incrt = pi/size(mx,2) ;
incrp = 2*pi/size(mx,1);
t = 0:incrt:pi;
p = 0:incrp:2*pi;

t = t(1:size(mx,2));
p = p(1:size(mx,1));
[t p] = meshgrid(t,p);


%%%% 
[L, K] = shp_surface.indices_gen(1:(s.L_max + 1)^2); M = length(L);N = numel(t);

if size(Y_LK_sf,1)~= N || size(Y_LK_sf,2)~= M
%     disp([N M]);
%     disp(size(Y_LK_sf));
%     disp('Generating new basis ...');
    Y_LK_sf  = zeros(N, M, 'single');
    for S = 1:length(L),
        %disp([num2str(S) ' of ' num2str(length(L))]);
        Y_LK_sf(:,S) = sh_basis.ylk_bosh(L(S),K(S),p(:)',t(:)')'; % uses bosh version
    end
    [U_sf, S_sf, V_sf] = svd(Y_LK_sf, 'econ');
    invS_sf = 1./(S_sf);invS_sf(invS_sf==inf) = 0;
end

[s.xc] = (V_sf*invS_sf) * (U_sf'*double(mx(:)));
      
    cos_mx_basis.Y_LK_sf = Y_LK_sf;
    cos_mx_basis.U_sf = U_sf;
    cos_mx_basis.S_sf = S_sf;
    cos_mx_basis.V_sf = V_sf;
    cos_mx_basis.invS_sf = invS_sf;
    
% sf = Y_LK_sf(:,1:length(s.xc))*s.xc(:);%% generate the  scalar field
% tgrid = t; %reshape(t,size(ct));
% pgrid = p;%reshape((p-pi), size(ct));
% 
% surf((tgrid), (pgrid),double(reshape(sf, size(t))), 'EdgeColor', 'none');
% xlim([0 1]);ylim([-1 1]);axis tight;
% view(0,90);
% sf = reshape(sf,size(t));