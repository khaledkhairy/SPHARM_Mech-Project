function plot_H_2d(obj)
%% 
% dfig;
incr = 0.01;
t = 0:incr/2:pi;
p = -pi/2:incr:3*pi/2;
[t p] = meshgrid(t,p);
L_max = obj.basis.L_max;
%%%%
[Y P]                       = obj.basis.ylk_cos_sin_bosh(p, t, L_max);
[Y_P]                       = obj.basis.ylk_cos_sin_dphi_bosh(p, t, L_max, P);
[Y_T P_T]                   = obj.basis.ylk_cos_sin_dtheta_bosh(p, t, L_max, P);
Y_PP                    	= obj.basis.ylk_cos_sin_dphiphi_bosh(p, t, L_max, P);
Y_TT                        = obj.basis.ylk_cos_sin_dthetatheta_bosh(p, t, L_max,P_T);
Y_TP                        = obj.basis.ylk_cos_sin_dthetaphi_bosh(p, t, L_max, P_T);


%%%
lb = (L_max + 1)*(L_max + 1);
gdimp = length(p);gdimt = length(t);
c = obj.xc(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
x = sum(c.*(Y(:,:,1:lb)),3);
xp =(sum(c.*Y_P(:,:,1:lb),3));
xt = (sum(c.*Y_T(:,:,1:lb),3));
xpp = (sum(c.*Y_PP(:,:,1:lb),3));
xtt = (sum(c.*Y_TT(:,:,1:lb),3));
xtp = (sum(c.*Y_TP(:,:,1:lb),3));

c = obj.yc(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
y = sum(c.*Y(:,:,1:lb),3);
yp =(sum(c.*Y_P(:,:,1:lb),3));
yt = (sum(c.*Y_T(:,:,1:lb),3));
ypp = (sum(c.*Y_PP(:,:,1:lb),3));
ytt = (sum(c.*Y_TT(:,:,1:lb),3));
ytp = (sum(c.*Y_TP(:,:,1:lb),3));

c = obj.zc(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
z = sum(c.*Y(:,:,1:lb),3);
zp =(sum(c.*Y_P(:,:,1:lb),3));
zt = (sum(c.*Y_T(:,:,1:lb),3));
zpp = (sum(c.*Y_PP(:,:,1:lb),3));
ztt = (sum(c.*Y_TT(:,:,1:lb),3));
ztp = (sum(c.*Y_TP(:,:,1:lb),3));

%%%%%%%%%%%%%%%%%%%%%%%% Calculate first and second fundamental forms
X  =[x(:) y(:) z(:)];
Xt = [xt(:) yt(:) zt(:)];
Xp = [xp(:) yp(:) zp(:)];
Xpp = [xpp(:) ypp(:) zpp(:)];
Xtp = [xtp(:) ytp(:) ztp(:)];
Xtt = [xtt(:) ytt(:) ztt(:)];

E = dot(Xt,Xt,2);
F = dot(Xt,Xp,2);
G = dot(Xp,Xp,2);
SS = (cross(Xt,Xp,2));SSn = sqrt(E.*G-F.*F);
n = SS./SSn(:,ones(1,3));
L = dot(Xtt,n,2);
M = dot(Xtp,n,2);
N = dot(Xpp,n,2);


%%%%%%%%%%%%%%%%%%%%%%%% update geometrical properties
H = -(E.*N + G.*L - 2*F.*M)./(2 * (E.*G-F.*F));

%% %
tgrid = t; %reshape(t,size(ct));
pgrid = p;%reshape((p-pi), size(ct));

surf((tgrid), (pgrid),double(reshape(H, size(t))), 'EdgeColor', 'none');
xlim([0 1]);ylim([-1 1]);axis tight;
view(0,90);


