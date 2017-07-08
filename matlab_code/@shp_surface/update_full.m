function obj = update_full(obj)
if obj.needs_updating
    lb = (obj.L_max + 1)*(obj.L_max + 1);
    gdimp = length(obj.basis.p);gdimt = length(obj.basis.t);
    c = obj.xc(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
    obj.x = sum(c.*(obj.basis.Y(:,:,1:lb)),3);
    xp =(sum(c.*obj.basis.Y_P(:,:,1:lb),3));
    xt = (sum(c.*obj.basis.Y_T(:,:,1:lb),3));
    xpp = (sum(c.*obj.basis.Y_PP(:,:,1:lb),3));
    xtt = (sum(c.*obj.basis.Y_TT(:,:,1:lb),3));
    xtp = (sum(c.*obj.basis.Y_TP(:,:,1:lb),3));
    
    c = obj.yc(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
    obj.y = sum(c.*obj.basis.Y(:,:,1:lb),3);
    yp =(sum(c.*obj.basis.Y_P(:,:,1:lb),3));
    yt = (sum(c.*obj.basis.Y_T(:,:,1:lb),3));
    ypp = (sum(c.*obj.basis.Y_PP(:,:,1:lb),3));
    ytt = (sum(c.*obj.basis.Y_TT(:,:,1:lb),3));
    ytp = (sum(c.*obj.basis.Y_TP(:,:,1:lb),3));
    
    c = obj.zc(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
    obj.z = sum(c.*obj.basis.Y(:,:,1:lb),3);
    zp =(sum(c.*obj.basis.Y_P(:,:,1:lb),3));
    zt = (sum(c.*obj.basis.Y_T(:,:,1:lb),3));
    zpp = (sum(c.*obj.basis.Y_PP(:,:,1:lb),3));
    ztt = (sum(c.*obj.basis.Y_TT(:,:,1:lb),3));
    ztp = (sum(c.*obj.basis.Y_TP(:,:,1:lb),3));
    
    %%%%%%%%%%%%%%%%%%%%%%%% Calculate first and second fundamental forms
    X  =[obj.x(:) obj.y(:) obj.z(:)];
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
    obj.V = abs(1/3.*sum(obj.basis.w.*(dot(X,n,2)).*SSn));
    obj.A = sum(obj.basis.w.*SSn);
    obj.H = -(E.*N + G.*L - 2*F.*M)./(2 * (E.*G-F.*F));
    obj.KG = (L.*N - M.*M)./(E.*G-F.*F);
    %% other interesting quantities that we get practically for free
    obj.h = 1./obj.A.*sum(obj.H(:).*obj.basis.w.*SSn);                % total mean curvature
    obj.T = sum(sum(obj.KG(:).*obj.basis.w.*SSn))/4/pi;            % total Gaussian curvature (constant for topology)
    obj.S = (2.*obj.H.^2-obj.KG).^(0.5);                     % curvedness
    obj.Eb = 1/2*sum((2.*obj.H(:)).^2.*obj.basis.w.*SSn)./8/pi;       % bending energy relative to that of the sphere
    %%% calculating reduced volume
    r = sqrt(obj.A/4/pi);% Area of a sphere = 4 * pi * r^2
    V_sphere = 4/3 * pi *r^3;% Volume of that sphere = 4/3 * pi * r^3
    obj.v = obj.V/V_sphere;
    
    obj.needs_updating = 0;
end