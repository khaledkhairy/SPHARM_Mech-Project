function [lambda1 lambda2] = principal_stretches(obj,obj_o)
%% given a reference surface calculate the principal stretches
%%% yielding shear and stretch measures

% % %% calculate coefficients of the 1st fundamental form (components of g_lij)
% % E = kk_dot(Xt,Xt);          % g_l11
% % F = kk_dot(Xt,Xp);          % g_l12 and g_l21
% % G = kk_dot(Xp,Xp);          % g_l22
% % %% calculate differential area element and surface normal
% % SS = (kk_cross(Xt,Xp));
% % SSn = sqrt(E.*G-F.*F);
% % n = SS./SSn(:,ones(1,3));       % m_l3
% % %% set up covariant basis (with letter l), and contravariant basis with the letter u
% % m_l1 = Xt;
% % m_l2 = Xp;
% % bb = kk_norm(kk_cross(m_l1,m_l2));                  % bb is essentially SSn
% % bb = bb(:,ones(1,3));
% % m_l3 = kk_cross(Xt,Xp)./bb;
% % 
% % % bb = kk_dot(m_l1, kk_cross(m_l2,m_l3));bb = bb(:);  % formulation as in AMoS p.658
% % m_u1 = kk_cross(Xp,m_l3)./bb;
% % m_u2 = kk_cross(m_l3,Xt)./bb;
% % m_u3 = n;
% % %% generate unit tangent vectors -- for treatment as in AoPaS
% % Xtnorm = sqrt(kk_dot(Xt,Xt));
% % Xpnorm = sqrt(kk_dot(Xp,Xp));
% % 
% % tt = Xt./Xtnorm(:,ones(1,3));
% % tp = Xp./Xpnorm(:,ones(1,3));
% % tn = kk_cross(Xt,Xp)./Xpnorm(:,ones(1,3))./Xtnorm(:,ones(1,3));
% % %% set up the metric tensor components
% % % covariant components
% % g_l11 = E;
% % g_l22 = G;
% % g_l12 = F;
% % g_l13 = kk_dot(m_l1, m_l3);
% % g_l23 = kk_dot(m_l2, m_l3);
% % g_l33 = kk_dot(m_l3, m_l3);
% % % contravariant metric tensor components g_uij
% % g_u11 = kk_dot(m_u1, m_u1);
% % g_u12 = kk_dot(m_u1, m_u2);
% % g_u22 = kk_dot(m_u2, m_u2);
% % % g_u33 = kk_dot(m_u3, m_u3);
% % % g_u13 = kk_dot(m_u1, m_u3);
% % % g_u23 = kk_dot(m_u2, m_u3);
% % %% calculate the coefficients of the second fundamental form
% % L = kk_dot(Xtt,n);
% % M = kk_dot(Xtp,n);
% % N = kk_dot(Xpp,n);
% % SO = [(G.*L-F.*M)./(E.*G-F.*F) (E.*M-F.*L)./(E.*G-F.*F) (G.*M-F.*N)./(E.*G-F.*F) (E.*N-F.*M)./(E.*G-F.*F)]; % the shape operator
% % % SO = [(G.*L-F.*M) (E.*M-F.*L) (G.*M-F.*N) (E.*N-F.*M)]; % the shape operator
% % %% calculate  geometric properties: Area, Volume, reduced volume, Curvatures and self-check total Gaussian curvature
% % V = abs(1/3.*sum(w.*(dot(X,n,2)).*SSn));
% % A = sum(w.*SSn);
% % Vo = 4/3*pi*(A/4/pi)^(3/2);
% % v_red = V/Vo;
% % H = -(E.*N + G.*L - 2*F.*M)./(2 * (E.*G-F.*F));
% % K = (L.*N - M.*M)./(E.*G-F.*F);
% % h = 1./A.*sum(H(:).*w.*SSn);
% % T = sum(sum(K(:).*w.*SSn))/4/pi;       % total curvature (constant for topology)

% % % %% WORKS!! Calculate midplane lagrange strain tensor as in http://www.scribd.com/doc/27473029/5/Strain-in-Orthogonal-Curvilinear-Coordinates
% % % %%%% tested to be equal to the expression from AMoS p667 (see below)
% % % lambda = zeros(length(Xt),2);
% % % zer = zeros(length(Xt),1);
% % % gamma11 = 1/2*(g_l11-phys.g_l11)./sqrt(phys.g_l11)./sqrt(phys.g_l11);
% % % gamma12 = 1/2*(g_l12-phys.g_l12)./sqrt(phys.g_l11)./sqrt(phys.g_l22);
% % % gamma21 = 1/2*(g_l12-phys.g_l12)./sqrt(phys.g_l11)./sqrt(phys.g_l22);
% % % gamma22 = 1/2*(g_l22-phys.g_l22)./sqrt(phys.g_l22)./sqrt(phys.g_l22);
% % % eps_vec = [gamma11 gamma12 zer gamma21 gamma22 zer zer zer zer];    % construct the strain tensor
% % % for ix = 1:length(Xt)
% % %     eps = reshape(eps_vec(ix,:), 3,3);% the lagrangian strain tensor 3x3 matrix
% % %     FTF = 2*eps+eye(3,3);             % Right Cauchy-Green tensor
% % %     eFTF= sort(eig(FTF));             % calculate eigenvalues of right Cauchy-Green tensor (square root of which provides the principal stretches AMoS p.29)
% % %     tol = 1e-4;indx = find(eFTF<(1+tol) & eFTF>(1-tol)); eFTF(indx(1)) = [];
% % %     lambda(ix,:) = [sqrt(eFTF(1)) sqrt(eFTF(2))];       % store the principal stretches
% % % end
% % % b4 = (lambda(:,1)-lambda(:,2)).^2./2./lambda(:,1)./lambda(:,2); % shear invariant
% % % ai = real(lambda(:,1).*lambda(:,2) - 1);     % area dilation invariant (as function of principal stretches
% % % E_stretch4 = phys.k_stretch/2 * ...
% % %     sum((ai.*ai + phys.a3 * ai .*ai.*ai +phys.a4 * ai.*ai.*ai.*ai).*phys.SSn_o.*w) ;%% multiplied by undeformed area
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %% area invariant DG
% % % ai = (SSn./phys.SSn_o - 1);     % area_invariant
% % % E_stretch1 = phys.k_stretch/2 * ...
% % %     sum((ai.*ai + phys.a3 * ai .*ai.*ai +phys.a4 * ai.*ai.*ai.*ai).*phys.SSn_o.*w) ;%% multiplied by undeformed area
% % % %% % calculate the shear and stretch energy based on triangles
% % % X = phys.Y_LK_self(:,1:length(xclks))* [xclks(:) yclks(:) zclks(:)];
% % % [E_shear4 E_stretch3 E1 E2 b_tri] = shear_calc(X,phys.C_slave,phys.pre,phys.a3, phys.a4,phys.b1, phys.b2, phys.k_shear);
% % % %% WORKS!! calculate Strain tensor components from first term of F in AMoS: p 665 
% % % %%%  using loops
% % % lambda = zeros(length(Xt),2);
% % % FTF_vec = zeros(length(Xt),9);
% % % beta = zeros(length(Xt),1);
% % % E1v = zeros(length(Xt),1);
% % % E2v = zeros(length(Xt),1);
% % % eigFTF = ones(length(Xt), 3);
% % % DGT_vec = zeros(length(Xt), 9);
% % % for ix = 1:length(Xt)
% % %     DGT = kk_kron(m_l1(ix,:), phys.m_u1(ix,:)) +...
% % %           kk_kron(m_l2(ix,:), phys.m_u2(ix,:));
% % %     DGT_vec(ix,:) = DGT(:)';
% % %     DGT = reshape(DGT,3,3);%     DGT = DGT(1:2,1:2);
% % %     FTF = DGT'*DGT;     % right Cauchy-Green tensor
% % %     FTF_vec(ix,:) = reshape(FTF,1,9);
% % %     eFTF= eig(FTF);eFTF = eFTF(:)';
% % %     tol = 1e-4;indx = find(eFTF(1,:)<(tol) & eFTF(1,:)>(-tol)); eFTF(indx(1)) = [];
% % %     lambda(ix,:) = [sqrt(eFTF(1)) sqrt(eFTF(2))];
% % %     eigFTF(ix,1:2) = eFTF(:)';
% % %     eps = 1/2*(FTF-eye(size(FTF)));     % lagrange strain tensor
% % %     Eig_vec = eig(eps, 'nobalance');
% % %     tol = 1e-4;indx = find(Eig_vec<(-0.5+tol) & Eig_vec>(-0.5-tol)); Eig_vec(indx(1)) = [];
% % %     E1 = Eig_vec(1);
% % %     E2 = Eig_vec(2);
% % %     beta(ix) = (E1 + E2 - ai(ix))./(1 + ai(ix));    % shear component
% % %     E1v(ix) = E1;
% % %     E2v(ix) = E2;
% % % end
% % % b1 = beta;
% % % b2 = (lambda(:,1)-lambda(:,2)).^2./2./lambda(:,1)./lambda(:,2); % expression for beta as function of principal stretches
% % % 
% % % % stretch 2
% % % ai = lambda(:,1).*lambda(:,2) - 1;     % area_invariant (as function of principal stretches
% % % E_stretch2 = phys.k_stretch/2 * ...
% % %     sum((ai.*ai + phys.a3 * ai .*ai.*ai +phys.a4 * ai.*ai.*ai.*ai).*phys.SSn_o.*w) ;%% multiplied by undeformed area
% % % %% WORKS!! first term of DGT: calculated Cauchy-Green deformation tensor -> principal stretches
% % % %%% first term corresponds to in-plane deformations (stretch and shear)
% % % %%% vectorized version
% % % DGT1 = kk_kron(m_l1, phys.m_u1) + ...
% % %       kk_kron(m_l2, phys.m_u2);
% % % FT  = kk_transpose(DGT1);
% % % C = kk_mx_mult(FT,DGT1);% right Cauchy-Green deformation tensor
% % % [res] = real(eig3(reshape(C',3,3,size(C,1)))); % uses vectorized Cardan's formula for root finding
% % % res = sort(res',2);
% % % tol = 1e-4;indx = find(res(1,:)<(tol) & res(1,:)>(-tol));res(:,indx(1)) = [];
% % % % tol = 1e-4;indx = find(res(1,:)<(1+tol) & res(1,:)>(1-tol));res(:,indx(1)) = [];
% % % lambda = sqrt(res); % take square root for finding the principal stretches
% % % b3  = (lambda(:,1)-lambda(:,2)).^2./2./lambda(:,1)./lambda(:,2);
% % % ai = real(lambda(:,1).*lambda(:,2) - 1);     % area_invariant (as function of principal stretches
% % % E_stretch5 = phys.k_stretch/2 * ...
% % %     sum((ai.*ai + phys.a3 * ai .*ai.*ai +phys.a4 * ai.*ai.*ai.*ai).*phys.SSn_o.*w) ;
% % % %% %%%Calculate shear energy
% % % ai = (SSn./phys.SSn_o - 1);     % area_invariant
% % % beta = real(b1);E_shear1 = phys.k_shear *sum((beta + phys.b1 .* ai .*beta + phys.b2 .* beta.*beta).*phys.SSn_o.*w) ;
% % % beta = real(b2);E_shear2 = phys.k_shear *sum((beta + phys.b1 .* ai .*beta + phys.b2 .* beta.*beta).*phys.SSn_o.*w) ;
% % % beta = real(b3);E_shear3 = phys.k_shear *sum((beta + phys.b1 .* ai .*beta + phys.b2 .* beta.*beta).*phys.SSn_o.*w) ;
% % % beta = real(b4);E_shear5 = phys.k_shear *sum((beta + phys.b1 .* ai .*beta + phys.b2 .* beta.*beta).*phys.SSn_o.*w) ;
% % % beta = real(b5);E_shear6 = phys.k_shear *sum((beta + phys.b1 .* ai .*beta + phys.b2 .* beta.*beta).*phys.SSn_o.*w) ;
% % % %%
% % % tol = 1e-10;
% % % if E_shear1<tol, E_shear1 = 0;end
% % % if E_shear2<tol, E_shear2 = 0;end
% % % if E_shear3<tol, E_shear3 = 0;end
% % % if E_shear4<tol, E_shear4 = 0;end
% % % if E_shear5<tol, E_shear5 = 0;end
% % % if E_shear6<tol, E_shear6 = 0;end
% % % 
% % % if E_stretch1<tol, E_stretch1 = 0;end
% % % if E_stretch2<tol, E_stretch2 = 0;end
% % % if E_stretch3<tol, E_stretch3 = 0;end
% % % if E_stretch4<tol, E_stretch4 = 0;end
% % % if E_stretch5<tol, E_stretch5 = 0;end
% % % %% reporting
% % % disp('---------Stretch energy based on------------');
% % % disp(['triangles                            : ' num2str(E_stretch3)]);
% % % disp(['original area invariant DG           : ' num2str(E_stretch1)]);
% % % disp(['AMoS p665 DGT1 (using loops) prin str: ' num2str(E_stretch2)]);
% % % disp(['AMoS p665 DGT1 (vectorized)  prin str: ' num2str(E_stretch5)]);
% % % disp(['midplane lagrange strain tensor      : ' num2str(E_stretch4)]);
% % % 
% % % disp('-------- Shear energy based on --------------');
% % % disp(['triangles                                   :  ' num2str(E_shear4)]);
% % % disp(['AMoS p665 (loops)beta from strain-ai from DG:  ' num2str(E_shear1)]);
% % % disp(['AMoS p665 (loops)beta from principal str    :  ' num2str(E_shear2)]);
% % % disp(['AMoS p665 (vectorized) beta from prin str   :  ' num2str(E_shear3)]);
% % % disp(['midplane lagrange strain tensor (loops)     :  ' num2str(E_shear5)]);
% % % disp(['in-plane shear strain as in AoPaS page 214  :  ' num2str(E_shear6)]);