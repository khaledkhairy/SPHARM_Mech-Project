function [sf, t, p] = plot_cos_view(s, incr)
%%   generate a 2D plot of the sh field
% Usage [sf, t, p] = plot_cos_view(s, incr);
persistent Y_LK_cv
if nargin==1, incr = 0.05;end
t = 0:incr:pi;
p = -pi:incr:pi;

[t p] = meshgrid(t,p);

%%%% 
[L, K] = shp_surface.indices_gen(1:(s.L_max + 1)^2); M = length(L);N = numel(t);
if size(Y_LK_cv,1)~= N || size(Y_LK_cv,2)~= M
    %disp('Generating new basis ...');
    Y_LK_sf  = zeros(N, M, 'single');
    for S = 1:length(L),
        %disp([num2str(S) ' of ' num2str(length(L))]);
        Y_LK_cv(:,S) = sh_basis.ylk_bosh(L(S),K(S),p(:)',t(:)')'; % uses bosh version
    end
end
sf = Y_LK_cv(:,1:length(s.xc))*s.xc(:);%% generate the  scalar field
tgrid = t; %reshape(t,size(ct));
pgrid = p;%reshape((p-pi), size(ct));
sf = double(reshape(sf,size(t)));
imagesc(t(:), p(:),sf,[50 255]);

% surf((tgrid), (pgrid),sf, 'EdgeColor', 'none');
% %xlim([0 1]);ylim([-1 1]);%
% axis tight;
% view(0,90);
