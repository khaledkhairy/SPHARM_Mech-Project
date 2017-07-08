function [sf, t, p] = plot_cos_view(obj, nsf, incr)
% nsf = the number of the scalar field to be plotted
if nargin==1,nsf = 1;end
if nargin ==2, nsf = 1;incr = 0.05;end
sf = obj.sf{nsf};disp(['Plotting scalar field (cos view):  ' sf{1}]);
s = sf{2};%
[sf, t, p] = plot_cos_view(s, incr);



% %%%%%%%%%
% gdimp = length(obj.basis.p);gdimt = length(obj.basis.t);
% c = s.xc(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h);
% sf = sum(c.*obj.basis.Y,3);
% surf(cos(obj.basis.t),cos(obj.basis.p),sf, 'EdgeColor', 'none');
% xlim([-1 1]);ylim([-1 1]);
% view(0,90);





%% % % %% % % % % %% project the scalar field
% % % % lb = (obj.L_max + 1)*(obj.L_max + 1);
% % % % gdimp = length(obj.basis.p);gdimt = length(obj.basis.t);
% % % % sf = obj.sf{nsf}{2}.xc;
% % % % c = sf(:)';h = length(c);c = c(ones(gdimp*gdimt,1),:);c = reshape(c,gdimp, gdimt, h); % %
% % % % s = sum(c.*(obj.basis.Y(:,:,1:lb)),3);
% % % % 
% % % % %% plot
% % % % surf(cos(obj.basis.t),cos(obj.basis.p),double(s), 'EdgeColor', 'none');
% % % % xlim([-1 1]);ylim([-1 1]);
% % % % view(0,90);
