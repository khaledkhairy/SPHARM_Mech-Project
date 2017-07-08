function Y = ylk_cos_sin(L,K,phi, theta)
% here we calclulate a real combination of the YLK's for  a specific value of L and K
% this function negates the Condon-Shortly phase factor introduced through
% Matlab\s lengendre function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if K ==0 && L ==0,
	Y = 1/2/sqrt(pi)* ones(size(theta));
else
	NLK = sh_basis.N_LK_nocs((L),(K));
	P_LK = legendre(L,cos(theta(:)'));                     %returns a matrix of dimentions (l+1) x dim1(theta) x dim2(theta)
	P_LK = squeeze(P_LK(abs(K)+1,:,:));
	P_LK = reshape(P_LK,size(theta,1),size(theta,2));
    CS = (-1)^K;
	if K >= 0 
		Y = CS*NLK * P_LK.*cos(K * phi);     %%% 
	else	% K < 0
		Y = CS*NLK * P_LK.*sin(abs(K) * phi);
	end
end



%%%%%%%%%%%% THis is the old version without CS / compatible with the old
%%%%%%%%%%%% shapes calculated and fitted
% % if K ==0 && L ==0,
% % 	Y = 1/2/sqrt(pi)* ones(size(theta));
% % else
% % 	NLK = sqrt((2 * L + 1)/(2*pi*(1+isequal(K,0))) * factorial(L - abs(K))/factorial(L + abs(K)));
% % 	P_LK = legendre(L,cos(theta(:)'));                     %returns a matrix of dimentions (l+1) x dim1(theta) x dim2(theta)
% % 	P_LK = squeeze(P_LK(abs(K)+1,:,:));
% % 	P_LK = reshape(P_LK,size(theta,1),size(theta,2));
% % 	if K >= 0 
% % 		Y = NLK * P_LK.*cos(K * phi);     %%% 
% % 	else	% K < 0
% % 		Y = NLK * P_LK.*sin(abs(K) * phi);
% % 	end
% % end