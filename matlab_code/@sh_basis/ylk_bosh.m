function Y = ylk_bosh(L,K,phi, theta)
% here we calclulate a real combination of the YLK's for  a specific value of L and K
% this function negates the Condon-Shortly phase factor introduced through
% Matlab's lengendre function, and uses the Bosh 2000 normalization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if K ==0 && L ==0,
	Y = sh_basis.N_LK_bosh(0,0) * ones(size(theta));
else
	NLK_bosh = sh_basis.N_LK_bosh((L),(K));
	P_LK = legendre(L,cos(theta(:)'));                     %returns a matrix of dimentions (l+1) x dim1(theta) x dim2(theta)
	P_LK = squeeze(P_LK(abs(K)+1,:,:));
	P_LK = reshape(P_LK,size(theta,1),size(theta,2));
    CS = 1;%(-1)^K;
	if K >= 0 
		Y = CS*NLK_bosh * P_LK.*cos(K * phi);     %%% 
	else	% K < 0
		Y = CS*NLK_bosh * P_LK.*sin(abs(K) * phi);
	end
end
