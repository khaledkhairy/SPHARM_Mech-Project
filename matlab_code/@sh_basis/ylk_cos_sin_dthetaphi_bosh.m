function Y_TP = ylk_cos_sin_dthetaphi_bosh(phi, theta, L_max, P_T)
%% Pre-calculate the normalized associated Legendre functions
%% up to L = L_max
gdimp = size(phi,2);
gdimt = size(theta,1);

 if (min(size(phi))==1 && max(size(theta))==1),
     phi = phi(:);
 end
%     theta = theta(:);
%     Y_TP = zeros(gdimt,1, (L_max+1)^2);
% else
    Y_TP 	= zeros(gdimp, gdimt, (L_max+1)^2);
% end



counter 	= 0;%
for L = 0:L_max
    for K = -L:L
        counter = counter + 1;
        if K == 0,
            Y_TP(:,:,counter) = zeros(size(Y_TP(:,:,counter)));    
        elseif K > 0,
            Y_TP(:,:,counter) =  P_T(:,:,get_LK_index_P(L,K)).* (-K).*sin(K * phi);%
        elseif K < 0,
            Y_TP(:,:,counter) =  P_T(:,:,get_LK_index_P(L,K)).* abs(K).* cos(abs(K) * phi);%
        end
     end
end



