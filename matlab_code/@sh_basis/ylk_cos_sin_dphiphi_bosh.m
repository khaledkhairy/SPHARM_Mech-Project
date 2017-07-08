function Y_PP = ylk_cos_sin_dpp_bosh(p, t, L_max, P)
%% Pre-calculate the normalized associated Legendre functions
%% up to L = L_max
gdimp = length(p);
gdimt = length(t);
 if (min(size(p))==1 && max(size(t))==1),
     p = p(:);
 end
     %     t = t(:);
%     Y_PP = zeros(gdimt,1, (L_max+1)^2);
% 
% else
    Y_PP = zeros(gdimp, gdimt, (L_max+1)^2); 
% end


counter = 0;%
for L = 0:L_max
    for K = -L:L
        counter = counter + 1;
        if K == 0,
            Y_PP(:,:,counter) = zeros(size(Y_PP(:,:,counter)));
        elseif K > 0,
            Y_PP(:,:,counter) = P(:,:,get_LK_index_P(L,K)).*-(K^2).*cos(K.* p);%
        elseif K < 0,
            Y_PP(:,:,counter) = P(:,:,get_LK_index_P(L,K)).*-(K^2).*sin(abs(K).* p);%
        end
     end
end




