function [Y_TT, P_TT] = ylk_cos_sin_dthetatheta_bosh(p, t, L_max, P_T)
%% Pre-calculate the normalized associated Legendre functions
%% up to L = L_max

gdimp = size(p,2);
gdimt = size(t,1);
 if (min(size(p))==1 && max(size(t))==1),
     p = p(:);
 end
%     t = t(:);
%     Y_TT = zeros(gdimp,1, (L_max+1)^2);
% else
    Y_TT = zeros(gdimp, gdimt, (L_max+1)^2);
% end

P_TT = plkt(P_T, L_max); 
counter 	= 0;%
for L = 0:L_max

    for K = -L:L
        counter = counter + 1;
        if K >= 0,
            Y_TT(:,:,counter)    =  P_TT(:,:,get_LK_index_P(L,K)).* cos(K * p);%
        elseif K < 0,
            Y_TT(:,:,counter)    =  P_TT(:,:,get_LK_index_P(L,K)).* sin(abs(K) * p);
        end
     end
end

%%%%%%%%%%%%%%%%%%%%%%%%%
function P_T = plkt(P, L_max)
%%% Calculate the derivative of the associated legendre functions using the
%%% recursion relations from Bosh 2000
P_T = zeros(size(P));


ia = 1;
dpnm2 = P(:,:,2);

for L = 1:L_max
    ia = ia + L;
    temp = P(:,:,ia);
    P_T(:,:,ia) = -sqrt(ia-1).*P(:,:,ia+1);
    fac1 = sqrt(2*L*(L+1));
    for K = 1:L-1
        ix = ia+K;
        fac2 = sqrt((L-K)*(L+K+1));
        P_T(:,:,ix) = 1/2*(fac1*temp-fac2*P(:,:,ix+1));
        temp = P(:,:,ix);
        fac1 = fac2;
    end
    P_T(:,:,ia+L) = sqrt(L/2)*temp;
end
P_T(:,:,3) = dpnm2;

% % % % 
% % %%% unnormalized version
% % ia = 1;
% % for L = 1:L_max
% %     ia = ia + L;
% %     temp = P(:,:,ia);
% %     P_T(:,:,ia) = -P(:,:,ia+1);
% %     for K = 1:L-1
% %         ix = ia+K;
% %         fac = (L+K)*(L-K+1);
% %         P_T(:,:,ix) = 1/2*(fac*temp-P(:,:,ix+1));
% %         temp = P(:,:,ix);
% %     end
% %     P_T(:,:,ia+L) = L*temp;
% % end
