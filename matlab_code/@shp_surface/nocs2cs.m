function S_tr = nocs2cs(S)
%%% convert the "new" CLKs to the old ones
%%% nocs2cs stands for conversion from coefficients that were calculated
%%% whithout the use of the Condon-Shortly phase factor to coefficients
%%% with. The conversion is necessary for compliance with the
%%% derivatives of the basis functions.
%%% Author: Khaled Khairy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S_tr = nocs2cs_single(S(:));


%%
function Xnew = nocs2cs_single(X_o)


nc = length(X_o)/3;yclks = X_o(nc+1:2*nc);zclks = X_o(2*nc+1:end);xclks = X_o(1:nc);
L_max = round(sqrt(length(xclks))-1);
counter = 0;
for L = 0:L_max,
    for K = -L:L,
        counter = counter + 1;
        CS = (-1)^K;
        NLK_old = sqrt((2 * L + 1)/(2*pi*(1+isequal(K,0))) * factorial(L - abs(K))/factorial(L + abs(K)));
        NLK = sh_basis.N_LK_nocs(L,K);
        fac = CS*NLK_old/NLK;
        xclks(counter) = xclks(counter) * 1/fac;
        yclks(counter) = yclks(counter) * 1/fac;
        zclks(counter) = zclks(counter) * 1/fac;
    end
end
Xnew = [xclks(:)' yclks(:)' zclks(:)'];