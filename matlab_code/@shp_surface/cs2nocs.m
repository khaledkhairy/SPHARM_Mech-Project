function S_tr = cs2nocs(S)
%%% convert the "old" CLKs to the new ones
%%% cs2nocs stands for conversion from coefficients that were calculated
%%% whith the use of the Condon-Shortly phase factor to coefficients
%%% without. The conversion is necessary for compliance with the
%%% derivatives of the basis functions.
%%% Author: Khaled Khairy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(S,2)==1,S = S(:)';
    
elseif mod(size(S,2),3) ~= 0,
    S = S(:,1:end-1);
end

    for ix = 1:size(S,1),       % loop over the shapes
        S_tr(ix,:) = cs2nocs_single(S(ix,:));
    end


%%
function Xnew = cs2nocs_single(X_o)
nc = length(X_o)/3;yclks = X_o(nc+1:2*nc);zclks = X_o(2*nc+1:end);xclks = X_o(1:nc);
L_max = round(sqrt(length(xclks))-1);
counter = 0;
for L = 0:L_max,
    for K = -L:L,
        counter = counter + 1;
        CS = (-1)^K;
        NLK_old = sqrt((2 * L + 1)/(2*pi*(1+isequal(K,0))) * factorial(L - abs(K))/factorial(L + abs(K)));
        NLK = sh_basis.N_LK_nocs(L,K);
        xclks(counter) = xclks(counter) * CS*NLK_old/NLK;
        yclks(counter) = yclks(counter) * CS*NLK_old/NLK;
        zclks(counter) = zclks(counter) * CS*NLK_old/NLK;
    end
end
Xnew = [xclks(:)' yclks(:)' zclks(:)'];