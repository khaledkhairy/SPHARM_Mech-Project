function NLK = N_LK_nocs(L,K)
%K = abs(K);
if abs(K)>L, NLK = 0;
else
NLK = sqrt((2*L+1)/4/pi*factorial(L-K)/factorial(L+K));
end

%NLK = sqrt((2 * L + 1)/(2*pi*(1+isequal(K,0))) * factorial(L - abs(K))/factorial(L + abs(K)));