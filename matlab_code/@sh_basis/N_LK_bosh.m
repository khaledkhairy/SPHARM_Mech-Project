function NLK = N_LK_bosh(L,K)
K = abs(K);
if abs(K)>L, NLK = 0;
else
%NLK = sqrt((2*L+1)/4/pi*factorial(L-K)/factorial(L+K));
NLK = sqrt((2-isequal(K,0))*(2*L+1)*factorial(L-K)/factorial(L+K));
end