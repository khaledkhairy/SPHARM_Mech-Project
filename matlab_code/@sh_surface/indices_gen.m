function [l,m, flags] = indices_gen(c)
%prepare the indices for the potential calculation  
% We are assuming the same order for the L's and K's as given 
% in the c vector
counter = 1;
dim = length(c);
l = zeros(0,dim);
m = zeros(0,dim);
lval = 1;
mval = 1;
while counter <= dim
   mval = -1*lval;
   for i = 1:(2*lval + 1)
      l(counter) = lval;
      m(counter) = mval;
      mval = mval + 1;
      counter = counter + 1;
   end
   lval = lval + 1;
end
RNDOFF = 1e-10;
flags = ones(1,length(c));
for i = 1:length(c)
   if abs(c(i)) <RNDOFF
      flags(i) = 0;
   end
end
l = [0 l];m = [0 m];
l = l(1:length(c));

m = m(1:length(c));
