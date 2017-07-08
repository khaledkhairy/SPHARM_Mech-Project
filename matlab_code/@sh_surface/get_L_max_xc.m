function L_max = get_L_max_xc(xc)
%%% Returns  L_max based on the length of the shape vector
L_max = round(sqrt(length(xc))-1);
if L_max<0, warning('Wrong L_max');end