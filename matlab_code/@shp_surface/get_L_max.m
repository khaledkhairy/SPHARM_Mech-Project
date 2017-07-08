function L_max = get_L_max(X_o)
%%% Returns  L_max based on the length of the shape vector
L_max = round(sqrt(length(X_o)/3)-1);
if L_max<0, warning('Wrong L_max');end