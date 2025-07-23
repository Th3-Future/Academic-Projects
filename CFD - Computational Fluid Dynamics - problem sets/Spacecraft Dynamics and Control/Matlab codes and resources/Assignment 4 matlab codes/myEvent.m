function [value, isterminal, direction] = myEvent1(t, X)
value      = 1.075-X(1);
isterminal = 1;   % will stop the integration if value=0
direction  = 0;  
end