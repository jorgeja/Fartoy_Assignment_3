function [beta] = sideslip(v)
U = sqrt(v(1)^2+v(2)^2);
beta = asin(v(2)/U);
end

