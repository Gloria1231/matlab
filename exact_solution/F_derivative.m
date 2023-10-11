function r = F_derivative(phi)
global epsilon S2
    r = 1/(epsilon^2)*phi.*(phi.^2 - 1-S2);
end
