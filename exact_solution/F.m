%势函数
function r = F(phi)
global epsilon S2
  r =1/(4*epsilon^2)*((phi.^2-1).^2-2*S2*phi.^2);
end
