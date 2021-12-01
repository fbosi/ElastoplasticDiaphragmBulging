% clc;
% clear;
% close all;
function [w_max,er_center,er_edge,phi,s,er,eth,eh,sr,sth,w,dw,dphi,P] = Elastic_range (so)

% This code is based on strain increment rather than stress increment for
% the apex of the membrane. 

global v E ho ro dr

r = 0.0000001:dr:ro;

e(1) = so/E;   
phi(1) = atan(1/sqrt(3));
er(1) = e(1)*((2-v)/sqrt(3)*cos(phi(1)) - v*sin(phi(1)));
eth(1) = e(1)*((1-2*v)/sqrt(3)*cos(phi(1)) + sin(phi(1)));

p = 10^(round(log10(eth(1)))-2); 
p_inc = 10^(round(log10(eth(1)))-3);

eth(length(r))=1;
count=1;
while eth(length(r))>eth(1)*1e-2
p = p + p_inc;  
count = count+1;
   
for i=1:1:length(r)-1 
  
    ri = r(i);
dphi(i) = ((1/2).*exp((-1).*eth(i)).*ri.^(-1).*e(i).^(-1).*((-1).*3.^(1/2)+2.*v.*(1+v).*cos(phi(i)).*e(i)).^(-1) ...
  .*(2.*3.^(1/2).*cos(phi(i)).*(exp(eth(i))+(-1).*exp(er(i)).* ...
  (1+(-3/4).*e(i).^(-2).*exp((-2).*er(i)+2.*(1+(-1).*v).^(-1).*( ...
  er(i)+eth(i))).*p.^2.*ri.^2.*ro.^(-2).*sec(phi(i)).^2).^(1/2))+e(i) ...
  .*((-3).*exp(eth(i)).*v+exp(er(i)).*(1+(-3/4).*e(i).^(-2).*exp((-2).*er(i)+2.*(1+(-1).*v).^(-1).*(er(i)+eth(i))).*p.^2.* ...
  ri.^2.*ro.^(-2).*sec(phi(i)).^2).^(1/2)+4.*exp(er(i)).*v.*(1+( ...
  -3/4).*e(i).^(-2).*exp((-2).*er(i)+2.*(1+(-1).*v).^(-1).*(er(i)+ ...
  eth(i))).*p.^2.*ri.^2.*ro.^(-2).*sec(phi(i)).^2).^(1/2)+cos(2.*phi(i)).*((-3).*exp(eth(i)).*v+2.*exp(er(i)).*((-1)+2.*v).*(1+( ...
  -3/4).*e(i).^(-2).*exp((-2).*er(i)+2.*(1+(-1).*v).^(-1).*(er(i)+ ...
  eth(i))).*p.^2.*ri.^2.*ro.^(-2).*sec(phi(i)).^2).^(1/2))+(-1).*3.^( ...
  1/2).*exp(eth(i)).*v.*sin(2.*phi(i)))));

de(i) = ((1/4).*exp((-1).*eth(i)).*ri.^(-1).*( ...
  3.^(1/2)+(-2).*v.*(1+v).*cos(phi(i)).*e(i)).^(-1).*(3.^(1/2).*(( ...
  -1)+2.*v).*cos(phi(i))+(-3).*sin(phi(i))).^(-1).*((-1).*e(i).*(( ...
  -6).*v.*cos(phi(i)).*((-3).*exp(eth(i)).*v+exp(er(i)).*((-1) ...
  +2.*v).*(4+(-3).*e(i).^(-2).*exp((-2).*((-1)+v).^(-1).*(v.*er(i) ...
  +eth(i))).*p.^2.*ri.^2.*ro.^(-2).*sec(phi(i)).^2).^(1/2))+3.*cos( ...
  3.*phi(i)).*(2.*exp(eth(i)).*((-2)+v).*v+exp(er(i)).*((-1)+ ...
  2.*v).*(4+(-3).*e(i).^(-2).*exp((-2).*((-1)+v).^(-1).*(v.*er(i)+ ...
  eth(i))).*p.^2.*ri.^2.*ro.^(-2).*sec(phi(i)).^2).^(1/2))+3.^(1/2).* ...
  ((-12).*exp(eth(i)).*v.^2+exp(er(i)).*((-1)+(-2).*v+8.*v.^2) ...
  .*(4+(-3).*e(i).^(-2).*exp((-2).*((-1)+v).^(-1).*(v.*er(i)+eth(i))).*p.^2.*ri.^2.*ro.^(-2).*sec(phi(i)).^2).^(1/2)+2.*cos(2.*phi(i)).*((-6).*exp(eth(i)).*v.^2+exp(er(i)).*(1+(-2).*v).^2.*( ...
  4+(-3).*e(i).^(-2).*exp((-2).*((-1)+v).^(-1).*(v.*er(i)+eth(i))) ...
  .*p.^2.*ri.^2.*ro.^(-2).*sec(phi(i)).^2).^(1/2))).*sin(phi(i)))+6.* ...
  cos(phi(i)).*(2.*exp(eth(i))+(-1).*exp(er(i)).*(4+(-3).*e(i).^( ...
  -2).*exp((-2).*((-1)+v).^(-1).*(v.*er(i)+eth(i))).*p.^2.* ...
  ri.^2.*ro.^(-2).*sec(phi(i)).^2).^(1/2)).*sin(phi(i)).*(1+(-2).*v+ ...
  3.^(1/2).*tan(phi(i)))));

    e(i+1) = e(i) + de(i)*dr;
    phi(i+1) = phi(i) + dphi(i)*dr;
    er(i+1) = e(i+1)*((2-v)/sqrt(3)*cos(phi(i+1)) - v*sin(phi(i+1)));
    eth(i+1) = e(i+1)*((1-2*v)/sqrt(3)*cos(phi(i+1)) + sin(phi(i+1)));

end
end

s=e*E;

sr=s.*2/sqrt(3).*cos(phi);
sth=s.*(1/sqrt(3)*cos(phi) + sin(phi));

eh = (er + eth)*v/(v-1);
dw = (-1/2).*3.^(1/2).*e.^(-1).*exp((1+(-1).*v).^(-1).*(er+eth)).*p.*r.*ro.^(-1).*sec(phi);
w = cumtrapz(r,dw) - trapz(r,dw);
w_max = w(1);
P=2*p*ho*E/ro;

er_center = er(1);
er_edge = er(length(r));

end
% figure;
% plot(r,er)
% title('Strain in radial direction')
% xlabel('r (m)')
% ylabel('e_r')
% 
% figure;
% plot(r,eth)
% title('Strain in \theta direction')
% xlabel('r (m)')
% ylabel('e_\theta')
% 
% figure;
% plot(r,eh)
% title('Strain in thickness direction')
% xlabel('r (m)')
% ylabel('e_h')
% 
% figure;
% plot(r,sr)
% title('Stress in radial direction')
% xlabel('r (m)')
% ylabel('\sigma_r')
% 
% figure;
% plot(r,sth)
% title('Stress in \theta direction')
% xlabel('r (m)')
% ylabel('\sigma_\theta')
% 
% figure;
% plot(r,phi)
% title('\phi')
% xlabel('r (m)')
% ylabel('\phi')
% 
% figure;
% plot(r,w)
% title('Vertical displacement')
% xlabel('r (m)')
% ylabel('W (m)')