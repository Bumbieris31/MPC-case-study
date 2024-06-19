function dxdt = vanHentenODEs(x, u, d, p)
p1 = p(1); p2 = p(2); p3 = p(3); p4 = p(4); p5 = p(5); p6 = p(6);
p7 = p(7); p8 = p(8); p9 = p(9); p10 = p(10); p11 = p(11); p12 = p(12);
p13 = p(13); p14 = p(14); p15 = p(15); p16 = p(16); p17 = p(17);
p18 = p(18); p19 = p(19); p20 = p(20); p21 = p(21); p22 = p(22);
p23 = p(23); p24 = p(24); p25 = p(25); p26 = p(26); p27 = p(27);
p28 = p(28);

w1 = 0; w2 = 0; w3 = 0; w4 = 0; % No disturbances

x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4);

u1        = u(1); % Supply rate of CO2 in mg*m^{-2}*s^{-1}
u2        = u(2); % Ventilation rate through the vents in mm s^{-1}
u3        = u(3); % Energy supply by the heating system in W*m^{-2}

d1        = d(1); % Incoming radiation in W m^{-2}
d2        = d(2); % Outside CO2 in kg m^{-3}
d3        = d(3); % Outdoor tempetarure in ^oC
d4        = d(4); % Outdoor humidity content C_{H2O} in kg m^{-3}

phi       = p4*d1 + (-p5*x3.^2 + p6*x3 - p7)*(x2 - p8);
PhiPhot_c = (1-exp(-p3*x1))*(p4*d1*(-p5*x3.^2 + p6*x3 - p7)*(x2-p8))/phi;         % gross canopy phootsynthesis rate
PhiVent_c = (u2*10^(-3) + p11)*(x2-d2);                                           % mass exhcnage of CO2 thorought the vents
PhiVent_h = (u2*10^(-3) + p11)*(x4 - d4);                                         % canopy transpiration
PhiTransp_h = p21*(1 - exp(-p3*x1))*(p22/(p23*(x3+p24))*exp(p25*x3/(x3+p26))-x4); % mass exchange of H2) through the vents

dx1dt = (p1*PhiPhot_c - p2*x1*2^(x3/10 - 5/2))*(1+w1);
dx2dt = 1/p9*(-PhiPhot_c + p10*x1*2^(x3/10 - 5/2) + u1*10^(-6) - PhiVent_c)*(1+w2);
dx3dt = 1/p16*(u3 - (p17*u2*10^(-3) + p18)*(x3 - d3) + p19*d1)*(1+w3);
dx4dt = 1/p20*(PhiTransp_h - PhiVent_h)*(1+w4);

dxdt = [dx1dt; dx2dt; dx3dt; dx4dt];
end
