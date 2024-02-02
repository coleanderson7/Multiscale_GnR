function [Pfo,so,Pfh,sh] = PreStr(PAR,parT0,parT4,riio,emco,riih,emch)
%
%*** computes mechanobiologically equilibrated stresses, pressure, and
%    transducer axial force at original (o) and evolved (h) states
%

%
%** PAR
%
c   = PAR(1);						% c elastin
Get = PAR(2);						% circumferential deposition stretch elastin
Gez = PAR(3);						% axial deposition stretch elastin
Bt  = PAR(4);						% fraction of circumferential collagen within the adventitia
Bz  = PAR(5);						% fraction of axial collagen within the adventitia
alp = PAR(6);						% orientation of diagonal collagen wrt axial direction
%
betaM = [Bz 1-Bz];					% medial betas [bzM 2*bdM]
betaA = [Bt Bz 1-Bt-Bz];			% adventitial betas [btA bzA 2*bdA]
%
%-------------------------------  STATE o  --------------------------------
%
%** parT0
%
c1m = parT0(1);						% c1t muscle
c2m = parT0(2);						% c2t muscle
c1c = parT0(3);						% c1t collagen
c2c = parT0(4);						% c2t collagen
Gm  = parT0(5);						% circumferential deposition stretch (combined medial collagen and smc)
Gc  = parT0(6);						% deposition stretch (collagen)
%
%** emco and riio
%
phiMo = [emco(1:2) emco(3)*betaM];	% local mass fractions of medial [e mt cz 2*cd]
phiAo = [emco(4)   emco(5)*betaA];	% local mass fractions of adventitial [e ct cz 2*cd]
%
rio  = riio(1);						% inner radius
rMAo = riio(2);						% M-A radius
roo  = riio(3);						% outer radius
%
lzoo = 1;							% axial stretch from o
%
hMo = rMAo-rio;						% medial thickness
hAo = roo-rMAo;						% adventitial thickness
%
SMo = pi/lzoo*(rMAo^2-rio^2);		% medial cross-sectional area
SAo = pi/lzoo*(roo^2-rMAo^2);		% adventitial cross-sectional area
%
Ge = [1/Get/Gez 1/Get/Gez Get Get Gez Gez];	% [GerM GerA GetM GetA GezM GezA]
%
stMo = phiMo(1)*c*(Ge(3)^2-Ge(1)^2) + ...
	   phiMo(2)*c1m*(Gm^2-1)*exp(c2m*(Gm^2-1)^2)*Gm^2 + ...
	   phiMo(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*sin(alp)^2;
%
stAo = phiAo(1)*c*(Ge(4)^2-Ge(2)^2) + ...
	   phiAo(2)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2 + ...
	   phiAo(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*sin(alp)^2;
%
szMo = phiMo(1)*c*(Ge(5)^2-Ge(1)^2) + ...
	   phiMo(3)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2 + ...
	   phiMo(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*cos(alp)^2;
%
szAo = phiAo(1)*c*(Ge(6)^2-Ge(2)^2) + ...
	   phiAo(3)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2 + ...
	   phiAo(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*cos(alp)^2;
%
sto = (stMo*hMo+stAo*hAo)/(hMo+hAo);
szo = (szMo*SMo+szAo*SAo)/(SMo+SAo);
%
Pfo(1,1) = (stMo*hMo + stAo*hAo)/rio;			% P at o
Pfo(1,2) = szMo*SMo + szAo*SAo;					% f at o
%
so = [stMo stAo szMo szAo sto szo sto+szo];		% stresses at o
%
%-------------------------------  STATE h  --------------------------------
%
%** parT4
%
c1m = parT4(1);						% c1t muscle
c2m = parT4(2);						% c2t muscle
c1c = parT4(3);						% c1t collagen
c2c = parT4(4);						% c2t collagen
Gm  = parT4(5);						% circumferential deposition stretch (combined medial collagen and smc)
Gc  = parT4(6);						% deposition stretch (collagen)
%
%** emch and riih
%
phiMh = [emch(1:2) emch(3)*betaM];	% local mass fractions of medial [e mt cz 2*cd]
phiAh = [emch(4)   emch(5)*betaA];	% local mass fractions of adventitial [e ct cz 2*cd]
%
rih  = riih(1);						% inner radius
rMAh = riih(2);						% M-A radius
roh  = riih(3);						% outer radius
%
lzoh = 1;							% axial stretch from o
%
hMh = rMAh-rih;						% medial thickness
hAh = roh-rMAh;						% adventitial thickness
%
SMh = pi/lzoh*(rMAh^2-rih^2);		% medial cross-sectional area
SAh = pi/lzoh*(roh^2-rMAh^2);		% adventitial cross-sectional area
%
lrMoh = (rMAh-rih)/(rMAo-rio);		% radial stretch from o to h, media
lrAoh = (roh-rMAh)/(roo-rMAo);		% radial stretch from o to h, adventitia
%
ltMoh = (rMAh+rih)/(rMAo+rio);		% circum. stretch from o to h, media
ltAoh = (roh+rMAh)/(roo+rMAo);		% circum. stretch from o to h, adventitia
%
Ge = [1/Get/Gez*lrMoh 1/Get/Gez*lrAoh Get*ltMoh Get*ltAoh Gez*lzoh Gez*lzoh];	%* Ge*Foh
%
stMh = phiMh(1)*c*(Ge(3)^2-Ge(1)^2) + ...
	   phiMh(2)*c1m*(Gm^2-1)*exp(c2m*(Gm^2-1)^2)*Gm^2 + ...
	   phiMh(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*sin(alp)^2;
%
stAh = phiAh(1)*c*(Ge(4)^2-Ge(2)^2) + ...
	   phiAh(2)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2 + ...
	   phiAh(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*sin(alp)^2;
%
szMh = phiMh(1)*c*(Ge(5)^2-Ge(1)^2) + ...
	   phiMh(3)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2 + ...
	   phiMh(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*cos(alp)^2;
%
szAh = phiAh(1)*c*(Ge(6)^2-Ge(2)^2) + ...
	   phiAh(3)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2 + ...
	   phiAh(4)*c1c*(Gc^2-1)*exp(c2c*(Gc^2-1)^2)*Gc^2*cos(alp)^2;
%
sth = (stMh*hMh+stAh*hAh)/(hMh+hAh);
szh = (szMh*SMh+szAh*SAh)/(SMh+SAh);
%
Pfh(1,1) = (stMh*hMh + stAh*hAh)/rih;			% P at h
Pfh(1,2) = szMh*SMh + szAh*SAh;					% f at h
%
sh = [stMh stAh szMh szAh sth szh sth+szh];		% stresses at h
%
end