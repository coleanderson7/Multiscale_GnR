%
%**** This script computes G&R for a bilayered thin-walled artery
%
%*** G&R-CMM base formulation from <Latorre and Humphrey - APLBE 2018>
%								   (without inflammation)
%
%*** DTA WT model data from <Latorre and Humphrey - BMMB 2018>
%							 (with additional inflammation)
%
%** marcos.latorre@yale.edu

clearvars
%
%***** INITIALIZATION AT ORIGINAL HOMEOSTATIC STATE (o)
%
load('DTA_cmm.mat')					% load values needed by CMM for G&R
%
rio  = riio(1);						% inner radius at o (s = 0)
rMAo = riio(2);						% M-A radius at o (s = 0)
roo  = riio(3);						% outer radius at o (s = 0)
%
hMo = rMAo-rio;						% medial thickness at o (s = 0)
hAo = roo-rMAo;						% adventitial thickness at o (s = 0)
%
%** PAR
%
ce  = PAR(1);						% elastin modulus
Get = PAR(2);						% circumferential deposition stretch elastin
Gez = PAR(3);						% axial deposition stretch elastin
Bt  = PAR(4);						% fraction of circumferential collagen within the adventitia
Bz  = PAR(5);						% fraction of axial collagen within the adventitia
alp = PAR(6);						% orientation of diagonal collagen wrt axial direction
%
Ge = [1/Get/Gez Get Gez];			% elastin deposition stretches [r t z] in media and adventitia
%
betaM = [Bz 1-Bz];					% medial betas [BzM 2*BdM]
betaA = [Bt Bz 1-Bt-Bz];			% adventitial betas [BtA BzA 2*BdA]
%
%** parT0
%
c1m0 = parT0(1);					% c1t muscle at day 0
c2m0 = parT0(2);					% c2t muscle
c1c0 = parT0(3);					% c1t collagen
c2c0 = parT0(4);					% c2t collagen
Gm0  = parT0(5);					% circumferential deposition stretch (combined medial collagen and smc)
Gc0  = parT0(6);					% deposition stretch (collagen)
%
%** parT4
%
c1m4 = parT4(1);					% c1t muscle at day 28 (week 4)
c2m4 = parT4(2);					% c2t muscle
c1c4 = parT4(3);					% c1t collagen
c2c4 = parT4(4);					% c2t collagen
Gm4  = parT4(5);					% circumferential deposition stretch (combined medial collagen and smc)
Gc4  = parT4(6);					% deposition stretch (collagen)
%
Tqoc = 7;							% collagen half-life (hypertensive)
Tqom = (1/etaq)*Tqoc;				% smooth muscle half-life
kmco = [1/Tqom,1/Tqoc,1/Tqoc,1/Tqoc,...	% ko circ m M, ax c M, diag c M, diag c M
		1/Tqoc,1/Tqoc,1/Tqoc,1/Tqoc];	% ko circ c A, ax c A, diag c A, diag c A
%
kmc = [kmco; kmco];					% rate parameters
%
K = [KicM KscM KfcM];				% medial collagen gains [KicM KscM KfcM]
%
Kim = etaUm*KicM;					% smooth muscle gains [Kim Ksm Kfm]
Ksm = etaUm*KscM;
Kfm = etaUm*KfcM;
%
KicA = etaUc*KicM;					% adventitial collagen gains [KicA KscA KfcA]
KscA = etaUc*KscM;
KfcA = etaUc*KfcM;
%
Ds = Tqom/20;						% G&R time increment
%
numHl = 20;							% number of Half Lives (Tqo) to consider for integration
domain = numHl*Tqom;				% total integration domain
invals = round(domain/Ds);		% Ds-length intervals within integration domain
breaks = invals+1;					% breaks within integration domain
%
t = (0:Ds:domain)';					% integration times for past histories
%
q = exp(-t*kmco);					% survival function
%
SimpsErr = simpsons(q,t(1),t(end),[])./(1./kmco) - ones(size(kmco)); %* ... should go to zero if properly integrated
%
phioM = [emco(1:2) emco(3)*betaM];	% local mass fractions of medial [e mt cz 2*cd]
phioA = [emco(4)   emco(5)*betaA];	% local mass fractions of adventitial [e ct cz 2*cd]
%
rho = 1050;							% arterial mass density
%
rhooM = phioM*rho;					% medial [e mt cz 2*cd] densities
rhooA = phioA*rho;					% advent [e ct cz 2*cd] densities
%
rhoReM = rhooM(1);					% referential medial elastin density (constant)
rhoReA = rhooA(1);					% referential advent elastin density (constant)
%
%* original homeostatic referential density of fibers (smc and collagen)
%
rhoRf = ones(breaks,1)*[rhooM([2,3,4,4]),...	% circ m M, ax c M, diag c M, diag c M
						rhooA([2,3,4,4])];		% circ c A, ax c A, diag c A, diag c A
%
mN = (ones(breaks,1)*kmco).*rhoRf;	% nominal mass production rates
%
Ups = ones(size(rhoRf));			% stimulus function
mR = mN.*Ups;						% referential mass production rate
%
ltTauInvM = ones(breaks,1);			% history of medial circum stretches
ltTauInvA = ones(breaks,1);			% history of advent circum stretches
lzTauInv  = ones(breaks,1);			% history of axial stretches
%
lTauInv = [ltTauInvM ltTauInvA lzTauInv];	% stretches history
%
c1m = c1m0*ones(breaks,1);			% past history of material properties
c2m = c2m0*ones(breaks,1);
c1c = c1c0*ones(breaks,1);
c2c = c2c0*ones(breaks,1);
Gm  =  Gm0*ones(breaks,1);
Gc  =  Gc0*ones(breaks,1);
%
[Pfo,so,Pfh,sh] = PreStr(PAR,parT0,parT4,riio,emco,riih,emch);	% P, f, stresses at o and h
%
Po    = Pfo(1);						% inner pressure at o
sgmIo = so(7);						% tr(sgm) at o
fo    = Pfo(2);						% vessel axial force at o
%
Ph   = Pfh(1);						% inner pressure at h
fctr = Ph/Po;						% relative increase in P
%
Epso = 1;							% relative cardiac output at o
%
days = 28*4;						% total simulation time (days)
%
steps = round(days/Ds+1);			% total G&R steps
%
ltM = 1; ltA = 1;					% circum stretches at o
JM  = 1; JA  = 1;					% Jacobians at o
%
fctxi = 1/3;						% factor for inflammation evolution
%
%** initialize solution arrays
%
s    = 0:Ds:days;					% G&R times
ri   = rio*ones(steps,1);			% inner radius
hM   = hMo*ones(steps,1);			% medial thickness
hA   = hAo*ones(steps,1);			% adventitial thickness
rhoR = ones(steps,1)*[rhooM([2,3,4,4]),...	% circ m M, ax c M, diag c M, diag c M
					  rhooA([2,3,4,4])];	% circ c A, ax c A, diag c A, diag c A
Upsi = ones(steps,8);				% stimulus function
Eps  = Epso*ones(steps,1);			% relative cardiac output Q/Qo
lz   = ones(steps,1);				% axial stretch relative to state o
f    = fo*ones(steps,1);			% axial force
xi   = zeros(steps,1);				% prescribed inflammatory burden
%
%** prescribe inner pressure
%
wks4 = 28;							% 4 weeks = 28 days
SP = wks4/4;						% auxiliary time
%
P(s<=SP) = (1 + 1/2*(fctr-1)*(1-cos(pi*s(s<=SP)/SP)))*Po;
P(s> SP) = fctr*Po;
%
%** prescribe inflammatory burden
%
SX = wks4/4; SX2 = wks4; SX3 = SX2+2*SX;	% auxiliary times
fctrXi = 1/2;								% remaining inflammation
%
xi(s<=SX) = 1/2*(1-cos(pi*s(s<=SX)/SX));
xi(s>SX & s<=SX2) = 1;
xi(s>SX2 & s<=SX3) = 1 - fctrXi/2*(1-cos(pi*(s(s>SX2 & s<=SX3)-SX2)/(2*SX)));
xi(s>SX3) = 1 - fctrXi;
%
%***** COMPUTE GROWTH & REMODELING FOR s > 0
%
options = optimoptions(@fsolve,'Display','none','FunctionTolerance',1e-9,'StepTolerance',1e-9,'OptimalityTolerance',1e-9);
%
j = 0;								% maximum number of local iterations performed
%
for k = 2:steps
	%
	s(k) = s(k-1) + Ds;					% current G&R time
	%
	%** update inflammation dependent properties for hereditary integrals
	%
	c1m = [c1m0 + (xi(k)^fctxi)*(c1m4-c1m0)
		   c1m(1:end-1)];
	c2m = [c2m0 + (xi(k)^fctxi)*(c2m4-c2m0)
		   c2m(1:end-1)];
	c1c = [c1c0 + (xi(k)^fctxi)*(c1c4-c1c0)
		   c1c(1:end-1)];
	c2c = [c2c0 + (xi(k)^fctxi)*(c2c4-c2c0)
		   c2c(1:end-1)];
	Gm  = [ Gm0 + (xi(k)^fctxi)*(Gm4-Gm0)
		    Gm(1:end-1)];
	Gc  = [ Gc0 + (xi(k)^fctxi)*(Gc4-Gc0)
		    Gc(1:end-1)];
	%
	%** shift past history arrays so that current state is located first
	%
	qaux = q;
	qaux(2:end,:) = qaux(1:end-1,:);
	kmc(2,:) = kmc(1,:);
	%
	rhoRf(2:end,:) = rhoRf(1:end-1,:);
	mN(2:end,:) = mN(1:end-1,:);
	mR(2:end,:) = mR(1:end-1,:);
	lTauInv(2:end,:) = lTauInv(1:end-1,:);
	Ups(2:end,:) = Ups(1:end-1,:);
	%
	rhoRfCor = zeros(1,size(rhoRf,2));	% trial referential densities
	%
	i = 0;								% local iteration
	%
	f(k) = f(k-1);						% trial axial force
	%
	%** perform local iterations until desired tolerance is attained
	%
	while norm(rhoRf(1,:)-rhoRfCor)/norm(rhoRf(1,:)) > eps(norm(rhoRf(1,:)))
		%
		%* compute circumferential stretches and axial force
		%
		parL = [P(k),lz(k),sgmIo,JM,JA,rio,hMo,hAo,ce,Ge,rhoReM,rhoReA,alp,rho,t(1),t(end),Ds];
		%
		[ltf,fval,eF] = fsolve(@(ltf) LaplaceCMM(ltf,parL,mR,qaux,kmco,kmc,lTauInv,c1m,c2m,c1c,c2c,Gm,Gc),[ltM ltA f(k)],options);
		%
		ltM  = ltf(1);									% medial circum stretch
		ltA  = ltf(2);									% advent circum stretch
		f(k) = ltf(3);									% vessel axial force
		%
		%* additional variables from converged solution
		%
		lrM = JM/ltM/lz(k);								% medial radial strech
		lrA = JA/ltA/lz(k);								% advent radial strech
		%
		hM(k) = lrM*hMo;								% medial thickness
		hA(k) = lrA*hAo;								% advent thickness
		ri(k) = (ltM*(2*rio+hMo)-hM(k))/2;				% inner radius
		%
		lTauInv(1,:) = [1/ltM	1/ltA	1/lz(k)];		% update current stretches
		%
		%* intramural and shear stress stimuli
		%
		Dsgm = ( P(k)*ri(k)/(hM(k)+hA(k)) + f(k)/(pi*(hM(k)+hA(k))*(2*ri(k)+hM(k)+hA(k))) )/sgmIo-1;
		Dtau = Eps(k)*rio^3/ri(k)^3-1;
		%
		%* rate parameters, survival functions, and stimulus functions
		%
		kmc(1,:)   = kmco*(1+Dsgm^2);
		q(2:end,:) = qaux(2:end,:).*( ones(breaks-1,1) * exp(-Ds*mean(kmc)) );
		Ups(1,:)   = [ 1 +  Kim*Dsgm -  Ksm*Dtau +  Kfm*xi(k) ...
					  (1 + KicM*Dsgm - KscM*Dtau + KfcM*xi(k))*ones(1,3) ...
					  (1 + KicA*Dsgm - KscA*Dtau + KfcA*xi(k))*ones(1,4)];    %* a row array
		%
		mR(1,:) = mN(1,:).*Ups(1,:);					% mass production before mass update
		%
		%* update/integrate mass densities and recompute mass production
		%
		rhoRfCor = simpsons(mR.*q,t(1),t(end),[]);		% corrected referential mass densities
		%
		mN(1,:) = kmc(1,:).*rhoRfCor;					% update nominal rates
		mR(1,:) = mN(1,:).*Ups(1,:);					% update mass production
		%
		rhoRf(1,:) = simpsons(mR.*q,t(1),t(end),[]);	% update referential mass densities
		%
		JM = ( rhoReM + sum(rhoRf(1,1:3)) ) / rho;		% update medial [e mt cz 2*cd] volume ratios
		JA = ( rhoReA + sum(rhoRf(1,5:7)) ) / rho;		% update advent [e ct cz 2*cd] volume ratios
		%
		i = i+1;										% update iteration counter
		%
	end
	%
	%** converged values needed below for plots
	%
	rhoR(k,:) = rhoRf(1,:);				% referential mass densities
	Upsi(k,:) = Ups(1,:);				% stimulus functions
	%
	j = max(i,j);						% update max number of iterations
	%
end
%
%% a few plots
%
scrsz = get(0,'ScreenSize');
figure
set(gcf,'position',[0.2*scrsz(3) 0.25*scrsz(4) 0.6*scrsz(3) 0.5*scrsz(4)])
%
%** inner pressure & cytokine dimensionless factor
%
subplot(4,4,1)
hold on
grid on
plot(s,P/P(1),'k','linew',1)
ylabel('$P/P_o$','interpreter','latex')
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days],'XTickLabel',{'','','','',''})
set(gca,'fontsize',13)
%
subplot(4,4,5)
hold on
grid on
plot(s,xi,'linew',1)
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days],'ylim',[-0.2 1.2])
xlabel('$s$ [days]','interpreter','latex')
ylabel('$\varrho_f/\varrho_{fm}$','interpreter','latex')
set(gca,'fontsize',13)
%
%** stimulus function Upsilon for medial smooth muscle
%
subplot(242)
hold on
grid on
plot(s,Upsi(:,1),'linew',1)
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days])
xlabel('$s$ [days]','interpreter','latex')
ylabel('$\Upsilon^m_M$','interpreter','latex')
set(gca,'fontsize',13)
%
%** referential mass density of medial smc relative to homeostatic
%
subplot(243)
hold on
grid on
plot(s,rhoR(:,1)/rhoR(1,1),'linew',1)
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days])
xlabel('$s$ [days]','interpreter','latex')
ylabel('$\rho^m_{MR}/\rho^m_{Mo}$','interpreter','latex')
set(gca,'fontsize',13)
%
%** medial wall thickness
%
subplot(244)
hold on
grid on
plot(s,hM/hMo,'linew',1)
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days])
xlabel('$s$ [days]','interpreter','latex')
ylabel('$h_M/h_{Mo}$','interpreter','latex')
set(gca,'fontsize',13)
%
%** luminal radius
%
subplot(245)
hold on
grid on
plot(s,ri/rio,'linew',1)
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days])
xlabel('$s$ [days]','interpreter','latex')
ylabel('$a/a_o$','interpreter','latex')
set(gca,'fontsize',13)
%
%** stimulus function Upsilon for adventitial collagen
%
subplot(246)
hold on
grid on
plot(s,Upsi(:,5),'linew',1)
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days])
xlabel('$s$ [days]','interpreter','latex')
ylabel('$\Upsilon^c_A$','interpreter','latex')
set(gca,'fontsize',13)
%
%** referential mass density of adventitial collagen relative to homeostatic
%
subplot(247)
hold on
grid on
plot(s,rhoR(:,5)/rhoR(1,5),'linew',1)
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days])
xlabel('$s$ [days]','interpreter','latex')
ylabel('$\rho^c_{AR}/\rho^c_{Ao}$','interpreter','latex')
set(gca,'fontsize',13)
%
%** adventitial wall thickness
%
subplot(248)
hold on
grid on
plot(s,hA/hAo,'linew',1)
set(gca,'xlim',[0 days],'XTick',[0 days/4 days/2 3*days/4 days])
xlabel('$s$ [days]','interpreter','latex')
ylabel('$h_A/h_{Ao}$','interpreter','latex')
set(gca,'fontsize',13)
%