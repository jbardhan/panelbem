Elprot_solv = Elprot - Erefprot;
Eldeprot_solv = Eldeprot - Erefdeprot;
dE_vac = Erefprot - Erefdeprot;
dE_Coul = E_prot_Coul - E_deprot_Coul;

Enl_prot_l2_solv = squeeze(Enlprot(:,1,:))-Erefprot;
Enl_prot_l5_solv = squeeze(Enlprot(:,2,:))-Erefprot;
Enl_prot_l10_solv = squeeze(Enlprot(:,3,:))-Erefprot;
Enl_prot_l15_solv = squeeze(Enlprot(:,4,:))-Erefprot;

Enl_deprot_l2_solv = squeeze(Enldeprot(:,1,:))-Erefdeprot;
Enl_deprot_l5_solv = squeeze(Enldeprot(:,2,:))-Erefdeprot;
Enl_deprot_l15_solv = squeeze(Enldeprot(:,3,:))-Erefdeprot;
Enl_deprot_l15_solv = squeeze(Enldeprot(:,4,:))-Erefdeprot;

dE_l2_total_nl = -Enl_deprot_l2_solv + (dE_vac+dE_Coul) + ...
	 Enl_prot_l2_solv;
dE_l5_total_nl = -Enl_deprot_l5_solv + (dE_vac+dE_Coul) + ...
	 Enl_prot_l5_solv;
dE_l10_total_nl = -Enl_deprot_l10_solv + (dE_vac+dE_Coul) + ...
	 Enl_prot_l10_solv;
dE_l15_total_nl = -Enl_deprot_l15_solv + (dE_vac+dE_Coul) + ...
	 Enl_prot_l15_solv;
dE_total_l  = -Eldeprot_solv + (dE_vac + dE_Coul) + Elprot_solv;

Enl_l2_res_prot_solv = estimateDetails.prot_nl(1,:)-estimateDetails.prot_ref(1,:);
Enl_l5_res_prot_solv = estimateDetails.prot_nl(2,:)-estimateDetails.prot_ref(2,:);
Enl_l10_res_prot_solv = estimateDetails.prot_nl(3,:)-estimateDetails.prot_ref(3,:);
Enl_l15_res_prot_solv = estimateDetails.prot_nl(4,:)-estimateDetails.prot_ref(4,:);

Enl_l2_res_deprot_solv = estimateDetails.deprot_nl(1,:)-estimateDetails.deprot_ref(1,:);
Enl_l5_res_deprot_solv = estimateDetails.deprot_nl(2,:)-estimateDetails.deprot_ref(2,:);
Enl_l10_res_deprot_solv = estimateDetails.deprot_nl(3,:)-estimateDetails.deprot_ref(3,:);
Enl_l15_res_deprot_solv = estimateDetails.deprot_nl(4,:)-estimateDetails.deprot_ref(4,:);

Elres_prot_solv  = estimateDetails.prot_l(1,:)-estimateDetails.prot_ref(1,:);
Elres_deprot_solv  = estimateDetails.deprot_l(1,:)-estimateDetails.deprot_ref(1,:);

dE_res_vac = estimateDetails.prot_ref(1,:)  - estimateDetails.deprot_ref(1,:);
dE_res_Coul = estimateDetails.dE_Coul(1,:);
dE_res_total_l  = -Elres_deprot_solv + (dE_res_vac+dE_res_Coul) ...
	 + Elres_prot_solv;

dE_l2_res_total_nl = -Enl_l2_res_deprot_solv + (dE_res_vac+dE_res_Coul) ...
	 + Enl_l2_res_prot_solv;
dE_l5_res_total_nl = -Enl_l5_res_deprot_solv + (dE_res_vac+dE_res_Coul) ...
	 + Enl_l5_res_prot_solv;
dE_l10_res_total_nl = -Enl_l10_res_deprot_solv + (dE_res_vac+dE_res_Coul) ...
	 + Enl_l10_res_prot_solv;
dE_l15_res_total_nl = -Enl_l15_res_deprot_solv + (dE_res_vac+dE_res_Coul) ...
	 + Enl_l15_res_prot_solv;

T = 300;
RT = 0.0098721*T;
scaleFact = 1/(joulesPerCalorie*log(10)*RT);

figure;
LH1 = plot(E_omega_vec,-(dE_total_l(1,:)-dE_res_total_l)*scaleFact, ...
			  'bo-','linewidth',1.8);
hold on;
set(gca,'FontSize',16);
LH2 = plot(E_omega_vec,-(dE_total_l(2,:)-dE_res_total_l)*scaleFact, ...
			  'ro--','linewidth',1.8);
NLH1 = plot(E_omega_vec,-(dE_l15_total_nl(1,:)-dE_l15_res_total_nl)*scaleFact,...
				'bs-','linewidth',1.8);
NLH2 = plot(E_omega_vec,-(dE_l15_total_nl(2,:)-dE_l15_res_total_nl)*scaleFact,...
				'rs--','linewidth',1.8);
XH1 = xlabel('Solute Dielectric Constant   \epsilon_\Omega');
YH1 = ylabel('\Delta pKa');
LegH1 = legend('Local, Buried','Local, Surface','Nonlocal, Buried',...
					'Nonlocal, Surface','location','northwest');
axis([1.8 40 -1.5 5])
figure;
LH1 = plot(E_omega_vec,-(dE_total_l(1,:)-dE_res_total_l)*scaleFact, ...
			  'bo-','linewidth',1.8);
hold on;
set(gca,'FontSize',16);
LH2 = plot(E_omega_vec,-(dE_total_l(2,:)-dE_res_total_l)*scaleFact, ...
			  'ro--','linewidth',1.8);
NLH1 = plot(E_omega_vec,-(dE_l15_total_nl(1,:)-dE_l15_res_total_nl)*scaleFact,...
				'bs-','linewidth',1.8);
NLH2 = plot(E_omega_vec,-(dE_l15_total_nl(2,:)-dE_l15_res_total_nl)*scaleFact,...
				'rs--','linewidth',1.8);
XH1 = xlabel('Solute Dielectric Constant   \epsilon_\Omega');
YH1 = ylabel('\Delta pKa');
axis([2 5 -0.8 0.8]);

figure;
LH1 = plot(E_omega_vec,-(dE_total_l(1,:)-dE_total_l(2,:))*scaleFact, ...
			  'bo-','linewidth',1.8);
hold on;
NLH1 = plot(E_omega_vec,-(dE_l15_total_nl(1,:)-dE_l15_total_nl(2,:))*scaleFact,...
				'rs--','linewidth',1.8);
set(gca,'FontSize',16);
XH1 = xlabel('Solute Dielectric Constant   \epsilon_\Omega');
YH1 = ylabel('\Delta \Delta pKa');
LegH2 = legend('Local','Nonlocal','location','northeast');
axis([1.8 20 -0.1 +0.4]);

figure;
Lres = plot(E_omega_vec,-(dE_res_total_l-dE_Coul(1,:))/ ...
			  joulesPerCalorie,'r-','linewidth',1.8);
hold on; set(gca,'Fontsize',16);
Lbur = plot(E_omega_vec,-(dE_total_l(1,:)-dE_Coul(1,:))/ ...
			  joulesPerCalorie,'k--','linewidth',1.8);
Lsrf = plot(E_omega_vec,-(dE_total_l(2,:)-dE_Coul(1,:))/ ...
			  joulesPerCalorie,'b-.','linewidth',1.8);
LEGH = legend('Local, Residue', 'Local, Buried', 'Local, Surface', ...
				  'location','southeast');
XH1 = xlabel('Solute Dielectric Constant   \epsilon_\Omega');
YH1 = ylabel('\Delta G_{solv}(MH^+\rightarrow M)  (kcal/mol)');

figure;
Lres = plot(E_omega_vec,-(dE_l15_res_total_nl-dE_Coul(1,:))/ ...
			  joulesPerCalorie,'r-','linewidth',1.8);
hold on; set(gca,'Fontsize',16);
Lbur = plot(E_omega_vec,-(dE_l15_total_nl(1,:)-dE_Coul(1,:))/ ...
			  joulesPerCalorie,'k--','linewidth',1.8);
Lsrf = plot(E_omega_vec,-(dE_l15_total_nl(2,:)-dE_Coul(1,:))/ ...
			  joulesPerCalorie,'b-.','linewidth',1.8);
LEGH = legend('Nonlocal, Residue', 'Nonlocal, Buried', 'Nonlocal, Surface', ...
				  'location','southeast');
XH1 = xlabel('Solute Dielectric Constant   \epsilon_\Omega');
YH1 = ylabel('\Delta G_{solv}(MH^+\rightarrow M)  (kcal/mol)');
axis([1.8 30 -15 11]);
