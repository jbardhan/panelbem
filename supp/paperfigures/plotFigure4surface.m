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
Enl_deprot_l10_solv = squeeze(Enldeprot(:,3,:))-Erefdeprot;
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

surfOrBuried = 2; % 1=buried, 2=surface
plotFigure4generic

