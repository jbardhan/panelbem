load fig2

El_pair_solv = (El_pair(:,1) - Eref_pair(:,1))/joulesPerCalorie;
El_pair_highdiel_solv = (El_pair_highdiel(:,1) - Eref_pair_highdiel(:,1))/joulesPerCalorie;
Enl_2_pair_solv = (Enl_pair(:,1) -Eref_pair(:,1))/joulesPerCalorie;
Enl_11_pair_solv = (Enl_pair(:,2) -Eref_pair(:,2))/joulesPerCalorie;
Enl_15_pair_solv = (Enl_pair(:,3) -Eref_pair(:,3))/joulesPerCalorie;

El_self_solv = (El_self(:,1) - Eref_self(:,1))/joulesPerCalorie;
El_self_highdiel_solv = (El_self_highdiel(:,1) - Eref_self_highdiel(:,1))/joulesPerCalorie;
Enl_2_self_solv = (Enl_self(:,1) - Eref_self(:,1))/joulesPerCalorie;
Enl_11_self_solv = (Enl_self(:,2) - Eref_self(:,2))/joulesPerCalorie;
Enl_15_self_solv = (Enl_self(:,3) - Eref_self(:,3))/joulesPerCalorie;

local_burial = (El_pair(:,1)-E_l_part1_ionpair(1))/joulesPerCalorie;
local_highdiel_burial = (El_pair_highdiel(:,1)-E_l_part1_highdiel_ionpair(1))/joulesPerCalorie;
nonlocal_2_burial = (Enl_pair(:,1)-E_nl_part1_ionpair(1))/joulesPerCalorie;
nonlocal_11_burial = (Enl_pair(:,2)-E_nl_part1_ionpair(2))/joulesPerCalorie;
nonlocal_15_burial = (Enl_pair(:,3)-E_nl_part1_ionpair(3))/joulesPerCalorie;

figure;
LH  = plot(distance, local_burial,'r','linewidth',1.8);
set(gca,'fontsize',16);
hold on;
NLH_2 = plot(distance, nonlocal_2_burial,'b--','linewidth',1.8);
NLH_15 = plot(distance, nonlocal_15_burial,'k-.','linewidth',1.8);
L_HD_H = plot(distance, local_highdiel_burial,'m:','linewidth',1.8);
NLH_11 = plot(distance(1:2:end), nonlocal_11_burial(1:2:end),'go','linewidth',1.8);
XH = xlabel('Ion pair position (Angstrom)');
YH = ylabel('Burial Penalty (kcal/mol)');
LegH = legend('Local',...
				  'Nonlocal, \lambda=2',...
				  'Nonlocal, \lambda=15',...
				  'Local, \epsilon =15',...
				  'Nonlocal, \lambda=11',...
				  'location','southwest');

print -depsc2 fig2-pair-burial-penalty.eps
print -dpng fig2-pair-burial-penalty.png
print -dill fig2-pair-burial-penalty.ai

figure;
LH  = plot(distance, El_pair_solv,'r','linewidth',1.8);
set(gca,'fontsize',16);
hold on;
NLH_2 = plot(distance, Enl_2_pair_solv,'b--','linewidth',1.8);
NLH_15 = plot(distance, Enl_15_pair_solv,'k-.','linewidth',1.8);
L_HD_H = plot(distance, El_pair_highdiel_solv,'m:','linewidth',1.8);
XH = xlabel('Ion pair position (Angstrom)');
YH = ylabel('Solvation Free Energy (kcal/mol)');
LegH = legend('Local',...
				  'Nonlocal, \lambda=2',...
				  'Nonlocal, \lambda=15',...
				  'Local, \epsilon =15',...
				  'location','southwest');
axis([0 12 -25 0]);
print -depsc2 fig2-pair-solvation.eps
print -dpng fig2-pair-solvation.png
print -dill fig2-pair-solvation.ai




El_pairwise_interaction = (El_pair_solv(:,1) - 2 * El_self_solv)/2;
El_pairwise_interaction_highdiel = (El_pair_highdiel_solv(:,1) - 2 * El_self_highdiel_solv)/2;
Enl_2_pairwise_interaction = (Enl_2_pair_solv - 2 * Enl_2_self_solv)/2;
Enl_15_pairwise_interaction = (Enl_15_pair_solv - 2 * Enl_15_self_solv)/2;

figure;
LH2 = plot(distance,El_pairwise_interaction, 'r','linewidth',1.8);
set(gca,'fontsize',16);
hold on;
NLH2 = plot(distance, Enl_2_pairwise_interaction, 'b--','linewidth',1.8);
NLH15 = plot(distance, Enl_15_pairwise_interaction, 'k-.','linewidth',1.8);
L_HD_H2 = plot(distance,El_pairwise_interaction_highdiel, ...
					'm:','linewidth',1.8);
XH2 = xlabel('Ion pair position (Angstrom)');
YH2 = ylabel('Screened Pairwise Energy (kcal/mol)');
LegH2 = legend('Local',...
					'Nonlocal, \lambda=2',...
					'Nonlocal, \lambda=15',...
					'Local, \epsilon =15',...
					'location','southeast');
print -depsc2 fig2-pairwise.eps
print -dpng fig2-pairwise.png
print -dill fig2-pairwise.ai

