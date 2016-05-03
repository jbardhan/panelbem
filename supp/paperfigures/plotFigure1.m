%load fig1

local_burial_lowres = (El_part0(:,1) - E_l_part0_singleIon(1))/joulesPerCalorie;
nonlocal_burial_2   = (Enl_part0(:,1) - E_nl_part0_singleIon(1))/joulesPerCalorie;
nonlocal_burial_15   = (Enl_part0(:,2) - E_nl_part0_singleIon(2))/joulesPerCalorie;

local_burial_highres = (El_part0_highres(:,1) - E_l_part0_singleIon(1))/joulesPerCalorie;
nonlocal_burial_2_highres   = (Enl_part0_highres(:,1) - E_nl_part0_singleIon(1))/joulesPerCalorie;
nonlocal_burial_15_highres   = (Enl_part0_highres(:,2) - E_nl_part0_singleIon(2))/joulesPerCalorie;

figure;
L_low  = plot(distance,local_burial_lowres,'r','linewidth',1.8);
hold on;
set(gca,'FontSize',16);
XH = xlabel('Charge Position (Angstrom)');
YH = ylabel('Burial Penalty (kcal/mol)');
NL_2_low = plot(distance,nonlocal_burial_2,'b--','linewidth',1.8);
NL_15_low = plot(distance,nonlocal_burial_15,'k-.','linewidth',1.8);

L_high = plot(distance2,local_burial_highres,'rs','linewidth',1.8);
NL_2_high = plot(distance2,nonlocal_burial_2_highres,'bs','linewidth',1.8);
NL_15_high = plot(distance2,nonlocal_burial_15_highres,'ks','linewidth',1.8);

%LegendH = legend('Local', 'Nonlocal, \lambda = 2',...
%					  'Nonlocal, \lambda = 15','location','southwest');
print -depsc2 fig1-charge-burial-penalty.eps
print -dpng fig1-charge-burial-penalty.png
print -dill fig1-charge-burial-penalty.ai

figure;
E_l_solv = (El_part0(:,1)-Eref_part0(:,1))/joulesPerCalorie;
E_nl_2_solv = (Enl_part0(:,1)-Eref_part0(:,1))/joulesPerCalorie;
E_nl_15_solv = (Enl_part0(:,2)-Eref_part0(:,1))/joulesPerCalorie;
plot(distance, E_l_solv,'r','linewidth',1.8);
hold on;
set(gca,'FontSize',16);
plot(distance, E_nl_2_solv,'b--','linewidth',1.8);
plot(distance, E_nl_15_solv,'k-.','linewidth',1.8);
LegendH = legend('Local', 'Nonlocal, \lambda = 2',...
					  'Nonlocal, \lambda = 15','location','southwest');
XH = xlabel('Charge Position (Angstrom)');
YH = ylabel('Solvation Free Energy (kcal/mol)');
print -depsc2 fig1-solvation.eps
print -dpng fig1-solvation.png
print -dill fig1-solvation.ai
