@#include "Initialize.mod"
@#define EndoVariables = EndoVariables + [ "PI", "0", "theta^(1/(1-varepsilon))" ]
@#define EndoVariables = EndoVariables + [ "L", "0", "((1+vartheta)/( psi_M_1Pvartheta + 1 + vartheta ))^(1/(1+vartheta))" ]
@#define EndoVariables = EndoVariables + [ "NU", "0", "Inf" ]
@#define EndoVariables = EndoVariables + [ "AUX1", "0", "Inf" ]
@#define ShockProcesses = ShockProcesses + [ "A", "0", "Inf", "A_STEADY", "rho_a", "sigma_a" ]
@#define ShockProcesses = ShockProcesses + [ "Sg", "0", "1", "Sg_STEADY", "rho_g", "sigma_g" ]
@#define ShockProcesses = ShockProcesses + [ "beta", "0", "1", "beta_STEADY", "rho_b", "sigma_b" ]
@#include "CreateShocks.mod"
@#include "ClassifyDeclare.mod"

parameters beta_STEADY, A_STEADY, Sg_STEADY, PI_STEADY_PROP, varepsilon, theta, phi_pi, phi_y, rho_a, rho_g, sigma_g, sigma_a, sigma_m, vartheta, psi_M_1Pvartheta, rho_b, sigma_b;

beta_STEADY = 0.994;
Sg_STEADY = 0.2;
A_STEADY = 1;
theta = 0.75;
varepsilon = 6;
PI_STEADY_PROP = ( 1.005 - 1 ) / ( theta ^ ( - 1 / varepsilon ) - 1 );
vartheta = 1;
psi_M_1Pvartheta = 2 - ( 1 + vartheta );
phi_pi = 1.5;
phi_y = 0.25;
rho_r = 0;
rho_a = 0.9;
rho_g = 0.8;
rho_b = 0.8;
sigma_a = 0.0025;
sigma_m = 0.0025;
sigma_b = 0.26;
sigma_g = 0.0032;

varexo epsilon_m;

var r;

model;
	#PI_STEADY_ = 1 + PI_STEADY_PROP * ( theta ^ ( - 1 / varepsilon ) - 1 );
	#psi_ = psi_M_1Pvartheta + 1 + vartheta;
	@#include "InsertNewModelEquations.mod"
	#Y = (A/NU) * L;
	#Y_LEAD = (A_LEAD/NU_LEAD) * L_LEAD;
	#G = Sg*Y;
	#G_LEAD = Sg_LEAD*Y_LEAD;
	#PI_STAR = (( 1 - theta * (PI^(varepsilon-1)) ) / (1 - theta))^(1/(1-varepsilon));
	#PI_STAR_LEAD = (( 1 - theta * (PI_LEAD^(varepsilon-1)) ) / (1 - theta))^(1/(1-varepsilon));
	#C = Y - G;
	#C_LEAD = Y_LEAD - G_LEAD;
	#W = psi_ * L^vartheta*C / ( 1 - psi_ * L^(1+vartheta) / ( 1+vartheta) ) * ( 1 - psi_ * STEADY_STATE(L)^(1+vartheta) / ( 1+vartheta) );
	#MC = W/A;
	#M = exp(-sigma_m * epsilon_m);
	#R = ( PI_STEADY_ / beta_STEADY ) * ( ((PI/STEADY_STATE(PI))^phi_pi) * ((Y/STEADY_STATE(Y))^phi_y) ) * M;
	#AUX2 = varepsilon / (varepsilon - 1) * AUX1;
	#AUX2_LEAD = varepsilon / (varepsilon - 1) * AUX1_LEAD;
	r = log( R );
	1 = R * beta_LEAD * ( C / C_LEAD ) / PI_LEAD;
	AUX1 = MC * (Y/C) + theta * beta_LEAD * PI_LEAD^(varepsilon) * AUX1_LEAD;
	AUX2 = PI_STAR * ((Y/C) + theta * beta_LEAD * ((PI_LEAD^(varepsilon-1))/PI_STAR_LEAD) * AUX2_LEAD);
	log( NU ) = log( theta * (PI^varepsilon) * NU_LAG + (1 - theta) * PI_STAR^(-varepsilon) );
	#CObs = log( C );
	#GObs = log( G );
	#PIObs = log( PI );
	#RObs = r;
	#LObs = log( L );
	#WObs = log( W );
end;

steady_state_model;
	PI_STEADY = 1 + PI_STEADY_PROP * ( theta ^ ( - 1 / varepsilon ) - 1 );
	psi = psi_M_1Pvartheta + 1 + vartheta;
	@#include "NKTransSteadyState.mod"
	@#include "InsertNewSteadyStateEquations.mod"
	r = log( R_ );
end;

shocks;
	@#include "InsertNewShockBlockLines.mod"
	var epsilon_m = 1;
end;

steady;
check;

stoch_simul( order = 3, pruning, irf = 0, periods = 0 );