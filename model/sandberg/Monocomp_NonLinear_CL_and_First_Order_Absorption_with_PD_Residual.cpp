//////////////////////////////
///         Summary        ///
//////////////////////////////
$PROB
- Author : Michael BURTIN
- Date: 01/04/2025

//////////////////////////////
/// Parameters definition  ///
//////////////////////////////
$PARAM @annotated
V1     : 1.31   : Plamatic volume (L/kg)
Vm     : 57.99  : Maximal elimination rate (mg/h/kg)
Km     : 46.51  : Michaelis constant (mg/L)
Ka     : 9.37   : Absorption rate constant (h^-1)

$PARAM @annotated   // PD parameters
B0   : 6.05 : Initial bacterial count (log10(CFU/ml))
Bmax : 9.3  : Maximal bacterial count (log10(CFU/mL))
Kg   : 0.83 : Bacterial growth rate (h-1)
Emax : 1.56 : Maximal rate of effect (h-1)
EC50 : 1.34 : Necessary concentration to obtain 50% of the maximum effect (mg/L)
Kd   : 0.0  : Natural death rate of bacteria (h-1) // Not estimated as insignificant
Kon  : 0.0  : Rate of adaptive resistance development (h-1) // For now disable
Koff : 0.0  : Rate of adaptive resistance disappearance (h-1) // For now disable

$PARAM @annotated   // Other parameters
MIC  : 2.0  : Minimal inhibitory concentration (mg/L)

////////////////////////////////
/// Compartments definition  ///
////////////////////////////////
$CMT @annotated
A : Depot
CENTRAL    : Plasma with total drug amount (mg/h/kg)

// PD compartments specific to CSF
S         : Active bacteria population (CFU/mL)
Rp        : Resting bacteria population (CFU/mL)
ARoff     : Fictive adaptive resistance compartment
ARon      : Fictive adaptive resistance compartment

AUC_CENTRAL         : Total AUC in central compartment (mg/L*h)
TOVER_MIC_CENTRAL   : T>MIC in central compartment (h)
////////////////////////////////
///       Main function      ///
////////////////////////////////
$GLOBAL
double Cmax_CENTRAL;

$MAIN
S_0 = pow(10, B0);

if(TIME == 0.0) {
  Cmax_CENTRAL = 0;
}

//////////////////////////////
///  Equations definition  ///
//////////////////////////////
$ODE    
// PK Equations
dxdt_A = -Ka*A;
dxdt_CENTRAL = Ka*A - (Vm*(CENTRAL/V1)/(Km + (CENTRAL/V1)));

// PD Equations
double B    = S + Rp;
double E    = (Emax * (1 - ARon) * C_CENTRAL) / (C_CENTRAL + EC50);
double Ksr  = ((Kg - Kd) * B) / pow(10, Bmax);

dxdt_S      = Kg*S - (E + Kd)*S - Ksr*S;
dxdt_Rp     = Ksr*S - Kd*Rp;
dxdt_ARoff  = -Kon*ARoff + Koff*ARon;
dxdt_ARon   = -Koff*ARon + Kon*ARoff;

// PKPD Indexes
Cmax_CENTRAL            = (C_CENTRAL > Cmax_CENTRAL) ? C_CENTRAL : Cmax_CENTRAL;
dxdt_AUC_CENTRAL        = C_CENTRAL;
dxdt_TOVER_MIC_CENTRAL  = (C_CENTRAL > MIC) ? 1 : 0;

$SIGMA @labels PD_RES
0.35

$TABLE
double C_CENTRAL  = CENTRAL/V1;
double Log10CFU   = log10(S + Rp) + PD_RES;

//////////////////////////////
///   Outputs definition   ///
//////////////////////////////
$CAPTURE @annotated
C_CENTRAL     : Concentration in central compartment (mg/L)
Log10CFU      : Variation of the bacterial count (log10, CFU/mL)
Cmax_CENTRAL  : Maximal concentration in central compartment (mg/L)
MIC           : Minimal inhibitory concentration (mg/L)
B0            : Initial bacterial count (log10(CFU/ml))