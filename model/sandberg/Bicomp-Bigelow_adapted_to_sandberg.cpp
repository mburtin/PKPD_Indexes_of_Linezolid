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
V1     : 0.27  : Plamatic volume (L/kg)
V2     : 0.402  : Peripheral volume (L/kg)
Q      : 0.504  : Blood flow rate (L/h/kg)
CL     : 0.0649 : Clearance (L/h/kg)
Vmax   : 3.26   : Maximum rate of metabolism (mg/h/kg)
Km     : 26.4   : Michaelis-Menten constant (mg/L)
Ka     : 5.42   : Absorption rate constant (h^-1)
FU     : 0.75   : Unbound fraction

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
A : depot compartment (mg)
CENTRAL    : Plasma with total drug amount (mg/h/kg)
PERIPHERAL : Peripheral compartment with total drug amount (mg/h/kg)

// PD compartments specific to CSF
S         : Active bacteria population in CSF (CFU/mL)
Rp        : Resting bacteria population in CSF (CFU/mL)
ARoff     : Fictive adaptive resistance compartment for CSF
ARon      : Fictive adaptive resistance compartment for CSF

AUC_CENTRAL        : AUC of muscle concentrations (mg/L*h)
TOVER_MIC_CENTRAL   : Time over MIC in Muscle (h)

////////////////////////////////
///       Main function      ///
////////////////////////////////
$MAIN
F_A = 0.38;

S_0 = pow(10, B0);

if(TIME == 0.0) {
  Cmax_CENTRAL = 0;
}

//////////////////////////////
///  Equations definition  ///
//////////////////////////////
$ODE    
// PK Equations
double k10   = CL/V1; // (h)
double k12   = Q/V1;  // (h)
double k21   = Q/V2;  // (h)

dxdt_A  = -Ka*A;
dxdt_CENTRAL = Ka*A - k12*CENTRAL + k21*PERIPHERAL - k10*CENTRAL - (Vmax*(CENTRAL/V1)/(Km + (CENTRAL/V1)));
dxdt_PERIPHERAL = k12*CENTRAL - k21*PERIPHERAL;

// PD Equations
double B    = S + Rp;
double E    = (Emax * (1 - ARon) * (C_PERIPHERAL*FU)) / ((C_PERIPHERAL*FU) + EC50);
double Ksr  = ((Kg - Kd) * B) / pow(10, Bmax);

dxdt_S      = Kg*S - (E + Kd)*S - Ksr*S;
dxdt_Rp     = Ksr*S - Kd*Rp;
dxdt_ARoff  = -Kon*ARoff + Koff*ARon;
dxdt_ARon   = -Koff*ARon + Kon*ARoff;

// PKPD Indexes
Cmax_CENTRAL           = (C_CENTRAL > Cmax_CENTRAL) ? C_CENTRAL : Cmax_CENTRAL;
dxdt_AUC_CENTRAL       = C_CENTRAL;
dxdt_TOVER_MIC_CENTRAL = (C_CENTRAL > MIC) ? 1 : 0;

$TABLE
double C_CENTRAL  = CENTRAL/V1;
double C_PERIPHERAL  = PERIPHERAL/V2;

//////////////////////////////
///   Outputs definition   ///
//////////////////////////////
$CAPTURE @annotated
C_CENTRAL     : Concentration in central compartment (mg/L)
C_PERIPHERAL  : Concentration in peripheral compartment (mg/L)
Log10CFU      : Variation of the bacterial count (log10, CwFU/mL)
Cmax_CENTRAL  : Maximal concentration in central compartment (mg/L)
MIC           : Minimal inhibitory concentration (mg/L)
B0            : Initial bacterial count (log10(CFU/ml))