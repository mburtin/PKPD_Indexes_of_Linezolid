//////////////////////////////
///         Summary        ///
//////////////////////////////
$PROB
- Author : Michael BURTIN
- Date: 06/11/2024
- Description: PKPD model for antibiotics effect on simulated bacteria populations. 
               Calculate also PKPD index, and outputs the PD-Index PKPD graph corresponding.

//////////////////////////////
/// Parameters definition  ///
//////////////////////////////
$PARAM @annotated   // PK parameters
PV1  : 36.4    : Volume of central compartment (L)
PCL  : 8.38    : Clearance (L/h)
V2   : 0.15    : Volume of CSF compartment (L)
PQIN : 0.04    : Rate of drug transfer from plasma to CSF (L/h)
QOUT : 0.0367  : Rate of drug transfer from CSF to plasma (L/h)
QEVD : 0.00764 : Flow rate of the external ventricular drain (L/h) 
FU   : 0.69    : Unbound fraction

$OMEGA @annotated 
EV1  : 0.0534  : ETA on Volume of central compartment
ECL  : 0.206   : ETA on Clearance
EQIN : 0.0225  : ETA on Rate of drug transfer from plasma to CSF

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

$GLOBAL
double Cmax_CENTRAL;
double Cmax_CSF;

////////////////////////////////
/// Compartments definition  ///
////////////////////////////////
$CMT @annotated     
// PK compartments
CENTRAL    : Plasma with total drug amount (mg)
CSF        : Cerebrospinal fluid with drug amount (mg)

// PD compartments specific to CSF
S_CSF     : Active bacteria population in CSF (CFU/mL)
Rp_CSF    : Resting bacteria population in CSF (CFU/mL)
ARoff_CSF : Fictive adaptive resistance compartment for CSF
ARon_CSF  : Fictive adaptive resistance compartment for CSF

// PD compartments specific to Plasma
S_CENTRAL     : Active bacteria population in CENTRAL (CFU/mL)
Rp_CENTRAL    : Resting bacteria population in CENTRAL (CFU/mL)
ARoff_CENTRAL : Fictive adaptive resistance compartment for CENTRAL
ARon_CENTRAL  : Fictive adaptive resistance compartment for CENTRAL

// PKPD Indexes (Cmax calculate in Quarto)
AUCCENTRAL        : AUC of Unbound plasma concentrations (mg/L*h)
AUCCSF            : AUC of CSF concentrations (mg/L*h)
TOVER_MIC_CENTRAL : Time over MIC in Central compartment (h)
TOVER_MIC_CSF     : Time over MIC in CSF compartment (h)
TOVER_4MIC_CENTRAL : Time over 4x MIC in Central compartment (h)
TOVER_4MIC_CSF     : Time over 4x MIC in CSF compartment (h)
TOVER_10MIC_CENTRAL : Time over 10x MIC in Central compartment (h)
TOVER_10MIC_CSF     : Time over 10x MIC in CSF compartment (h)
  
////////////////////////////////
///       Main function      ///
////////////////////////////////
$MAIN
S_CSF_0 = pow(10, B0);
S_CENTRAL_0 = pow(10, B0);

ARoff_CSF_0 = 1;
ARoff_CENTRAL_0 = 1;

double CL  = PCL*exp(ECL);
double V1  = PV1*exp(EV1);
double QIN = PQIN*exp(EQIN);

if(TIME == 0.0) {
  Cmax_CENTRAL = 0;
  Cmax_CSF = 0;
}

//////////////////////////////
///  Equations definition  ///
//////////////////////////////
$ODE    
// PK Equations
double k10   = CL/V1;
double k12   = QIN/V1;
double k21   = QOUT/V2;
double k23   = QEVD/V2;

dxdt_CENTRAL = - k10*CENTRAL - k12*CENTRAL*FU + k21*CSF;
dxdt_CSF     = - k21*CSF - k23*CSF + k12*CENTRAL*FU;

double C_CENTRAL  = (CENTRAL/V1)*FU;  // Unbound concentrations
double C_CSF      = CSF/V2;

if(C_CENTRAL > Cmax_CENTRAL) Cmax_CENTRAL = C_CENTRAL;
if(C_CSF > Cmax_CSF) Cmax_CSF = C_CSF;

// PD Equations for CSF target
double B_CSF    = S_CSF + Rp_CSF;
double E_CSF    = (Emax * (1 - ARon_CSF) * C_CSF) / (C_CSF + EC50);
double Ksr_CSF  = ((Kg - Kd) * B_CSF) / pow(10, Bmax);

dxdt_S_CSF      = Kg*S_CSF - (E_CSF + Kd)*S_CSF - Ksr_CSF*S_CSF;
dxdt_Rp_CSF     = Ksr_CSF*S_CSF - Kd*Rp_CSF;
dxdt_ARoff_CSF  = -Kon*ARoff_CSF + Koff*ARon_CSF;
dxdt_ARon_CSF   = -Koff*ARon_CSF + Kon*ARoff_CSF;

// PD Equations for Plasma target
double B_CENTRAL    = S_CENTRAL + Rp_CENTRAL;
double E_CENTRAL    = (Emax * (1 - ARon_CENTRAL) * C_CENTRAL) / (C_CENTRAL + EC50);
double Ksr_CENTRAL  = ((Kg - Kd) * B_CENTRAL) / pow(10, Bmax);

dxdt_S_CENTRAL      = Kg*S_CENTRAL - (E_CENTRAL + Kd)*S_CENTRAL - Ksr_CENTRAL*S_CENTRAL;
dxdt_Rp_CENTRAL     = Ksr_CENTRAL*S_CENTRAL - Kd*Rp_CENTRAL;
dxdt_ARoff_CENTRAL  = -Kon*ARoff_CENTRAL + Koff*ARon_CENTRAL;
dxdt_ARon_CENTRAL   = -Koff*ARon_CENTRAL + Kon*ARoff_CENTRAL;

// Track AUC in Central & CSF compartment
dxdt_AUCCENTRAL = C_CENTRAL;
dxdt_AUCCSF     = C_CSF;

// Track fraction of time above MIC in Central compartment
if (C_CENTRAL > MIC)
  dxdt_TOVER_MIC_CENTRAL = 1;
else
  dxdt_TOVER_MIC_CENTRAL = 0;

// Track fraction of time above MIC in CSF compartment
if (C_CSF > MIC)
  dxdt_TOVER_MIC_CSF = 1;
else
  dxdt_TOVER_MIC_CSF = 0;

// Track fraction of time above 4x MIC in Central compartment
if (C_CENTRAL > (4*MIC))
  dxdt_TOVER_4MIC_CENTRAL = 1;
else
  dxdt_TOVER_4MIC_CENTRAL = 0;

// Track fraction of time above 4x MIC in CSF compartment
if (C_CSF > (4*MIC))
  dxdt_TOVER_4MIC_CSF = 1;
else
  dxdt_TOVER_4MIC_CSF = 0;

// Track fraction of time above 10x MIC in Central compartment
if (C_CENTRAL > (10*MIC))
  dxdt_TOVER_10MIC_CENTRAL = 1;
else
  dxdt_TOVER_10MIC_CENTRAL = 0;

// Track fraction of time above 10x MIC in CSF compartment
if (C_CSF > (10*MIC))
  dxdt_TOVER_10MIC_CSF = 1;
else
  dxdt_TOVER_10MIC_CSF = 0;

$SIGMA @labels PD_RES
0.38

$TABLE
double Log10CFU_CENTRAL = log10(S_CENTRAL + Rp_CENTRAL) + PD_RES;
double Log10CFU_CSF     = log10(S_CSF + Rp_CSF) + PD_RES;

//////////////////////////////
///   Outputs definition   ///
//////////////////////////////
$CAPTURE @annotated
C_CENTRAL         : Unbound concentration in central compartment (mg/L)
C_CSF             : Concentration in CSF (mg/L)
Log10CFU_CSF      : Total bacterial count in CSF (log10(CFU/ml))
Log10CFU_CENTRAL  : Total bacterial count in central compartment (log10(CFU/ml))
MIC               : Minimal inhibitory concentration (mg/L)
B0                : Initial bacterial count (log10(CFU/ml))
Cmax_CENTRAL      : Maximum free-concentration in central compartment (mg/L)
Cmax_CSF          : Maximum concentration in CSF compartment (mg/L)
