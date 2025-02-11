//////////////////////////////
///         Summary        ///
//////////////////////////////
$PROB
- Author : Michael BURTIN
- Date: 27/01/2025

//////////////////////////////
/// Parameters definition  ///
//////////////////////////////
$PARAM @annotated
V1     : 0.268  : Plamatic volume (L/kg)
V2     : 0.402  : Peripheral volume (L/kg)
V3     : 0.4    : Muscle volume (L/kg)
Q      : 0.504  : Blood flow rate (L/h/kg)
CL     : 0.0649 : Clearance (L/h/kg)
Vmax   : 3.26   : Maximum rate of metabolism (mg/h/kg)
Km     : 26.4   : Michaelis-Menten constant (mg/L)
Kp     : 0.63  : Partition coefficient between plasma and muscle
Q2     : 2.592 : Blood flow rate in muscle (L/h/kg)
RBP    : 0.78  : Relative blood flow to the muscle
Ka     : 4.23 : Absorption

////////////////////////////////
/// Compartments definition  ///
////////////////////////////////
$CMT @annotated
A : Depot
CENTRAL    : Plasma with total drug amount (mg/h/kg)
PERIPHERAL : Peripheral compartment with total drug amount (mg/h/kg)
MUSCLE     : Muscle compartment with total drug amount (mg/h/kg)

////////////////////////////////
///       Main function      ///
////////////////////////////////
$MAIN
F_A = 0.46;

//////////////////////////////
///  Equations definition  ///
//////////////////////////////
$ODE    
// PK Equations
double k10   = CL/V1; // (h)
double k12   = Q/V1;  // (h)
double k21   = Q/V2;  // (h)

dxdt_A = -Ka*A;
dxdt_CENTRAL = Ka*A - k12*CENTRAL + k21*PERIPHERAL - k10*CENTRAL - (Vmax*(CENTRAL/V1)/(Km + (CENTRAL/V1)));
dxdt_PERIPHERAL = k12*CENTRAL - k21*PERIPHERAL;
dxdt_MUSCLE = 0.0;

$TABLE
double C_CENTRAL  = CENTRAL/V1;

//////////////////////////////
///   Outputs definition   ///
//////////////////////////////
$CAPTURE @annotated
C_CENTRAL     : Concentration in central compartment (mg/L)