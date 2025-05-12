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
Q      : 0.504  : Blood flow rate (L/h/kg)
CL     : 0.0649 : Clearance (L/h/kg)
Vmax   : 3.26   : Maximum rate of metabolism (mg/h/kg)
Km     : 26.4   : Michaelis-Menten constant (mg/L)
Ka     : 7.32   : Absorption rate constant (h^-1)

////////////////////////////////
/// Compartments definition  ///
////////////////////////////////
$CMT @annotated
A : depot compartment (mg)
CENTRAL    : Plasma with total drug amount (mg/h/kg)
PERIPHERAL : Peripheral compartment with total drug amount (mg/h/kg)

////////////////////////////////
///       Main function      ///
////////////////////////////////
$MAIN
if(self.time == 0.0)
{
  if(self.amt <= 10) 
    F_A = 0.184;
  else if(self.amt <= 30)
    F_A = 0.233;
  else
    F_A = 1.0;
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

$TABLE
double C_CENTRAL  = CENTRAL/V1;

//////////////////////////////
///   Outputs definition   ///
//////////////////////////////
$CAPTURE @annotated
C_CENTRAL     : Concentration in central compartment (mg/L)