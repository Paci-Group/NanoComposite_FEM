Group{
  Sur_Neu_Ele += Region[4];
  Vol_Ele += Region[1];
  Vol_Ele += Region[2];
  Vol_Ele += Region[3];
  Inc_Ele = Region[3];
  Int_Ele = Region[2];
  Mat_Ele = Region[1];
}
Include "/Users/bhenders/Desktop/MgO_Ag_Project/Figures/Epsilon_infty_FEM/Lib_Materials.pro";
Function {
  dn[Region[4]] = 0;
  epsr[Region[1]] = 3.04;
  epsr[Region[2]] = 3.04;
  epsr[Region[3]] = 191;
}
Constraint { { Name ElectricScalarPotential; Case { { Region Region[5]; Value 0; } { Region Region[6]; Value 0.5; } } } }
Constraint { { Name GlobalElectricPotential; Case { } } }
Constraint { { Name GlobalElectricCharge; Case { } } }
Include "/Users/bhenders/Desktop/MgO_Ag_Project/Figures/Epsilon_infty_FEM/Lib_Electrostatics_v.pro";
