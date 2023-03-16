// Lib_Electrostatics_v.pro
//
// Template library for electrostatics using a scalar electric potential (v)
// formulation, with floating potentials.

// Default definitions of constants, groups and functions that can/should be
// redefined from outside the template:

DefineConstant[
  modelPath = "", // default path of the model
  resPath = StrCat[modelPath, "results_cuboids/"], // path for post-operation files
  eps0 = 8.854187818e-12, // permittivity of vacuum
  Flag_Axi = 0, // axisymmetric model?
  Val_Rint = 0, // internal radius of Vol_Inf_Ele annulus
  Val_Rext = 0, // external radius of Vol_Inf_Ele annulus
  Val_Cx = 0, // x-coordinate of center of Vol_Inf_Ele
  Val_Cy = 0, // y-coordinate of center of Vol_Inf_Ele
  Val_Cz = 0, // z-coordinate of center of Vol_Inf_Ele
  modelDim = 3 // all of our models are 3D
];

Group {
  DefineGroup[
    // Ellipsoidal_Inclusion
    Inc_Ele,
    // Interface
    Int_Ele
    // Matrix
    Mat_Ele

    // Full dielectric domain:
    Vol_Ele,

    // Subsets of Vol_Ele:
    Vol_Q_Ele, // region with imposed free charge density rho[]
    Vol_Inf_Ele, // annulus where a infinite shell transformation is applied

    // Boundaries:
    Sur_Neu_Ele, // surfaces with Neumann boundary conditions (n . d = dn[])
    Sur_C_Ele // boundary of conductors (constant v)
  ];
  Dom_Ele = Region[ {Vol_Ele, Sur_Neu_Ele} ];
}

Function{
  DefineFunction[
    epsr, // relative permittivity (in Vol_Ele)
    rho, // free charge density (in Vol_Q_Ele)
    dn // normal displacement (on Sur_Neu_Ele)
  ];
}

// End of definitions.

Jacobian {
  { Name Vol;
    Case {
      If(Flag_Axi && modelDim < 3)
        { Region Vol_Inf_Ele;
          Jacobian VolAxiSquSphShell{Val_Rint, Val_Rext, Val_Cx, Val_Cy, Val_Cz}; }
        { Region All; Jacobian VolAxiSqu; }
      Else
        { Region Vol_Inf_Ele;
          Jacobian VolSphShell{Val_Rint, Val_Rext, Val_Cx, Val_Cy, Val_Cz}; }
        { Region All; Jacobian Vol; }
      EndIf
    }
  }
  { Name Sur;
    Case {
      If(Flag_Axi && modelDim < 3)
        { Region All; Jacobian SurAxi; }
      Else
        { Region All; Jacobian Sur; }
      EndIf
    }
  }
}

Integration {
  { Name Int;
    Case {
      { Type Gauss;
        Case {
          { GeoElement Point; NumberOfPoints  1; }
          { GeoElement Line; NumberOfPoints  3; }
          { GeoElement Triangle; NumberOfPoints  3; }
          { GeoElement Triangle2; NumberOfPoints  6; }
          { GeoElement Quadrangle; NumberOfPoints  4; }
          { GeoElement Tetrahedron; NumberOfPoints  4; }
          { GeoElement Tetrahedron2; NumberOfPoints  15; }
          { GeoElement Hexahedron; NumberOfPoints  6; }
          { GeoElement Prism; NumberOfPoints  9; }
          { GeoElement Pyramid; NumberOfPoints  8; }
	}
      }
    }
  }
}

FunctionSpace {
  { Name Hgrad_vf_Ele; Type Form0;
    BasisFunction {
      // v = v  s  + v    s
      //      n  n    c,k  c,k
      { Name sn; NameOfCoef vn; Function BF_Node;
        Support Dom_Ele; Entity NodesOf[ All, Not Sur_C_Ele ]; }
      { Name sck; NameOfCoef vck; Function BF_GroupOfNodes;
        Support Vol_Ele; Entity GroupsOfNodesOf[ Sur_C_Ele ]; }
    }
    SubSpace { // only for a PostOperation
      { Name vf; NameOfBasisFunction sck; }
    }
    GlobalQuantity {
      { Name GlobalElectricPotential; Type AliasOf; NameOfCoef vck; }
      { Name GlobalElectricCharge; Type AssociatedWith; NameOfCoef vck; }
    }
    Constraint {
      { NameOfCoef vn;
        EntityType NodesOf; NameOfConstraint ElectricScalarPotential; }

      { NameOfCoef GlobalElectricPotential;
        EntityType GroupsOfNodesOf; NameOfConstraint GlobalElectricPotential; }
      { NameOfCoef GlobalElectricCharge;
        EntityType GroupsOfNodesOf; NameOfConstraint GlobalElectricCharge; }
    }
  }
}

Formulation {
  { Name Electrostatics_vf; Type FemEquation;
    Quantity {
      { Name v; Type Local; NameOfSpace Hgrad_vf_Ele; }
      { Name Q; Type Global;
        NameOfSpace Hgrad_vf_Ele [GlobalElectricCharge]; }
      { Name V; Type Global;
        NameOfSpace Hgrad_vf_Ele [GlobalElectricPotential]; }

      // only for a PostOperation
      { Name vf; Type Local; NameOfSpace Hgrad_vf_Ele [vf]; }
    }
    Equation {
      Integral { [ epsr[] * eps0 * Dof{d v} , {d v} ];
        In Vol_Ele; Jacobian Vol; Integration Int; }

      Integral { [ - rho[], {v} ];
        In Vol_Q_Ele; Jacobian Vol; Integration Int; }

      Integral { [ dn[] , {v} ];
        In Sur_Neu_Ele; Jacobian Sur; Integration Int; }

      GlobalTerm { [ - Dof{Q}, {V} ]; In Sur_C_Ele; }
    }
  }
}

Resolution {
  { Name Electrostatics_v;
    System {
      { Name A; NameOfFormulation Electrostatics_vf; }
    }
    Operation {
      Generate[A]; Solve[A]; SaveSolution[A];
    }
  }
}

PostProcessing {
  { Name Electrostatics_v; NameOfFormulation Electrostatics_vf;
    PostQuantity {
      { Name v; Value {
          Term { [ {v} ]; In Vol_Ele; Jacobian Vol; }
        }
      }
      { Name e; Value {
          Term { [ -{d v} ]; In Vol_Ele; Jacobian Vol; }
        }
      }
      { Name d; Value {
          Term { [ -eps0*epsr[] * {d v} ]; In Vol_Ele; Jacobian Vol; }
        }
      }
      { Name Q; Value {
          Term { [ {Q} ]; In Sur_C_Ele; }
        }
      }
      { Name V; Value {
          Term { [ {V} ]; In Sur_C_Ele; }
        }
      }
      { Name C; Value {
          Term { [ {Q}/{V} ]; In Sur_C_Ele; }
        }
      }
      { Name vf; Value {
          Term { [ {vf} ]; In Vol_Ele; Jacobian Vol; }
        }
      }
      { Name force; Value {
          Integral { [ eps0*epsr[] / 2. * VirtualWork[{d v}] ];
            //In Vol_Ele; // restrict support to speed-up search
            In ElementsOf[Vol_Ele, OnOneSideOf Sur_C_Ele];
            Jacobian Vol; Integration Int;
      	  }
      	}
      }

      // Method from Markel
      { Name polarization; Value {
          Term { [ -eps0*(epsr[]-1) * {d v} ]; In Vol_Ele; Jacobian Vol; }
        }
      }
      { Name poli; Value {
          Term { [ -eps0*(epsr[] - 9.26) * {d v} ]; In Inc_Ele; Jacobian Vol; }
        }
      }
      { Name polint; Value {
          Term { [ -eps0*(epsr[] - 9.26) * {d v} ]; In Int_Ele; Jacobian Vol; }
        }
      }
      { Name polh; Value {
          Term { [ -eps0*(9.26 - 1) * {d v} ]; In Vol_Ele; Jacobian Vol; }
        }
      }
      { Name di; Value {
          Integral { [ -eps0*(epsr[] - 9.26) * CompZ[{d v}] ];
            In Inc_Ele; Jacobian Vol; Integration Int; }
        }
      }
      { Name dint; Value {
          Integral { [ -eps0*(epsr[] - 9.26) * CompZ[{d v}] ];
            In Int_Ele; Jacobian Vol; Integration Int; }
        }
      }
      { Name dh; Value {
          Integral { [ -eps0*(9.26 - 1) * CompZ[{d v}] ];
            In Vol_Ele; Jacobian Vol; Integration Int; }
        }
      }
      { Name dtot; Value {
          Integral { [ -eps0*(epsr[] - 1) * CompZ[{d v}] ];
            In Vol_Ele; Jacobian Vol; Integration Int; }
        }
      }
      // Just the integral of polarization density over the appropriate region
      { Name di_prop; Value {
          Integral { [ -eps0*(epsr[] - 1) * CompZ[{d v}] ];
            In Inc_Ele; Jacobian Vol; Integration Int; }
        }
      }
      { Name dint_prop; Value {
          Integral { [ -eps0*(epsr[] - 1) * CompZ[{d v}] ];
            In Int_Ele; Jacobian Vol; Integration Int; }
        }
      }
      { Name dh_prop; Value {
          Integral { [ -eps0*(epsr[] - 1) * CompZ[{d v}] ];
            In Mat_Ele; Jacobian Vol; Integration Int; }
        }
      }
      { Name energy; Value {
          Integral {  [ eps0*epsr[] / 2. * SquNorm[{d v}] ];
	          In Vol_Ele; Jacobian Vol; Integration Int;
          }
        }
      }
      { Name vol; Value {
          Term {  [ ElementVol[] ]; In Vol_Ele; Jacobian Vol;}
        }
      }
      { Name epv; Value {
          Integral {  [ eps0*epsr[] / 2. * SquNorm[{d v}] / ElementVol[] ];
            In Vol_Ele; Jacobian Vol; Integration Int;
          }
        }
      }
      { Name etet; Value {
          Integral {  [ eps0*epsr[] / 2. * SquNorm[{d v}] ];
            In Vol_Ele; Jacobian Vol; Integration Int;
          }
        }
      }
    }
  }
}

PostOperation {
  { Name Electrostatics_v; NameOfPostProcessing Electrostatics_v;
    Operation {
      CreateDir[resPath];
      // energy per unit volume of each element
      // Print[ epv, OnElementsOf Vol_Ele, File StrCat[resPath, "epv.pos"] ];
      // energy per unit volume of each element (local field only)
      // Print[ epv_local, OnElementsOf Vol_Ele, File StrCat[resPath, "epv_local.pos"] ];

      // Total Energy of model
      Print[ energy[Vol_Ele], OnGlobal,
      SendToServer "Output/Energy [J]"{0}, Color "Ivory" ];

      // Energy Per Tetrahedron
      Print[ etet, OnElementsOf Vol_Ele, File StrCat[resPath, "etet.pos"] ];

      // Dipoles
      // Print[ polarization, OnElementsOf Vol_Ele, File StrCat[resPath, "polarization.pos"] ];
      // Print[ polarization, OnGrid {$A, $B, $C} { -6.463:6.463:1/4, -6.463:6.463:1/4, -6.463:6.463:1/4 },
      // Format SimpleTable, File StrCat[resPath, "polarization.csv"] ];
      // Print[ polarization, OnElementsOf Inc_Ele, File StrCat[resPath, "polinc_prop.pos"] ];
      // Print[ polarization, OnElementsOf Int_Ele, File StrCat[resPath, "polint_prop.pos"] ];
      // Print[ polarization, OnElementsOf Mat_Ele, File StrCat[resPath, "polmat_prop.pos"] ];
      // Print[ poli, OnElementsOf Vol_Ele, File StrCat[resPath, "poli.pos"] ];
      // Print[ polint, OnElementsOf Vol_Ele, File StrCat[resPath, "polint.pos"] ];
      // Print[ polh, OnElementsOf Vol_Ele, File StrCat[resPath, "polh.pos"] ];

      Print[ dtot[Vol_Ele], OnGlobal,
      SendToServer "Output/Total Dipole Moment"{0}, Color "Ivory" ];
      // Print[ dh[Vol_Ele], OnGlobal,
      // SendToServer "Output/Matrix Dipole Moment"{0}, Color "Ivory" ];
      // Print[ di[Inc_Ele], OnRegion Inc_Ele,
      // SendToServer "Output/Inclusion Dipole Moment"{0}, Color "Ivory" ];
      // Print[ dint[Int_Ele], OnRegion Int_Ele,
      // SendToServer "Output/Interface Dipole Moment"{0}, Color "Ivory" ];

      Print[ dh_prop[Mat_Ele], OnGlobal,
      SendToServer "Output/Matrix Dipole Moment Ver. 2"{0}, Color "Ivory" ];
      Print[ di_prop[Inc_Ele], OnRegion Inc_Ele,
      SendToServer "Output/Inclusion Dipole Moment Ver. 2"{0}, Color "Ivory" ];
      Print[ dint_prop[Int_Ele], OnRegion Int_Ele,
      SendToServer "Output/Interface Dipole Moment Ver. 2"{0}, Color "Ivory" ];


      // Electric field
      // Print[ e, OnElementsOf Vol_Ele, File StrCat[resPath, "e.pos"] ];
      // Electrostatic Potential
      // Print[ v, OnElementsOf Vol_Ele, File StrCat[resPath, "v.pos"] ];
      // Volume of Each Element
      Print[ vol, OnElementsOf Inc_Ele, File StrCat[resPath, "vol.pos"] ];

      If(NbrRegions[Sur_C_Ele])
        Print[ Q, OnRegion Sur_C_Ele, File StrCat[resPath, "q.txt"],
          Format Table, SendToServer "}Output/Floating charge [C]" ];
        Print[ V, OnRegion Sur_C_Ele, File StrCat[resPath, "q.txt"],
          Format Table, SendToServer "}Output/Floating potential [V]" ];
      EndIf
    }
  }
}
