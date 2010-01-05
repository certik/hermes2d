#include "hermes2d.h"
#include "solver_umfpack.h"
#include <string>

// This example shows how to combine automatic adaptivity with the
// Newton's method for a nonlinear time-dependent PDE system.
// The time discretation is done using implicit Euler method.
// Some problem parameters can be changed below.
// The following PDE's are solved:
// Nernst-Planck (describes the diffusion and migration of charged particles):
// dC/dt - D*div[grad(C)] - K*C*div[grad(phi)]=0
// Poisson equation (describes the electric field):
// - div[grad(phi)] = L*(C - C0)
// The equation variables are phi and C and the system describes the
// migration/diffusion of charged particles due to electric field.
// The simulation domain looks as follows:
//      2
//  ____________
//  |          |
// 1|          |1
//  ____________
//      3
// For the Nernst-Planck equation, all the boundaries are natural i.e. Neumann.
// Which basically means that the normal derivative is 0:
// BC: -D*dC/dn - K*C*dphi/dn = 0
// For Poisson equation, 1 has natural boundary conditions (electric field
// derivative is 0). The voltage is applied to the boundaries 2 and 3.
// However, to make the equations stable, i.e. maintain the charge balance
// inside the domain, the positive voltages on the boundary 2 must be
// described by using Neumann boundary. Boundary 3 is an essential boundary, i.e:
// BC 2: dphi/dn = E_FIELD
// BC 3: phi = 0
// BC 1: dphi/dn = 0


#define SIDE_MARKER 1
#define TOP_MARKER 2
#define BOT_MARKER 3

#define VERBAL
#define NONCONT_OUTPUT

/*** Fundamental coefficients ***/
const double D = 1e-12; 	//[m^2/s]
const double R = 8.31; 		//[J/mol*K]
const double T = 293; 		//[K]
const double F = 96485.3415;	//[s * A / mol]
const double eps = 2.5e-2; 	//[F/m]
const double mu = D / (R * T);
const double z = 1;		//Charge number
const double K = z * mu * F;
const double L =  F / eps;	//[V/m^2]
const double VOLTAGE = 3;	//[V]
const double C_CONC = 1200;	//[mol/m^3]

/* For Neumann boundary */
const double height = 180e-6;	//[m]
const double E_FIELD = VOLTAGE / height;

const int NSTEP = 50;
const double TAU = 0.05;
const double NEWTON_TOL = 1e-2;

/* Simulation parameters */
const int P_INIT = 2;       	// initial uniform polynomial order
const int REF_INIT = 12;     	// number of initial refinements
const bool MULTIMESH = true;	// Multimesh?

/* Adaptivity parameters */
// Newton parameters
const double NEWTON_TOL_COARSE = 0.05;	// stopping criterion for Newton on coarse mesh
const double NEWTON_TOL_REF = 0.5;	// stopping criterion for Newton on fine mesh

const int UNREF_FREQ = 10;           	// every UNREF_FREQth time step the mesh
									// is unrefined
const double ERR_STOP = 0.5;		// stopping criterion for hp-adaptivity
const double THRESHOLD = 0.3;		// error threshold for element refinement
const int STRATEGY = 0;			// refinement strategy (0, 1, 2, 3 - see adapt_h1.cpp for explanation)
const int H_ONLY = 0;			// if H_ONLY == 0 then full hp-adaptivity takes place, otherwise
					// h-adaptivity is used. Use this parameter to check that indeed adaptive
					// hp-FEM converges much faster than adaptive h-FEM

const int MAX_NDOFS = 5000;		// To prevent adaptivity going on forever.
// Program params
const std::string USE_ADAPTIVE("adapt");


