#ifndef __HERMES2D_SEA_BREEZE_PARAMS_H
#define __HERMES2D_SEA_BREEZE_PARAMS_H

// Empirical relations of initial distributions valid for 0 <= z <= 10km
// "z" in p_z and T_z is in "m", so everything is in SI units
// presure p in Pascals (empirical formula)
#define p_z(z) {{ params.p_z }}
// temperature T in Kelvin (calculated by equations.py)
#define T_z(z) {{ params.T_z }}
// gas density (calculated by equations.py)
#define rho_z(z) {{ params.rho_z }}
#define R {{ params.R }}            // Gas constant
#define g {{ params.g }}           // gravitational acceleration
#define c_v {{ params.c_v }}            // specific heat capacity

#endif
