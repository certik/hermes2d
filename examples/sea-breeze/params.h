#ifndef __HERMES2D_SEA_BREEZE_PARAMS_H
#define __HERMES2D_SEA_BREEZE_PARAMS_H

// Empirical relations of initial distributions valid for 0 <= z <= 10km
// "z" in p_z and T_z is in "m", so everything is in SI units
// presure p in Pascals (empirical formula)
#define p_z(z) (100000 - 11.476*(z) + 0.00052954*pow((z),2) - 9.38e-9*pow((z),3))
// temperature T in Kelvin (calculated by equations.py)
#define T_z(z) (297.602407347392 - 0.00668816443401528*(z) + 2.28953333742903e-7*pow((z),2) + 9.61407194208881e-12*pow((z),3))
// gas density (calculated by equations.py)
#define rho_z(z) (1.17022632601347 - 0.000107996104684066*(z) + 2.86948142331989e-9*pow((z),2))
#define R 287.14            // Gas constant
#define g 9.80665           // gravitational acceleration
#define c_v 20.8            // specific heat capacity

#endif