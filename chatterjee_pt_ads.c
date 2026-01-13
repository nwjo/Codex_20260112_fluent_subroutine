/*
 * Based on Chatterjee et al., 2001
 * Write by NWJO, KATECH
 * Chat-GPT 5-pro used for basic structure
 * Date: 2025-10-29
 * coverage-dependent k^(s) applied
*/

/* Effectiveness Factor for Adsorption Reaction
 * r1   : O2    + Pt(s) + Pt(s) -> O(s)     + O(s)      S0= 7.00E-02
 * r2   : C3H6  + Pt(s)         -> C3H6(s)              S0= 9.80E-01
 * r3   : C3H6  + O(s)  + Pt(s) -> C3H5(s)  + OH(s)     S0= 5.00E-02
 *          Pt(s), mu=-0.9
 * r4   : H2    + Pt(s) + Pt(s) -> H(s)     + H(s)      S0= 4.60E-02
 * r5   : H2O   + Pt(s)         -> H2O(s)               S0= 7.50E-01
 * r6   : CO2   + Pt(s)         -> CO2(s)               S0= 5.00E-03
 * r7   : CO    + Pt(s)         -> CO(s)                S0= 8.40E-01
 * r48  : NO    + Pt(s)         -> NO(s)                S0= 8.50E-01
 * r53  : O2    + Rh(s) + Rh(s) -> O(Rh)    + O(Rh)     S0= 1.00E-02
 * r54  : CO    + Rh(s)         -> CO(Rh)               S0= 5.00E-01
 * r55  : NO    + Rh(s)         -> NO(Rh)               S0= 5.00E-01
*/

#include "udf.h"
#include <math.h>
#include <string.h>


/* Debug/print-gate + UDM utilities removed */

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef CLAMP01
#define CLAMP01(x) ((x) < 0.0 ? 0.0 : ((x) > 1.0 ? 1.0 : (x)))
#endif

/* string match check */
#ifndef STREQ
#define STREQ(a,b) (strcmp((a),(b))==0)
#endif

#define EPS     1.0E-5
#define EPS_FD  1.0e-3    /* C_CO increment (0.1 percent) */

/*** species index mapping ***/
/* indices order
 * gas species --> Pt site species --> Rh site species
 * start from 0
*/

// Gas species indices (at GUI panel, top of species has indices "0")
// Mass fraction at wall are provided in yi[]
#define IDX_O2 0
#define IDX_C3H6 1
#define IDX_H2 2
#define IDX_H2O 3
#define IDX_CO2 4
#define IDX_CO 5
#define IDX_NO 6
#define IDX_NO2 7
#define IDX_N2 8


// Site species index for vacant Pt site fraction theta_Pt(vac) in yi[]
// Fluent surface chemistry uses "site fractions" for site species
#define IDX_Pt_Vac 9
#define IDX_O_Pt 10
#define IDX_C3H6_Pt 11
#define IDX_H_Pt 12
#define IDX_H2O_Pt 13
#define IDX_CO2_Pt 14
#define IDX_CO_Pt 15
#define IDX_NO_Pt 16
#define IDX_N_Pt 17
#define IDX_C3H5_Pt 18
#define IDX_C2H3_Pt 19
#define IDX_CH2_Pt 20
#define IDX_CH3_Pt 21
#define IDX_OH_Pt 22
#define IDX_CH_Pt 23
#define IDX_C_Pt 24
#define IDX_CC2H5_Pt 25
#define IDX_CH3CO_Pt 26
#define IDX_Rh_Vac 27
#define IDX_O_Rh 28
#define IDX_CO_Rh 29
#define IDX_NO_Rh 30
#define IDX_N_Rh 31

/*** Constants ***/
//UNIVERSAL_GAS_CONSTANT: internal alias, no need to define 8314.46261815324 [J/kmol-K] */
#define SITE_DEN_TOT 2.72E-8                /* kmol/m2, =2.72E-9 mol/cm2 */
#define Pt_Frac 0.75                        /* Pt:Rh = 3:1, active surface ratio */
#define SITE_DEN_Pt (SITE_DEN_TOT*Pt_Frac)  /* #define SITE_DEN_Rh (SITE_DEN_TOT*(1.0-Pt_Frac)) */
#define SITE_DEN_Rh (SITE_DEN_TOT*(1.0-Pt_Frac))

// Molecular Weight [kg/kmol]
#define MW_O2   31.9988
#define MW_C3H6 42.081
#define MW_H2   2.01594 
#define MW_H2O  18.01534
#define MW_CO2  44.00995
#define MW_CO   28.01055
#define MW_NO   30.0061
#define MW_NO2  46.0055
#define MW_N2   28.0134


// Sticking Coefficients S0 [dimensionless]
#define S0_O2       7.0E-2  // Reaction #1
#define S0_C3H6     9.8E-1  // Reaction #2
#define S0_C3H6_O   5.00E-2 // Reaction #3
#define S0_H2       4.60E-2 // Reaction #4
#define S0_H2O      7.50E-1 // Reaction #5
#define S0_CO2      5.00E-3 // Reaction #6
#define S0_CO       8.40E-1 // Reaction #7
#define S0_NO       8.50E-1 // Reaction #48
#define S1_O2       1.00E-2 // Reaction #53
#define S1_CO       5.00E-1 // Reaction #54
#define S1_NO       5.00E-1 // Reaction #55

// site requirement = reaction order for adsorption reaction
#define q_R1    2.0
#define q_R2    2.0
#define q_R3    2.0
#define q_R4    2.0
#define q_R5    1.0
#define q_R6    1.0
#define q_R7    1.0
#define q_R48   1.0
#define q_R53   2.0
#define q_R54   1.0
#define q_R55   1.0

// Washcoat Factor
#define Wash_F 70.0                     // Washcoat Factor 

/* Thiele Constants */
#define D_EFF_FIXED (3.40776e-7)        // Effective diffusivity [m^2/s]
#define DELTA_W     (2.5e-5)            // washcoat thickness [m], 25 [um]
#define A_V         (2600)              // from Santos [m^2 catalytic surface/m^3 washcoat]

#define YCO_THIELE_MIN 1.0E-3           // Thiele turn-off threshold
#define CCO_THIELE_MIN 1.0E-5          // Thiele turn-off threshold


#define NU_GAS_R1   1.0             // for gas species order
#define NU_GAS_R2   1.0
#define NU_GAS_R3   1.0
#define NU_GAS_R4   1.0
#define NU_GAS_R5   1.0
#define NU_GAS_R6   1.0
#define NU_GAS_R7   1.0
#define NU_GAS_R48  1.0
#define NU_GAS_R53  1.0
#define NU_GAS_R54  1.0
#define NU_GAS_R55  1.0

static real MW_CO_G = -1.0;                 /* CO mole global storage */
#define CLAMP_POS(x) ((x) > 1e-30 ? (x) : 1e-30)


#define TARGET_WALL_ID 7



static inline real k_surface_covdep(real A, real beta, real Ea_J_per_kmol, real T, const int *idx_site, const real *mu, const real *eps_J_per_kmol, int Ns, const real yi[])
{
    const real Tpos  = MAX(EPS, T);
    const real invRT = 1.0 / (UNIVERSAL_GAS_CONSTANT * Tpos);
    real ln_k = log(MAX(EPS, A)) + beta*log(Tpos) - Ea_J_per_kmol*invRT;
    int i;

    /* coverage loop                */
    /* Ns=0 -> independent coverage */
    for (i = 0; i < Ns; ++i) {
        const int  si   = idx_site ? idx_site[i] : 0; /* Ns=0 -> no loop */
        const real th   = MAX(1.0e-20, CLAMP01(yi[si]));
        const real mui  = mu  ? mu[i]  : 0.0;
        const real epsi = eps_J_per_kmol ? eps_J_per_kmol[i] : 0.0;
        ln_k += mui * log(th) + (epsi * th) * invRT;
    }
    return exp(ln_k);
}

// gas-phase molar concentration at the wall-adjacent cell [kmol/m3]
// using local cell density and wall mass fraction from yi[]
static inline real gas_conc_cell(cell_t c0, Thread *t0, real yi_k, real MW_k)
{
    const real rho = C_R(c0, t0);             /* kg/m^3 */
    return (rho * yi_k) / MAX(EPS, MW_k);     /* kmol/m^3 */
}


/* calculate sticking coefficient */
static inline real A_from_sticking(real S0, real M_kg_per_kmol, real q_site)
{
    const real root = sqrt( UNIVERSAL_GAS_CONSTANT / MAX(EPS, 2.0*M_PI*M_kg_per_kmol) ); /*  [m/s]/sqrt(K) */
    const real invG = pow( 1.0 / MAX(EPS, SITE_DEN_TOT), q_site );
    return S0 * root * invG;  /* [m/s]*(1/gamma^q) */
}





/* ===== reaction parameters -> if no coverage dependecy, Ns=0/NULL ===== */

/* ======================= =================== ========================== */
/* ======================= =================== ========================== */
/* ======================= adsorption reaction ========================== */
/* R1: O2 adsorption on 2 Pt(s) */
static const int  idx_site_r1_ex[]   = {0};
static const real mu_r1_ex[]         = {0.0};
static const real eps_r1_ex[]        = {0.0};
#define idx_site_r1     idx_site_r1_ex
#define mu_r1           mu_r1_ex
#define eps_r1          eps_r1_ex
#define NS_R1   0       /* coverage denpendency species number -> 0 for independent */
#define A1_k    A_from_sticking(S0_O2, MW_O2, q_R1)
#define B1_beta 0.5
#define Ea1_Jpm 0.0     /* [J/kmol] */

/* R2: C3H6 adsorption on 2 Pt(s) */
static const int  idx_site_r2_ex[]   = {0};
static const real mu_r2_ex[]         = {0.0};
static const real eps_r2_ex[]        = {0.0};
#define idx_site_r2     idx_site_r2_ex
#define mu_r2           mu_r2_ex
#define eps_r2          eps_r2_ex
#define NS_R2   0
#define A2_k    A_from_sticking(S0_C3H6, MW_C3H6, q_R2)
#define B2_beta 0.5
#define Ea2_Jpm 0.0     /* [J/kmol] */

/* R3: C3H6 adsorption on 1 Pt(s), 1 O(s) */
static const int  idx_site_r3_ex[] = { IDX_Pt_Vac  };
static const real mu_r3_ex[]       = { -0.9        };
static const real eps_r3_ex[]      = { 0.0         };
#define idx_site_r3     idx_site_r3_ex
#define mu_r3           mu_r3_ex
#define eps_r3          eps_r3_ex
#define NS_R3   1       /* Pt(s) */
#define A3_k    A_from_sticking(S0_C3H6_O, MW_C3H6, q_R3)
#define B3_beta 0.5
#define Ea3_Jpm 0.0     /* [J/kmol] */

/* R4: H2 adsorption on 2 Pt(s) */
static const int  idx_site_r4_ex[]   = { IDX_Pt_Vac };
static const real mu_r4_ex[]         = { -1.0       };
static const real eps_r4_ex[]        = { 0.0        };
#define idx_site_r4     idx_site_r4_ex
#define mu_r4           mu_r4_ex
#define eps_r4          eps_r4_ex
#define NS_R4   1
#define A4_k    A_from_sticking(S0_H2, MW_H2, q_R4)
#define B4_beta 0.5
#define Ea4_Jpm 0.0     /* [J/kmol] */

/* R5: H2O adsorption on 1 Pt(s) */
static const int  idx_site_r5_ex[]   = {0};
static const real mu_r5_ex[]         = {0.0};
static const real eps_r5_ex[]        = {0.0};
#define idx_site_r5     idx_site_r5_ex
#define mu_r5           mu_r5_ex
#define eps_r5          eps_r5_ex
#define NS_R5   0
#define A5_k    A_from_sticking(S0_H2O, MW_H2O, q_R5)
#define B5_beta 0.5
#define Ea5_Jpm 0.0     /* [J/kmol] */

/* R6: CO2 adsorption on 1 Pt(s) */
static const int  idx_site_r6_ex[]   = {0};
static const real mu_r6_ex[]         = {0.0};
static const real eps_r6_ex[]        = {0.0};
#define idx_site_r6     idx_site_r6_ex
#define mu_r6           mu_r6_ex
#define eps_r6          eps_r6_ex
#define NS_R6   0
#define A6_k    A_from_sticking(S0_CO2, MW_CO2, q_R6)
#define B6_beta 0.5
#define Ea6_Jpm 0.0     /* [J/kmol] */

/* R7: CO adsorption on 1 Pt(s) */
static const int  idx_site_r7_ex[]   = {0};
static const real mu_r7_ex[]         = {0.0};
static const real eps_r7_ex[]        = {0.0};
#define idx_site_r7     idx_site_r7_ex
#define mu_r7           mu_r7_ex
#define eps_r7          eps_r7_ex
#define NS_R7   0
#define A7_k    A_from_sticking(S0_CO, MW_CO, q_R7)
#define B7_beta 0.5
#define Ea7_Jpm 0.0     /* [J/kmol] */

/* R48: NO adsorption on 1 Pt(s) */
static const int  idx_site_r48_ex[]   = {0};
static const real mu_r48_ex[]         = {0.0};
static const real eps_r48_ex[]        = {0.0};
#define idx_site_r48     idx_site_r48_ex
#define mu_r48           mu_r48_ex
#define eps_r48          eps_r48_ex
#define NS_R48   0
#define A48_k    A_from_sticking(S0_NO, MW_NO, q_R48)
#define B48_beta 0.5
#define Ea48_Jpm 0.0     /* [J/kmol] */

/* R53: O2 adsorption on 2 Rh(s) */
static const int  idx_site_r53_ex[]   = {IDX_Rh_Vac};
static const real mu_r53_ex[]         = {-1.0};
static const real eps_r53_ex[]        = {0.0};
#define idx_site_r53     idx_site_r53_ex
#define mu_r53           mu_r53_ex
#define eps_r53          eps_r53_ex
#define NS_R53   1      // competition species number
#define A53_k    A_from_sticking(S1_O2, MW_O2, q_R53)
#define B53_beta 0.5
#define Ea53_Jpm 0.0     /* [J/kmol] */

/* R54: CO adsorption on 1 Rh(s) */
static const int  idx_site_r54_ex[]   = {0};
static const real mu_r54_ex[]         = {0.0};
static const real eps_r54_ex[]        = {0.0};
#define idx_site_r54     idx_site_r54_ex
#define mu_r54           mu_r54_ex
#define eps_r54          eps_r54_ex
#define NS_R54   0
#define A54_k    A_from_sticking(S1_CO, MW_CO, q_R54)
#define B54_beta 0.5
#define Ea54_Jpm 0.0     /* [J/kmol] */

/* R55: NO adsorption on 1 Rh(s) */
static const int  idx_site_r55_ex[]   = {0};
static const real mu_r55_ex[]         = {0.0};
static const real eps_r55_ex[]        = {0.0};
#define idx_site_r55     idx_site_r55_ex
#define mu_r55           mu_r55_ex
#define eps_r55          eps_r55_ex
#define NS_R55   0
#define A55_k    A_from_sticking(S1_NO, MW_NO, q_R55)
#define B55_beta 0.5
#define Ea55_Jpm 0.0     /* [J/kmol] */


/* ======================= =================== ========================== */
/* ======================= =================== ========================== */
/* ======================= Desorption Reaction ========================== */

/* R08: 2 O(s) desorption on 2 Pt(s) */
static const int  idx_site_r8_ex[]   = {IDX_O_Pt};
static const real mu_r8_ex[]         = {0.0};
static const real eps_r8_ex[]        = {9.0E7};
#define idx_site_r8     idx_site_r8_ex
#define mu_r8           mu_r8_ex
#define eps_r8          eps_r8_ex
#define NS_R8   1
#define A8_k    3.7E20
#define B8_beta 0.0
#define Ea8_Jpm 2.322E8     /* [J/kmol] */

/* R09: 1 C3H6(s) desorption on 2 Pt(s) */
static const int  idx_site_r9_ex[]   = {0};
static const real mu_r9_ex[]         = {0.0};
static const real eps_r9_ex[]        = {0.0};
#define idx_site_r9     idx_site_r9_ex
#define mu_r9           mu_r9_ex
#define eps_r9          eps_r9_ex
#define NS_R9   0
#define A9_k    1E13
#define B9_beta 0.0
#define Ea9_Jpm 7.27E7     /* [J/kmol] */

/* R10: 1 C3H5(s) desorption on 1 Pt(s) */
static const int  idx_site_r10_ex[]   = {0};
static const real mu_r10_ex[]         = {0.0};
static const real eps_r10_ex[]        = {0.0};
#define idx_site_r10     idx_site_r10_ex
#define mu_r10           mu_r10_ex
#define eps_r10          eps_r10_ex
#define NS_R10   0
#define A10_k    3.7E20
#define B10_beta 0.0
#define Ea10_Jpm 3.1E7     /* [J/kmol] */

/* R11: 2 H(s) desorption on 2 Pt(s) */
static const int  idx_site_r11_ex[]   = {IDX_H_Pt};
static const real mu_r11_ex[]         = {0.0};
static const real eps_r11_ex[]        = {6.0E6};        //fixied 6.0E7 -> 6.0E6
#define idx_site_r11     idx_site_r11_ex
#define mu_r11           mu_r11_ex
#define eps_r11          eps_r11_ex
#define NS_R11   1
#define A11_k    3.7E20
#define B11_beta 0.0
#define Ea11_Jpm 6.74E7     /* [J/kmol] */

/* R12: 1 H2O(s) desorption on 1 Pt(s) */
static const int  idx_site_r12_ex[]   = {0};
static const real mu_r12_ex[]         = {0.0};
static const real eps_r12_ex[]        = {0.0};
#define idx_site_r12     idx_site_r12_ex
#define mu_r12           mu_r12_ex
#define eps_r12          eps_r12_ex
#define NS_R12   0
#define A12_k    1.0E13
#define B12_beta 0.0
#define Ea12_Jpm 4.03E7     /* [J/kmol] */

/* R13: 1 CO(s) desorption on 1 Pt(s) */
static const int  idx_site_r13_ex[]   = {IDX_CO_Pt};
static const real mu_r13_ex[]         = {0.0};
static const real eps_r13_ex[]        = {3.3E7};
#define idx_site_r13     idx_site_r13_ex
#define mu_r13           mu_r13_ex
#define eps_r13          eps_r13_ex
#define NS_R13   1
#define A13_k    1.0E13
#define B13_beta 0.0
#define Ea13_Jpm 1.364E8     /* [J/kmol] */

/* R14: 1 CO2(s) desorption on 1 Pt(s) */
static const int  idx_site_r14_ex[]   = {0};
static const real mu_r14_ex[]         = {0.0};
static const real eps_r14_ex[]        = {0.0};
#define idx_site_r14     idx_site_r14_ex
#define mu_r14           mu_r14_ex
#define eps_r14          eps_r14_ex
#define NS_R14   0
#define A14_k    1.0E13
#define B14_beta 0.0
#define Ea14_Jpm 2.71E7     /* [J/kmol] */

/* ======================= =================== ========================== */
/* ======================= =================== ========================== */
/* ======================= Surface Reaction Pt ========================== */

/* R15: 1 C3H5(s) oxidation on 5 O(s) */
static const int  idx_site_r15_ex[]   = {0};
static const real mu_r15_ex[]         = {0.0};
static const real eps_r15_ex[]        = {0.0};
#define idx_site_r15     idx_site_r15_ex
#define mu_r15           mu_r15_ex
#define eps_r15          eps_r15_ex
#define NS_R15   0
#define A15_k    3.7E16
#define B15_beta 0.0
#define Ea15_Jpm 9.5E7     /* [J/kmol] */

/* ======================= C3H6(s) Consumption ========================== */
/* R16: 1 C3H6(s) consumption */
static const int  idx_site_r16_ex[]   = {0};
static const real mu_r16_ex[]         = {0.0};
static const real eps_r16_ex[]        = {0.0};
#define idx_site_r16     idx_site_r16_ex
#define mu_r16           mu_r16_ex
#define eps_r16          eps_r16_ex
#define NS_R16   0
#define A16_k    1.0E12
#define B16_beta 0.0
#define Ea16_Jpm 7.54E7     /* [J/kmol] */

/* R17: 1 CC2H5(s) + 1 H(s) consumption */
static const int  idx_site_r17_ex[]   = {0};
static const real mu_r17_ex[]         = {0.0};
static const real eps_r17_ex[]        = {0.0};
#define idx_site_r17     idx_site_r17_ex
#define mu_r17           mu_r17_ex
#define eps_r17          eps_r17_ex
#define NS_R17   0
#define A17_k    3.7E20
#define B17_beta 0.0
#define Ea17_Jpm 4.88E7     /* [J/kmol] */

/* R18: 1 CC2H5(s) + 1 Pt(s) consumption */
static const int  idx_site_r18_ex[]   = {0};
static const real mu_r18_ex[]         = {0.0};
static const real eps_r18_ex[]        = {0.0};
#define idx_site_r18     idx_site_r18_ex
#define mu_r18           mu_r18_ex
#define eps_r18          eps_r18_ex
#define NS_R18   0
#define A18_k    3.7E20
#define B18_beta 0.0
#define Ea18_Jpm 1.082E8     /* [J/kmol] */

/* R19: 1 C2H3(s) + CH2(s) consumption */
static const int  idx_site_r19_ex[]   = {0};
static const real mu_r19_ex[]         = {0.0};
static const real eps_r19_ex[]        = {0.0};
#define idx_site_r19     idx_site_r19_ex
#define mu_r19           mu_r19_ex
#define eps_r19          eps_r19_ex
#define NS_R19   0
#define A19_k    3.7E20
#define B19_beta 0.0
#define Ea19_Jpm 3.2E6     /* [J/kmol] */

/* R20: 1 C2H3(s) + Pt(s) consumption */
static const int  idx_site_r20_ex[]   = {0};
static const real mu_r20_ex[]         = {0.0};
static const real eps_r20_ex[]        = {0.0};
#define idx_site_r20     idx_site_r20_ex
#define mu_r20           mu_r20_ex
#define eps_r20          eps_r20_ex
#define NS_R20   0
#define A20_k    3.7E20
#define B20_beta 0.0
#define Ea20_Jpm 4.6E7     /* [J/kmol] */

/* R21: 1 CH3(s) + C(s) consumption */
static const int  idx_site_r21_ex[]   = {0};
static const real mu_r21_ex[]         = {0.0};
static const real eps_r21_ex[]        = {0.0};
#define idx_site_r21     idx_site_r21_ex
#define mu_r21           mu_r21_ex
#define eps_r21          eps_r21_ex
#define NS_R21   0
#define A21_k    3.7E20
#define B21_beta 0.0
#define Ea21_Jpm 4.69E7     /* [J/kmol] */

/* ======================= CHx(s) Consumption ========================== */
/* R22: 1 CH3(s) + Pt(s) consumption */
static const int  idx_site_r22_ex[]   = {0};
static const real mu_r22_ex[]         = {0.0};
static const real eps_r22_ex[]        = {0.0};
#define idx_site_r22     idx_site_r22_ex
#define mu_r22           mu_r22_ex
#define eps_r22          eps_r22_ex
#define NS_R22   0
#define A22_k    1.26E21
#define B22_beta 0.0
#define Ea22_Jpm 7.04E7     /* [J/kmol] */

/* R23: 1 CH2(s) + H(s) consumption */
static const int  idx_site_r23_ex[]   = {0};
static const real mu_r23_ex[]         = {0.0};
static const real eps_r23_ex[]        = {0.0};
#define idx_site_r23     idx_site_r23_ex
#define mu_r23           mu_r23_ex
#define eps_r23          eps_r23_ex
#define NS_R23   0
#define A23_k    3.09E21
#define B23_beta 0.0
#define Ea23_Jpm 0.0     /* [J/kmol] */

/* R24: 1 CH2(s) + Pt(s) consumption */
static const int  idx_site_r24_ex[]   = {0};
static const real mu_r24_ex[]         = {0.0};
static const real eps_r24_ex[]        = {0.0};
#define idx_site_r24     idx_site_r24_ex
#define mu_r24           mu_r24_ex
#define eps_r24          eps_r24_ex
#define NS_R24   0
#define A24_k    7.0E21
#define B24_beta 0.0
#define Ea24_Jpm 5.92E7     /* [J/kmol] */

/* R25: 1 CH(s) + H(s) consumption */
static const int  idx_site_r25_ex[]   = {0};
static const real mu_r25_ex[]         = {0.0};
static const real eps_r25_ex[]        = {0.0};
#define idx_site_r25     idx_site_r25_ex
#define mu_r25           mu_r25_ex
#define eps_r25          eps_r25_ex
#define NS_R25   0
#define A25_k    3.09E21
#define B25_beta 0.0
#define Ea25_Jpm 0.0     /* [J/kmol] */

/* R26: 1 CH(s) + Pt(s) consumption */
static const int  idx_site_r26_ex[]   = {0};
static const real mu_r26_ex[]         = {0.0};
static const real eps_r26_ex[]        = {0.0};
#define idx_site_r26     idx_site_r26_ex
#define mu_r26           mu_r26_ex
#define eps_r26          eps_r26_ex
#define NS_R26   0
#define A26_k    3.09E21
#define B26_beta 0.0
#define Ea26_Jpm 0.0     /* [J/kmol] */

/* R27: 1 C(s) + H(s) consumption */
static const int  idx_site_r27_ex[]   = {0};
static const real mu_r27_ex[]         = {0.0};
static const real eps_r27_ex[]        = {0.0};
#define idx_site_r27     idx_site_r27_ex
#define mu_r27           mu_r27_ex
#define eps_r27          eps_r27_ex
#define NS_R27   0
#define A27_k    1.25E21
#define B27_beta 0.0
#define Ea27_Jpm 1.38E8     /* [J/kmol] */

/* ======================= C2Hx(s) Oxidation ========================== */
/* R28: 1 C2H3(s) + O(s) oxidation */
static const int  idx_site_r28_ex[]   = {0};
static const real mu_r28_ex[]         = {0.0};
static const real eps_r28_ex[]        = {0.0};
#define idx_site_r28     idx_site_r28_ex
#define mu_r28           mu_r28_ex
#define eps_r28          eps_r28_ex
#define NS_R28   0
#define A28_k    3.7E18
#define B28_beta 0.0
#define Ea28_Jpm 6.23E7     /* [J/kmol] */

/* R29: 1 CH3CO(s) + Pt(s) oxidation */
static const int  idx_site_r29_ex[]   = {IDX_O_Pt};
static const real mu_r29_ex[]         = {0.0};
static const real eps_r29_ex[]        = {-4.5E7};
#define idx_site_r29     idx_site_r29_ex
#define mu_r29           mu_r29_ex
#define eps_r29          eps_r29_ex
#define NS_R29   1
#define A29_k    3.7E20
#define B29_beta 0.0
#define Ea29_Jpm 1.967E8     /* [J/kmol] */

/* R30: 1 CH3(s) + CO(s) oxidation */
static const int  idx_site_r30_ex[]   = {0};
static const real mu_r30_ex[]         = {0.0};
static const real eps_r30_ex[]        = {0.0};
#define idx_site_r30     idx_site_r30_ex
#define mu_r30           mu_r30_ex
#define eps_r30          eps_r30_ex
#define NS_R30   0
#define A30_k    3.7E20
#define B30_beta 0.0
#define Ea30_Jpm 8.29E7     /* [J/kmol] */

/* R31: 1 CH3CO(s) + Pt(s) oxidation */
static const int  idx_site_r31_ex[]   = {0};
static const real mu_r31_ex[]         = {0.0};
static const real eps_r31_ex[]        = {0.0};
#define idx_site_r31     idx_site_r31_ex
#define mu_r31           mu_r31_ex
#define eps_r31          eps_r31_ex
#define NS_R31   0
#define A31_k    3.7E20
#define B31_beta 0.0
#define Ea31_Jpm 0.0     /* [J/kmol] */

/* ======================= CHx(s) Oxidation ========================== */
/* R32: 1 CH3(s) + O(s) oxidation */
static const int  idx_site_r32_ex[]   = {0};
static const real mu_r32_ex[]         = {0.0};
static const real eps_r32_ex[]        = {0.0};
#define idx_site_r32     idx_site_r32_ex
#define mu_r32           mu_r32_ex
#define eps_r32          eps_r32_ex
#define NS_R32   0
#define A32_k    3.7E20
#define B32_beta 0.0
#define Ea32_Jpm 3.66E7     /* [J/kmol] */

/* R33: 1 CH2(s) + OH(s) oxidation */
static const int  idx_site_r33_ex[]   = {0};
static const real mu_r33_ex[]         = {0.0};
static const real eps_r33_ex[]        = {0.0};
#define idx_site_r33     idx_site_r33_ex
#define mu_r33           mu_r33_ex
#define eps_r33          eps_r33_ex
#define NS_R33   0
#define A33_k    3.7E20
#define B33_beta 0.0
#define Ea33_Jpm 2.51E7     /* [J/kmol] */

/* R34: 1 CH2(s) + O(s) oxidation */
static const int  idx_site_r34_ex[]   = {0};
static const real mu_r34_ex[]         = {0.0};
static const real eps_r34_ex[]        = {0.0};
#define idx_site_r34     idx_site_r34_ex
#define mu_r34           mu_r34_ex
#define eps_r34          eps_r34_ex
#define NS_R34   0
#define A34_k    3.7E20
#define B34_beta 0.0
#define Ea34_Jpm 2.51E7     /* [J/kmol] */

/* R35: 1 CH(s) + OH(s) oxidation */
static const int  idx_site_r35_ex[]   = {0};
static const real mu_r35_ex[]         = {0.0};
static const real eps_r35_ex[]        = {0.0};
#define idx_site_r35     idx_site_r35_ex
#define mu_r35           mu_r35_ex
#define eps_r35          eps_r35_ex
#define NS_R35   0
#define A35_k    3.7E20
#define B35_beta 0.0
#define Ea35_Jpm 2.52E7     /* [J/kmol] */

/* R36: 1 CH(s) + O(s) oxidation */
static const int  idx_site_r36_ex[]   = {0};
static const real mu_r36_ex[]         = {0.0};
static const real eps_r36_ex[]        = {0.0};
#define idx_site_r36     idx_site_r36_ex
#define mu_r36           mu_r36_ex
#define eps_r36          eps_r36_ex
#define NS_R36   0
#define A36_k    3.7E20
#define B36_beta 0.0
#define Ea36_Jpm 2.51E7     /* [J/kmol] */

/* R37: 1 C(s) + OH(s) oxidation */
static const int  idx_site_r37_ex[]   = {0};
static const real mu_r37_ex[]         = {0.0};
static const real eps_r37_ex[]        = {0.0};
#define idx_site_r37     idx_site_r37_ex
#define mu_r37           mu_r37_ex
#define eps_r37          eps_r37_ex
#define NS_R37   0
#define A37_k    3.7E20
#define B37_beta 0.0
#define Ea37_Jpm 2.248E8     /* [J/kmol] */

/* ======================= H, OH, H2O Oxidation ========================== */
/* R38: 1 O(s) + 1 H(s) oxidation */
static const int  idx_site_r38_ex[]   = {0};
static const real mu_r38_ex[]         = {0.0};
static const real eps_r38_ex[]        = {0.0};
#define idx_site_r38     idx_site_r38_ex
#define mu_r38           mu_r38_ex
#define eps_r38          eps_r38_ex
#define NS_R38   0
#define A38_k    3.7E20
#define B38_beta 0.0
#define Ea38_Jpm 1.15E7     /* [J/kmol] */

/* R39: 1 OH(s) + 1 Pt(s) oxidation */
static const int  idx_site_r39_ex[]   = {0};
static const real mu_r39_ex[]         = {0.0};
static const real eps_r39_ex[]        = {0.0};
#define idx_site_r39     idx_site_r39_ex
#define mu_r39           mu_r39_ex
#define eps_r39          eps_r39_ex
#define NS_R39   0
#define A39_k    5.77E21
#define B39_beta 0.0
#define Ea39_Jpm 7.49E7     /* [J/kmol] */

/* R40: 2 O(s) desorption on 2 Pt(s) */
static const int  idx_site_r40_ex[]   = {0};
static const real mu_r40_ex[]         = {0.0};
static const real eps_r40_ex[]        = {0.0};
#define idx_site_r40     idx_site_r40_ex
#define mu_r40           mu_r40_ex
#define eps_r40          eps_r40_ex
#define NS_R40   0
#define A40_k    3.7E20
#define B40_beta 0.0
#define Ea40_Jpm 1.74E7     /* [J/kmol] */

/* R41: 2 O(s) desorption on 2 Pt(s) */
static const int  idx_site_r41_ex[]   = {0};
static const real mu_r41_ex[]         = {0.0};
static const real eps_r41_ex[]        = {0.0};
#define idx_site_r41     idx_site_r41_ex
#define mu_r41           mu_r41_ex
#define eps_r41          eps_r41_ex
#define NS_R41   0
#define A41_k    3.66E20
#define B41_beta 0.0
#define Ea41_Jpm 7.36E7     /* [J/kmol] */

/* R42: 2 O(s) desorption on 2 Pt(s) */
static const int  idx_site_r42_ex[]   = {0};
static const real mu_r42_ex[]         = {0.0};
static const real eps_r42_ex[]        = {0.0};
#define idx_site_r42     idx_site_r42_ex
#define mu_r42           mu_r42_ex
#define eps_r42          eps_r42_ex
#define NS_R42   0
#define A42_k    3.7E20
#define B42_beta 0.0
#define Ea42_Jpm 4.82E7     /* [J/kmol] */

/* R43: 2 O(s) desorption on 2 Pt(s) */
static const int  idx_site_r43_ex[]   = {0};
static const real mu_r43_ex[]         = {0.0};
static const real eps_r43_ex[]        = {0.0};
#define idx_site_r43     idx_site_r43_ex
#define mu_r43           mu_r43_ex
#define eps_r43          eps_r43_ex
#define NS_R43   0
#define A43_k    2.35E19
#define B43_beta 0.0
#define Ea43_Jpm 4.1E7     /* [J/kmol] */

/* R44: 2 O(s) desorption on 2 Pt(s) */
static const int  idx_site_r44_ex[]   = {IDX_CO_Pt, IDX_NO_Pt};
static const real mu_r44_ex[]         = {0.0, 0.0};
static const real eps_r44_ex[]        = {3.3E7, -9.0E7};
#define idx_site_r44     idx_site_r44_ex
#define mu_r44           mu_r44_ex
#define eps_r44          eps_r44_ex
#define NS_R44   2
#define A44_k    3.7E19
#define B44_beta 0.0
#define Ea44_Jpm 1.08E8     /* [J/kmol] */

/* R45: 2 O(s) desorption on 2 Pt(s) */
static const int  idx_site_r45_ex[]   = {IDX_O_Pt};
static const real mu_r45_ex[]         = {0.0};
static const real eps_r45_ex[]        = {-4.5E7};
#define idx_site_r45     idx_site_r45_ex
#define mu_r45           mu_r45_ex
#define eps_r45          eps_r45_ex
#define NS_R45   1
#define A45_k    3.7E20
#define B45_beta 0.0
#define Ea45_Jpm 1.651E8     /* [J/kmol] */

/* R46: 2 O(s) desorption on 2 Pt(s) */
static const int  idx_site_r46_ex[]   = {IDX_CO_Pt};
static const real mu_r46_ex[]         = {0.0};
static const real eps_r46_ex[]        = {-3.3E7};
#define idx_site_r46     idx_site_r46_ex
#define mu_r46           mu_r46_ex
#define eps_r46          eps_r46_ex
#define NS_R46   1
#define A46_k    3.7E20
#define B46_beta 0.0
#define Ea46_Jpm 0.0     /* [J/kmol] */

/* R47: 2 O(s) desorption on 2 Pt(s) */
static const int  idx_site_r47_ex[]   = {IDX_O_Pt};
static const real mu_r47_ex[]         = {0.0};
static const real eps_r47_ex[]        = {-4.5E7};
#define idx_site_r47     idx_site_r47_ex
#define mu_r47           mu_r47_ex
#define eps_r47          eps_r47_ex
#define NS_R47   1
#define A47_k    3.7E20
#define B47_beta 0.0
#define Ea47_Jpm 2.185E8     /* [J/kmol] */

/* ======================= =================== ========================== */
/* ======================= =================== ========================== */
/* =======================   NO Reduction Pt   ========================== */

/* R49: 2 O(s) desorption on 2 Pt(s) */
static const int  idx_site_r49_ex[]   = {0};
static const real mu_r49_ex[]         = {0.0};
static const real eps_r49_ex[]        = {0.0};
#define idx_site_r49     idx_site_r49_ex
#define mu_r49           mu_r49_ex
#define eps_r49          eps_r49_ex
#define NS_R49   0
#define A49_k    1.0E16
#define B49_beta 0.0
#define Ea49_Jpm 1.4E8     /* [J/kmol] */

/* R50: 2 O(s) desorption on 2 Pt(s) */
static const int  idx_site_r50_ex[]   = {IDX_CO_Pt};
static const real mu_r50_ex[]         = {0.0};
static const real eps_r50_ex[]        = {7.5E7};
#define idx_site_r50     idx_site_r50_ex
#define mu_r50           mu_r50_ex
#define eps_r50          eps_r50_ex
#define NS_R50   1
#define A50_k    3.7E20
#define B50_beta 0.0
#define Ea50_Jpm 1.139E8     /* [J/kmol] */

/* R51: 2 O(s) desorption on 2 Pt(s) */
static const int  idx_site_r51_ex[]   = {IDX_CO_Pt};
static const real mu_r51_ex[]         = {0.0};
static const real eps_r51_ex[]        = {-3.0E6};
#define idx_site_r51     idx_site_r51_ex
#define mu_r51           mu_r51_ex
#define eps_r51          eps_r51_ex
#define NS_R51   1
#define A51_k    5.0E19
#define B51_beta 0.0
#define Ea51_Jpm 1.078E8     /* [J/kmol] */

/* R52: 2 O(s) desorption on 2 Pt(s) */
static const int  idx_site_r52_ex[]   = {IDX_O_Pt};
static const real mu_r52_ex[]         = {0.0};
static const real eps_r52_ex[]        = {4.5E7};
#define idx_site_r52     idx_site_r52_ex
#define mu_r52           mu_r52_ex
#define eps_r52          eps_r52_ex
#define NS_R52   1
#define A52_k    3.7E20
#define B52_beta 0.0
#define Ea52_Jpm 1.281E8     /* [J/kmol] */

/* ======================= =================== ========================== */
/* ======================= =================== ========================== */
/* ======================= NO & CO Reduction Rh========================== */

/* R56: 2 O(s) desorption on 2 Pt(s) */
static const int  idx_site_r56_ex[]   = {0};
static const real mu_r56_ex[]         = {0.0};
static const real eps_r56_ex[]        = {0.0};
#define idx_site_r56     idx_site_r56_ex
#define mu_r56           mu_r56_ex
#define eps_r56          eps_r56_ex
#define NS_R56   0
#define A56_k    3.0E20
#define B56_beta 0.0
#define Ea56_Jpm 2.933E8     /* [J/kmol] */

/* R57: 2 O(s) desorption on 2 Pt(s) */
static const int  idx_site_r57_ex[]   = {IDX_CO_Rh, IDX_N_Rh};
static const real mu_r57_ex[]         = {0.0, 0.0};
static const real eps_r57_ex[]        = {1.88E7, 4.19E7};
#define idx_site_r57     idx_site_r57_ex
#define mu_r57           mu_r57_ex
#define eps_r57          eps_r57_ex
#define NS_R57   2
#define A57_k    1.0E14
#define B57_beta 0.0
#define Ea57_Jpm 1.323E8     /* [J/kmol] */

/* R58: 2 O(s) desorption on 2 Pt(s) */
static const int  idx_site_r58_ex[]   = {0};
static const real mu_r58_ex[]         = {0.0};
static const real eps_r58_ex[]        = {0.0};
#define idx_site_r58     idx_site_r58_ex
#define mu_r58           mu_r58_ex
#define eps_r58          eps_r58_ex
#define NS_R58   0
#define A58_k    5.0E13
#define B58_beta 0.0
#define Ea58_Jpm 1.089E8     /* [J/kmol] */

/* R59: 2 O(s) desorption on 2 Pt(s) */
static const int  idx_site_r59_ex[]   = {IDX_N_Rh};
static const real mu_r59_ex[]         = {0.0};
static const real eps_r59_ex[]        = {1.67E7};
#define idx_site_r59     idx_site_r59_ex
#define mu_r59           mu_r59_ex
#define eps_r59          eps_r59_ex
#define NS_R59   1
#define A59_k    1.11E18
#define B59_beta 0.0
#define Ea59_Jpm 1.369E8     /* [J/kmol] */

/* R60: 2 O(s) desorption on 2 Pt(s) */
static const int  idx_site_r60_ex[]   = {0};
static const real mu_r60_ex[]         = {0.0};
static const real eps_r60_ex[]        = {0.0};
#define idx_site_r60     idx_site_r60_ex
#define mu_r60           mu_r60_ex
#define eps_r60          eps_r60_ex
#define NS_R60   0
#define A60_k    3.7E19
#define B60_beta 0.0
#define Ea60_Jpm 5.99E7     /* [J/kmol] */

/* R61: 2 O(s) desorption on 2 Pt(s) */
static const int  idx_site_r61_ex[]   = {0};
static const real mu_r61_ex[]         = {0.0};
static const real eps_r61_ex[]        = {0.0};
#define idx_site_r61     idx_site_r61_ex
#define mu_r61           mu_r61_ex
#define eps_r61          eps_r61_ex
#define NS_R61   0
#define A61_k    2.22E21
#define B61_beta 0.0
#define Ea61_Jpm 7.95E7     /* [J/kmol] */

static inline void chatterjee_compute_R(cell_t c0, Thread *t0, real Tw, const real yi[], real R[62])
{
    int j;
    for (j = 0; j < 62; ++j) R[j] = 0.0;

    /* Thiele correction excluded (eta = 1.0 for all reactions) */
    const real eta_dummy = 1.0;

    const real theta_pt_vac = CLAMP01(yi[IDX_Pt_Vac]);
    const real theta_rh_vac = CLAMP01(yi[IDX_Rh_Vac]);

    const real Cv_R1 = pow(MAX(EPS, SITE_DEN_Pt * theta_pt_vac), q_R1);    // [kmol/m^2]
    const real Cv_R2 = pow(MAX(EPS, SITE_DEN_Pt * theta_pt_vac), q_R2);
    const real Cv_R3 = pow(MAX(EPS, SITE_DEN_Pt * theta_pt_vac), q_R3);
    const real Cv_R4 = pow(MAX(EPS, SITE_DEN_Pt * theta_pt_vac), q_R4);
    const real Cv_R5 = pow(MAX(EPS, SITE_DEN_Pt * theta_pt_vac), q_R5);
    const real Cv_R6 = pow(MAX(EPS, SITE_DEN_Pt * theta_pt_vac), q_R6);
    const real Cv_R7 = pow(MAX(EPS, SITE_DEN_Pt * theta_pt_vac), q_R7);
    const real Cv_R48 = pow(MAX(EPS, SITE_DEN_Pt * theta_pt_vac), q_R48);
    const real Cv_R53 = pow(MAX(EPS, SITE_DEN_Rh * theta_rh_vac), q_R53);
    const real Cv_R54 = pow(MAX(EPS, SITE_DEN_Rh * theta_rh_vac), q_R54);
    const real Cv_R55 = pow(MAX(EPS, SITE_DEN_Rh * theta_rh_vac), q_R55);

    /* ======================= reaction-01 ======================= */
    {
        const real k1 = k_surface_covdep(A1_k, B1_beta, Ea1_Jpm, Tw, idx_site_r1, mu_r1, eps_r1, NS_R1, yi);

                /* gas o2 concentration [kmol/m3] from wall mass fraction yi[] */
                const real Cg = gas_conc_cell(c0, t0, yi[IDX_O2], MW_O2);

                /* rate per area [kmol/m2 sec]: r = ks * Cg * (C_vac)^q */
                const real rate_base = k1 * Cg * Cv_R1; /* kmol/m^2-sec */

                const real rate = rate_base * Wash_F * eta_dummy;

                R[1] = rate;
    }

    /* ======================= reaction-02 ======================= */
    {
        const real k2 = k_surface_covdep(A2_k, B2_beta, Ea2_Jpm, Tw, idx_site_r2, mu_r2, eps_r2, NS_R2, yi);
                const real Cg = gas_conc_cell(c0, t0, yi[IDX_C3H6], MW_C3H6);
                const real rate_base = k2 * Cg * Cv_R2;

                const real rate = rate_base * Wash_F * eta_dummy;

                R[2] = rate;
    }

    /* ======================= reaction-03 ======================= */
    {
        const real k3 = k_surface_covdep(A3_k, B3_beta, Ea3_Jpm, Tw, idx_site_r3, mu_r3, eps_r3, NS_R3, yi);
                const real Cg = gas_conc_cell(c0, t0, yi[IDX_C3H6], MW_C3H6);
                const real rate_base = k3 * Cg * Cv_R3;

                const real rate = rate_base * Wash_F * eta_dummy;

                R[3] = rate;
    }

    /* ======================= reaction-04 ======================= */
    {
        const real k4 = k_surface_covdep(A4_k, B4_beta, Ea4_Jpm, Tw, idx_site_r4, mu_r4, eps_r4, NS_R4, yi);
                const real Cg = gas_conc_cell(c0, t0, yi[IDX_H2], MW_H2);
                const real rate_base = k4 * Cg * Cv_R4;

                const real rate = rate_base * Wash_F * eta_dummy;

                R[4] = rate;
    }

    /* ======================= reaction-05 ======================= */
    {
        const real k5 = k_surface_covdep(A5_k, B5_beta, Ea5_Jpm, Tw, idx_site_r5, mu_r5, eps_r5, NS_R5, yi);
                const real Cg = gas_conc_cell(c0, t0, yi[IDX_H2O], MW_H2O);
                const real rate_base = k5 * Cg * Cv_R5;

                const real rate = rate_base * Wash_F * eta_dummy;

                R[5] = rate;
    }

    /* ======================= reaction-06 ======================= */
    {
        const real k6 = k_surface_covdep(A6_k, B6_beta, Ea6_Jpm, Tw, idx_site_r6, mu_r6, eps_r6, NS_R6, yi);
                const real Cg = gas_conc_cell(c0, t0, yi[IDX_CO2], MW_CO2);
                const real rate_base = k6 * Cg * Cv_R6;

                const real rate = rate_base * Wash_F * eta_dummy;

                R[6] = rate;
    }

    /* ======================= reaction-07 ======================= */
    {
        const real k7 = k_surface_covdep(A7_k, B7_beta, Ea7_Jpm, Tw, idx_site_r7, mu_r7, eps_r7, NS_R7, yi);
                const real Cg = gas_conc_cell(c0, t0, yi[IDX_CO], MW_CO);
                const real rate_base = k7 * Cg  * Cv_R7;

                R[7] = rate_base * Wash_F * eta_dummy;
    }

    /* ======================= reaction-08 ======================= */
    {
        const real k8 = k_surface_covdep(A8_k, B8_beta, Ea8_Jpm, Tw,
                                                 idx_site_r8, mu_r8, eps_r8, NS_R8, yi);

                const real theta_O = CLAMP01(yi[IDX_O_Pt]);
                const real rate_base = k8 * theta_O * theta_O * SITE_DEN_Pt * SITE_DEN_Pt; /* kmol/m^2-s */

                /* desorption: pore diffusion effectiveness = 1.0 (no gas reactant) */
                const real rate = rate_base * Wash_F * eta_dummy;

                R[8] = rate;
    }

    /* ======================= reaction-09 ======================= */
    {
        const real k9 = k_surface_covdep(A9_k, B9_beta, Ea9_Jpm, Tw,
                                                 idx_site_r9, mu_r9, eps_r9, NS_R9, yi);
                const real theta = CLAMP01(yi[IDX_C3H6_Pt]);
                const real rate_base = k9 * theta * SITE_DEN_Pt;
                const real rate = rate_base * Wash_F * eta_dummy;
                R[9] = rate;
    }

    /* ======================= reaction-10 ======================= */
    {
        const real k10 = k_surface_covdep(A10_k, B10_beta, Ea10_Jpm, Tw,
                                                  idx_site_r10, mu_r10, eps_r10, NS_R10, yi);
                const real thC3H5 = CLAMP01(yi[IDX_C3H5_Pt]);
                const real thOH = CLAMP01(yi[IDX_OH_Pt]);
                const real rate_base = k10 * thC3H5 * SITE_DEN_Pt * thOH * SITE_DEN_Pt;
                const real rate = rate_base * Wash_F * eta_dummy;
                R[10] = rate;
    }

    /* ======================= reaction-11 ======================= */
    {
        const real k11 = k_surface_covdep(A11_k, B11_beta, Ea11_Jpm, Tw,
                                                  idx_site_r11, mu_r11, eps_r11, NS_R11, yi);
                const real thH = CLAMP01(yi[IDX_H_Pt]);
                const real rate = (k11 * thH * thH * SITE_DEN_Pt * SITE_DEN_Pt) * Wash_F * eta_dummy;
                R[11] = rate;
    }

    /* ======================= reaction-12 ======================= */
    {
        const real k12 = k_surface_covdep(A12_k, B12_beta, Ea12_Jpm, Tw,
                                                  idx_site_r12, mu_r12, eps_r12, NS_R12, yi);
                const real thH2O = CLAMP01(yi[IDX_H2O_Pt]);
                const real rate = (k12 * thH2O * SITE_DEN_Pt) * Wash_F * eta_dummy;
                R[12] = rate;
    }

    /* ======================= reaction-13 ======================= */
    {
        const real k13 = k_surface_covdep(A13_k, B13_beta, Ea13_Jpm, Tw,
                                                  idx_site_r13, mu_r13, eps_r13, NS_R13, yi);
                const real thCO = CLAMP01(yi[IDX_CO_Pt]);
                const real rate_base = (k13 * thCO * SITE_DEN_Pt);


                R[13] = rate_base * Wash_F * eta_dummy;
    }

    /* ======================= reaction-14 ======================= */
    {
        const real k14 = k_surface_covdep(A14_k, B14_beta, Ea14_Jpm, Tw,
                                                  idx_site_r14, mu_r14, eps_r14, NS_R14, yi);
                const real thCO2 = CLAMP01(yi[IDX_CO2_Pt]);
                const real rate = (k14 * thCO2 * SITE_DEN_Pt) * Wash_F * eta_dummy;
                R[14] = rate;
    }

    /* ======================= reaction-15 ======================= */
    {
        const real k15 = k_surface_covdep(A15_k, B15_beta, Ea15_Jpm, Tw,
                                                  idx_site_r15, mu_r15, eps_r15, NS_R15, yi);
                const real thC3H5 = CLAMP01(yi[IDX_C3H5_Pt]);
                const real thO = CLAMP01(yi[IDX_O_Pt]);
                const real cC3H5 = thC3H5 * SITE_DEN_Pt;
                const real cO = thO * SITE_DEN_Pt;
                R[15] = (k15 * cC3H5 * cO * cO * cO * cO * cO) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-16 ======================= */
    {
        const real k16 = k_surface_covdep(A16_k, B16_beta, Ea16_Jpm, Tw,
                                                  idx_site_r16, mu_r16, eps_r16, NS_R16, yi);
                const real thC3H6 = CLAMP01(yi[IDX_C3H6_Pt]);
                const real thH = CLAMP01(yi[IDX_H_Pt]);
                R[16] = (k16 * thC3H6 * SITE_DEN_Pt * thH * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-17 ======================= */
    {
        const real k17 = k_surface_covdep(A17_k, B17_beta, Ea17_Jpm, Tw,
                                                  idx_site_r17, mu_r17, eps_r17, NS_R17, yi);
                const real thCC2H5 = CLAMP01(yi[IDX_CC2H5_Pt]);
                const real thH = CLAMP01(yi[IDX_H_Pt]);
                R[17] = (k17 * thCC2H5 * SITE_DEN_Pt * thH * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-18 ======================= */
    {
        const real k18 = k_surface_covdep(A18_k, B18_beta, Ea18_Jpm, Tw,
                                                  idx_site_r18, mu_r18, eps_r18, NS_R18, yi);
                const real thCC2H5 = CLAMP01(yi[IDX_CC2H5_Pt]);
                const real thPt = CLAMP01(yi[IDX_Pt_Vac]);
                R[18] = (k18 * thCC2H5 * SITE_DEN_Pt * thPt * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-19 ======================= */
    {
        const real k19 = k_surface_covdep(A19_k, B19_beta, Ea19_Jpm, Tw,
                                                  idx_site_r19, mu_r19, eps_r19, NS_R19, yi);
                const real thC2H3 = CLAMP01(yi[IDX_C2H3_Pt]);
                const real thCH2 = CLAMP01(yi[IDX_CH2_Pt]);
                R[19] = (k19 * thC2H3 *  SITE_DEN_Pt * thCH2 *  SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-20 ======================= */
    {
        const real k20 = k_surface_covdep(A20_k, B20_beta, Ea20_Jpm, Tw,
                                                  idx_site_r20, mu_r20, eps_r20, NS_R20, yi);
                const real thC2H3 = CLAMP01(yi[IDX_C2H3_Pt]);
                const real thPt = CLAMP01(yi[IDX_Pt_Vac]);
                R[20] = (k20 * thC2H3 * SITE_DEN_Pt * thPt * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-21 ======================= */
    {
        const real k21 = k_surface_covdep(A21_k, B21_beta, Ea21_Jpm, Tw,
                                                  idx_site_r21, mu_r21, eps_r21, NS_R21, yi);
                const real thCH3 = CLAMP01(yi[IDX_CH3_Pt]);
                const real thC = CLAMP01(yi[IDX_C_Pt]);
                R[21] = (k21 * thCH3 * SITE_DEN_Pt * thC * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-22 ======================= */
    {
        const real k22 = k_surface_covdep(A22_k, B22_beta, Ea22_Jpm, Tw,
                                                  idx_site_r22, mu_r22, eps_r22, NS_R22, yi);
                const real thCH3 = CLAMP01(yi[IDX_CH3_Pt]);
                const real thPt = CLAMP01(yi[IDX_Pt_Vac]);
                R[22] = (k22 * thCH3 * SITE_DEN_Pt * thPt * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-23 ======================= */
    {
        const real k23 = k_surface_covdep(A23_k, B23_beta, Ea23_Jpm, Tw,
                                                  idx_site_r23, mu_r23, eps_r23, NS_R23, yi);
                const real thCH2 = CLAMP01(yi[IDX_CH2_Pt]);
                const real thH = CLAMP01(yi[IDX_H_Pt]);
                R[23] = (k23 * SITE_DEN_Pt * thCH2 * thH * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-24 ======================= */
    {
        const real k24 = k_surface_covdep(A24_k, B24_beta, Ea24_Jpm, Tw,
                                                  idx_site_r24, mu_r24, eps_r24, NS_R24, yi);
                const real thCH2 = CLAMP01(yi[IDX_CH2_Pt]);
                const real thPt = CLAMP01(yi[IDX_Pt_Vac]);
                R[24] = (k24 * thCH2 * thPt * SITE_DEN_Pt * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-25 ======================= */
    {
        const real k25 = k_surface_covdep(A25_k, B25_beta, Ea25_Jpm, Tw,
                                                  idx_site_r25, mu_r25, eps_r25, NS_R25, yi);
                const real thCH = CLAMP01(yi[IDX_CH_Pt]);
                const real thH = CLAMP01(yi[IDX_H_Pt]);
                R[25] = (k25 * thCH * thH * SITE_DEN_Pt * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-26 ======================= */
    {
        const real k26 = k_surface_covdep(A26_k, B26_beta, Ea26_Jpm, Tw,
                                                  idx_site_r26, mu_r26, eps_r26, NS_R26, yi);
                const real thCH = CLAMP01(yi[IDX_CH_Pt]);
                const real thPt = CLAMP01(yi[IDX_Pt_Vac]);
                R[26] = (k26 * thCH * thPt * SITE_DEN_Pt * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-27 ======================= */
    {
        const real k27 = k_surface_covdep(A27_k, B27_beta, Ea27_Jpm, Tw,
                                                  idx_site_r27, mu_r27, eps_r27, NS_R27, yi);
                const real thC = CLAMP01(yi[IDX_C_Pt]);
                const real thH = CLAMP01(yi[IDX_H_Pt]);
                R[27] = (k27 * thC * thH * SITE_DEN_Pt * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-28 ======================= */
    {
        const real k28 = k_surface_covdep(A28_k, B28_beta, Ea28_Jpm, Tw,
                                                  idx_site_r28, mu_r28, eps_r28, NS_R28, yi);
                const real thC2H3 = CLAMP01(yi[IDX_C2H3_Pt]);
                const real thO = CLAMP01(yi[IDX_O_Pt]);
                R[28] = (k28 * thC2H3 * SITE_DEN_Pt * thO * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-29 ======================= */
    {
        const real k29 = k_surface_covdep(A29_k, B29_beta, Ea29_Jpm, Tw,
                                                  idx_site_r29, mu_r29, eps_r29, NS_R29, yi);
                const real thCH3CO = CLAMP01(yi[IDX_CH3CO_Pt]);
                const real thPt = CLAMP01(yi[IDX_Pt_Vac]);
                R[29] = (k29 * thCH3CO * SITE_DEN_Pt * thPt * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-30 ======================= */
    {
        const real k30 = k_surface_covdep(A30_k, B30_beta, Ea30_Jpm, Tw,
                                                  idx_site_r30, mu_r30, eps_r30, NS_R30, yi);
                const real thCH3 = CLAMP01(yi[IDX_CH3_Pt]);
                const real thCO = CLAMP01(yi[IDX_CO_Pt]);
                R[30] = (k30 * thCH3* SITE_DEN_Pt * thCO * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-31 ======================= */
    {
        const real k31 = k_surface_covdep(A31_k, B31_beta, Ea31_Jpm, Tw,
                                                  idx_site_r31, mu_r31, eps_r31, NS_R31, yi);
                const real thCH3CO = CLAMP01(yi[IDX_CH3CO_Pt]);
                const real thPt = CLAMP01(yi[IDX_Pt_Vac]);
                R[31] = (k31 * thCH3CO * SITE_DEN_Pt * thPt * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-32 ======================= */
    {
        const real k32 = k_surface_covdep(A32_k, B32_beta, Ea32_Jpm, Tw,
                                                  idx_site_r32, mu_r32, eps_r32, NS_R32, yi);
                const real thCH3 = CLAMP01(yi[IDX_CH3_Pt]);
                const real thO = CLAMP01(yi[IDX_O_Pt]);
                R[32] = (k32 * thCH3 * SITE_DEN_Pt * thO * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-33 ======================= */
    {
        const real k33 = k_surface_covdep(A33_k, B33_beta, Ea33_Jpm, Tw,
                                                  idx_site_r33, mu_r33, eps_r33, NS_R33, yi);
                const real thCH2 = CLAMP01(yi[IDX_CH2_Pt]);
                const real thOH = CLAMP01(yi[IDX_OH_Pt]);
                R[33] = (k33 * thCH2 * SITE_DEN_Pt * thOH * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-34 ======================= */
    {
        const real k34 = k_surface_covdep(A34_k, B34_beta, Ea34_Jpm, Tw,
                                                  idx_site_r34, mu_r34, eps_r34, NS_R34, yi);
                const real thCH2 = CLAMP01(yi[IDX_CH2_Pt]);
                const real thO = CLAMP01(yi[IDX_O_Pt]);
                R[34] = (k34 * thCH2 * SITE_DEN_Pt * thO * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-35 ======================= */
    {
        const real k35 = k_surface_covdep(A35_k, B35_beta, Ea35_Jpm, Tw,
                                                  idx_site_r35, mu_r35, eps_r35, NS_R35, yi);
                const real thCH = CLAMP01(yi[IDX_CH_Pt]);
                const real thOH = CLAMP01(yi[IDX_OH_Pt]);
                R[35] = (k35 * thCH * SITE_DEN_Pt * thOH * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-36 ======================= */
    {
        const real k36 = k_surface_covdep(A36_k, B36_beta, Ea36_Jpm, Tw,
                                                  idx_site_r36, mu_r36, eps_r36, NS_R36, yi);
                const real thCH = CLAMP01(yi[IDX_CH_Pt]);
                const real thO = CLAMP01(yi[IDX_O_Pt]);
                R[36] = (k36 * thCH * SITE_DEN_Pt * thO * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-37 ======================= */
    {
        const real k37 = k_surface_covdep(A37_k, B37_beta, Ea37_Jpm, Tw,
                                                  idx_site_r37, mu_r37, eps_r37, NS_R37, yi);
                const real thC = CLAMP01(yi[IDX_C_Pt]);
                const real thOH = CLAMP01(yi[IDX_OH_Pt]);
                R[37] = (k37 * thC * SITE_DEN_Pt * thOH * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-38 ======================= */
    {
        const real k38 = k_surface_covdep(A38_k, B38_beta, Ea38_Jpm, Tw,
                                                  idx_site_r38, mu_r38, eps_r38, NS_R38, yi);
                const real thO = CLAMP01(yi[IDX_O_Pt]);
                const real thH = CLAMP01(yi[IDX_H_Pt]);
                R[38] = (k38 * thO * SITE_DEN_Pt * thH * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-39 ======================= */
    {
        const real k39 = k_surface_covdep(A39_k, B39_beta, Ea39_Jpm, Tw,
                                                  idx_site_r39, mu_r39, eps_r39, NS_R39, yi);
                const real thOH = CLAMP01(yi[IDX_OH_Pt]);
                const real thPt = CLAMP01(yi[IDX_Pt_Vac]);
                R[39] = (k39 * thOH * SITE_DEN_Pt * thPt * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-40 ======================= */
    {
        const real k40 = k_surface_covdep(A40_k, B40_beta, Ea40_Jpm, Tw,
                                                  idx_site_r40, mu_r40, eps_r40, NS_R40, yi);
                const real thH = CLAMP01(yi[IDX_H_Pt]);
                const real thOH = CLAMP01(yi[IDX_OH_Pt]);
                R[40] = (k40 * thH * SITE_DEN_Pt * thOH * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-41 ======================= */
    {
        const real k41 = k_surface_covdep(A41_k, B41_beta, Ea41_Jpm, Tw,
                                                  idx_site_r41, mu_r41, eps_r41, NS_R41, yi);
                const real thH2O = CLAMP01(yi[IDX_H2O_Pt]);
                const real thPt = CLAMP01(yi[IDX_Pt_Vac]);
                R[41] = (k41 * thH2O * SITE_DEN_Pt * thPt * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-42 ======================= */
    {
        const real k42 = k_surface_covdep(A42_k, B42_beta, Ea42_Jpm, Tw,
                                                  idx_site_r42, mu_r42, eps_r42, NS_R42, yi);
                const real thOH = CLAMP01(yi[IDX_OH_Pt]);
                R[42] = (k42 * thOH * thOH * SITE_DEN_Pt * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-43 ======================= */
    {
        const real k43 = k_surface_covdep(A43_k, B43_beta, Ea43_Jpm, Tw,
                                                  idx_site_r43, mu_r43, eps_r43, NS_R43, yi);
                const real thH2O = CLAMP01(yi[IDX_H2O_Pt]);
                const real thO = CLAMP01(yi[IDX_O_Pt]);
                R[43] = (k43 * thH2O * SITE_DEN_Pt * thO * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-44 ======================= */
    {
        const real k44 = k_surface_covdep(A44_k, B44_beta, Ea44_Jpm, Tw,
                                                  idx_site_r44, mu_r44, eps_r44, NS_R44, yi);

                const real theta_CO = CLAMP01(yi[IDX_CO_Pt]);
                const real theta_O  = CLAMP01(yi[IDX_O_Pt]);

                const real rate_base = k44 * theta_CO * SITE_DEN_Pt * theta_O * SITE_DEN_Pt; /* kmol/m^2-s */
                const real rate      = rate_base * Wash_F * eta_dummy;

                R[44] = rate;
    }

    /* ======================= reaction-45 ======================= */
    {
        const real k45 = k_surface_covdep(A45_k, B45_beta, Ea45_Jpm, Tw,
                                                  idx_site_r45, mu_r45, eps_r45, NS_R45, yi);
                const real thCO2 = CLAMP01(yi[IDX_CO2_Pt]);
                const real thPt = CLAMP01(yi[IDX_Pt_Vac]);
                R[45] = (k45 * thCO2 * SITE_DEN_Pt * thPt * SITE_DEN_Pt) * Wash_F * eta_dummy;
    }

    /* ======================= reaction-46 ======================= */
    {
        const real k46 = k_surface_covdep(A46_k, B46_beta, Ea46_Jpm, Tw,
                                                  idx_site_r46, mu_r46, eps_r46, NS_R46, yi);
                const real thC = CLAMP01(yi[IDX_C_Pt]);
                const real thO = CLAMP01(yi[IDX_O_Pt]);
                const real rate = k46 * thC * thO * SITE_DEN_Pt * SITE_DEN_Pt * Wash_F * eta_dummy;
                R[46] = rate;
    }

    /* ======================= reaction-47 ======================= */
    {
        const real k47 = k_surface_covdep(A47_k, B47_beta, Ea47_Jpm, Tw,
                                                  idx_site_r47, mu_r47, eps_r47, NS_R47, yi);
                const real thCO = CLAMP01(yi[IDX_CO_Pt]);
                const real thPt = CLAMP01(yi[IDX_Pt_Vac]);
                const real rate = k47 * thCO * SITE_DEN_Pt * thPt * SITE_DEN_Pt * Wash_F * eta_dummy;
                R[47] = rate;
    }

    /* ======================= reaction-48 ======================= */
    {
        const real k48 = k_surface_covdep(A48_k, B48_beta, Ea48_Jpm, Tw, idx_site_r48, mu_r48, eps_r48, NS_R48, yi);
                const real Cg = gas_conc_cell(c0, t0, yi[IDX_NO], MW_NO);
                const real rate_base = k48 * Cg * Cv_R48;

                const real rate = rate_base * Wash_F * eta_dummy;

                R[48] = rate;
    }

    /* ======================= reaction-49 ======================= */
    {
        const real k49 = k_surface_covdep(A49_k, B49_beta, Ea49_Jpm, Tw,
                                                  idx_site_r49, mu_r49, eps_r49, NS_R49, yi);
                const real thNO = CLAMP01(yi[IDX_NO_Pt]);
                const real rate = k49 * thNO * SITE_DEN_Pt * Wash_F * eta_dummy; 
                R[49] = rate;
    }

    /* ======================= reaction-50 ======================= */
    {
        const real k50 = k_surface_covdep(A50_k, B50_beta, Ea50_Jpm, Tw,
                                                  idx_site_r50, mu_r50, eps_r50, NS_R50, yi);
                const real thN = CLAMP01(yi[IDX_N_Pt]);
                const real rate = k50 * thN * thN * SITE_DEN_Pt * SITE_DEN_Pt * Wash_F * eta_dummy; 
                R[50] = rate;
    }

    /* ======================= reaction-51 ======================= */
    {
        const real k51 = k_surface_covdep(A51_k, B51_beta, Ea51_Jpm, Tw,
                                                  idx_site_r51, mu_r51, eps_r51, NS_R51, yi);
                const real thNO  = CLAMP01(yi[IDX_NO_Pt]);
                const real thVac = CLAMP01(yi[IDX_Pt_Vac]);
                const real rate = k51 * thNO * thVac * SITE_DEN_Pt * SITE_DEN_Pt * Wash_F * eta_dummy; 
                R[51] = rate;
    }

    /* ======================= reaction-52 ======================= */
    {
        const real k52 = k_surface_covdep(A52_k, B52_beta, Ea52_Jpm, Tw,
                                                  idx_site_r52, mu_r52, eps_r52, NS_R52, yi);
                const real thN = CLAMP01(yi[IDX_N_Pt]);
                const real thO = CLAMP01(yi[IDX_O_Pt]);
                const real rate = k52 * thN * thO * SITE_DEN_Pt * SITE_DEN_Pt * Wash_F * eta_dummy;
                R[52] = rate;
    }

    /* ======================= reaction-53 ======================= */
    {
        const real k53 = k_surface_covdep(A53_k, B53_beta, Ea53_Jpm, Tw, idx_site_r53, mu_r53, eps_r53, NS_R53, yi);
                const real Cg = gas_conc_cell(c0, t0, yi[IDX_O2], MW_O2);
                const real rate_base = k53 * Cg * Cv_R53;

                const real rate = rate_base * Wash_F * eta_dummy;

                R[53] = rate;
    }

    /* ======================= reaction-54 ======================= */
    {
        const real k54 = k_surface_covdep(A54_k, B54_beta, Ea54_Jpm, Tw, idx_site_r54, mu_r54, eps_r54, NS_R54, yi);
                const real Cg = gas_conc_cell(c0, t0, yi[IDX_CO], MW_CO);
                const real rate_base = k54 * Cg  * Cv_R54;

                R[54] = rate_base * Wash_F * eta_dummy;
    }

    /* ======================= reaction-55 ======================= */
    {
        const real k55 = k_surface_covdep(A55_k, B55_beta, Ea55_Jpm, Tw, idx_site_r55, mu_r55, eps_r55, NS_R55, yi);
                const real Cg = gas_conc_cell(c0, t0, yi[IDX_NO], MW_NO);
                const real rate_base = k55 * Cg * Cv_R55;

                const real rate = rate_base * Wash_F * eta_dummy;

                R[55] = rate;
    }

    /* ======================= reaction-56 ======================= */
    {
        const real k56 = k_surface_covdep(A56_k, B56_beta, Ea56_Jpm, Tw,
                                                  idx_site_r56, mu_r56, eps_r56, NS_R56, yi);
                const real thO_rh = CLAMP01(yi[IDX_O_Rh]);
                const real rate = k56 * thO_rh * thO_rh * SITE_DEN_Rh * SITE_DEN_Rh * Wash_F * eta_dummy; 
                R[56] = rate;
    }

    /* ======================= reaction-57 ======================= */
    {
        const real k57 = k_surface_covdep(A57_k, B57_beta, Ea57_Jpm, Tw,
                                                  idx_site_r57, mu_r57, eps_r57, NS_R57, yi);
                const real thCO_rh = CLAMP01(yi[IDX_CO_Rh]);
                const real rate_base = k57 * thCO_rh * SITE_DEN_Rh;

                R[57] = rate_base * Wash_F * eta_dummy;
    }

    /* ======================= reaction-58 ======================= */
    {
        const real k58 = k_surface_covdep(A58_k, B58_beta, Ea58_Jpm, Tw,
                                                  idx_site_r58, mu_r58, eps_r58, NS_R58, yi);
                const real thNO_rh = CLAMP01(yi[IDX_NO_Rh]);
                const real rate = k58 * thNO_rh * SITE_DEN_Rh * Wash_F * eta_dummy; 
                R[58] = rate;
    }

    /* ======================= reaction-59 ======================= */
    {
        const real k59 = k_surface_covdep(A59_k, B59_beta, Ea59_Jpm, Tw,
                                                  idx_site_r59, mu_r59, eps_r59, NS_R59, yi);
                const real thN_rh = CLAMP01(yi[IDX_N_Rh]);
                const real rate = k59 * thN_rh * thN_rh * SITE_DEN_Rh * SITE_DEN_Rh * Wash_F * eta_dummy; 
                R[59] = rate;
    }

    /* ======================= reaction-60 ======================= */
    {
        const real k60 = k_surface_covdep(A60_k, B60_beta, Ea60_Jpm, Tw,
                                                  idx_site_r60, mu_r60, eps_r60, NS_R60, yi);
                const real thCO_rh = CLAMP01(yi[IDX_CO_Rh]);
                const real thO_rh  = CLAMP01(yi[IDX_O_Rh]);
                const real rate = k60 * thCO_rh * thO_rh * SITE_DEN_Rh * SITE_DEN_Rh * Wash_F * eta_dummy;
                R[60] = rate;
    }

    /* ======================= reaction-61 ======================= */
    {
        const real k61 = k_surface_covdep(A61_k, B61_beta, Ea61_Jpm, Tw,
                                                  idx_site_r61, mu_r61, eps_r61, NS_R61, yi);
                const real thNO_rh  = CLAMP01(yi[IDX_NO_Rh]);
                const real thVac_rh = CLAMP01(yi[IDX_Rh_Vac]);
                const real rate = k61 * thNO_rh * thVac_rh * SITE_DEN_Rh * SITE_DEN_Rh * Wash_F * eta_dummy; 
                R[61] = rate;
    }

}

DEFINE_NET_REACTION_RATE(chatterjee_pt_ads_net, c, t, particle, pressure, temp, yi, rr, jac)
{
    int i;
    int ns = n_spe;

    (void)particle;
    (void)pressure;

    if (!temp || !yi || !rr) return;
    if (ns <= IDX_N_Rh) return;

    for (i = 0; i < ns; ++i) rr[i] = 0.0;
    if (jac) {
        int nn = ns * ns;
        for (i = 0; i < nn; ++i) jac[i] = 0.0;
    }

    {
        const cell_t c0 = c;
        Thread *t0 = t;
        const real Tw = *temp;
        real R[62];

        chatterjee_compute_R(c0, t0, Tw, yi, R);

    rr[IDX_O2] += -1.0 * R[1];  /* R1: O2 */
    rr[IDX_Pt_Vac] += -2.0 * R[1];  /* R1: Pt(S) */
    rr[IDX_O_Pt] += 2.0 * R[1];  /* R1: O(S) */
    rr[IDX_C3H6] += -1.0 * R[2];  /* R2: C3H6 */
    rr[IDX_Pt_Vac] += -2.0 * R[2];  /* R2: Pt(S) */
    rr[IDX_C3H6_Pt] += 1.0 * R[2];  /* R2: C3H6(S) */
    rr[IDX_C3H6] += -1.0 * R[3];  /* R3: C3H6 */
    rr[IDX_O_Pt] += -1.0 * R[3];  /* R3: O(S) */
    rr[IDX_Pt_Vac] += -1.0 * R[3];  /* R3: Pt(S) */
    rr[IDX_C3H5_Pt] += 1.0 * R[3];  /* R3: C3H5(S) */
    rr[IDX_OH_Pt] += 1.0 * R[3];  /* R3: OH(S) */
    rr[IDX_H2] += -1.0 * R[4];  /* R4: H2 */
    rr[IDX_Pt_Vac] += -2.0 * R[4];  /* R4: Pt(S) */
    rr[IDX_H_Pt] += 2.0 * R[4];  /* R4: H(S) */
    rr[IDX_H2O] += -1.0 * R[5];  /* R5: H2O */
    rr[IDX_Pt_Vac] += -1.0 * R[5];  /* R5: Pt(S) */
    rr[IDX_H2O_Pt] += 1.0 * R[5];  /* R5: H2O(S) */
    rr[IDX_CO2] += -1.0 * R[6];  /* R6: CO2 */
    rr[IDX_Pt_Vac] += -1.0 * R[6];  /* R6: Pt(S) */
    rr[IDX_CO2_Pt] += 1.0 * R[6];  /* R6: CO2(S) */
    rr[IDX_CO] += -1.0 * R[7];  /* R7: CO */
    rr[IDX_Pt_Vac] += -1.0 * R[7];  /* R7: Pt(S) */
    rr[IDX_CO_Pt] += 1.0 * R[7];  /* R7: CO(S) */
    rr[IDX_O_Pt] += -2.0 * R[8];  /* R8: O(S) */
    rr[IDX_Pt_Vac] += 2.0 * R[8];  /* R8: Pt(S) */
    rr[IDX_O2] += 1.0 * R[8];  /* R8: O2 */
    rr[IDX_C3H6_Pt] += -1.0 * R[9];  /* R9: C3H6(S) */
    rr[IDX_C3H6] += 1.0 * R[9];  /* R9: C3H6 */
    rr[IDX_Pt_Vac] += 2.0 * R[9];  /* R9: Pt(S) */
    rr[IDX_C3H5_Pt] += -1.0 * R[10];  /* R10: C3H5(S) */
    rr[IDX_OH_Pt] += -1.0 * R[10];  /* R10: OH(S) */
    rr[IDX_C3H6] += 1.0 * R[10];  /* R10: C3H6 */
    rr[IDX_O_Pt] += 1.0 * R[10];  /* R10: O(S) */
    rr[IDX_Pt_Vac] += 1.0 * R[10];  /* R10: Pt(S) */
    rr[IDX_H_Pt] += -2.0 * R[11];  /* R11: H(S) */
    rr[IDX_H2] += 1.0 * R[11];  /* R11: H2 */
    rr[IDX_Pt_Vac] += 2.0 * R[11];  /* R11: Pt(S) */
    rr[IDX_H2O_Pt] += -1.0 * R[12];  /* R12: H2O(S) */
    rr[IDX_Pt_Vac] += 1.0 * R[12];  /* R12: Pt(S) */
    rr[IDX_H2O] += 1.0 * R[12];  /* R12: H2O */
    rr[IDX_CO_Pt] += -1.0 * R[13];  /* R13: CO(S) */
    rr[IDX_CO] += 1.0 * R[13];  /* R13: CO */
    rr[IDX_Pt_Vac] += 1.0 * R[13];  /* R13: Pt(S) */
    rr[IDX_CO2_Pt] += -1.0 * R[14];  /* R14: CO2(S) */
    rr[IDX_CO2] += 1.0 * R[14];  /* R14: CO2 */
    rr[IDX_Pt_Vac] += 1.0 * R[14];  /* R14: Pt(S) */
    rr[IDX_C3H5_Pt] += -1.0 * R[15];  /* R15: C3H5(S) */
    rr[IDX_O_Pt] += -5.0 * R[15];  /* R15: O(S) */
    rr[IDX_Pt_Vac] += -2.0 * R[15];  /* R15: Pt(S) */
    rr[IDX_OH_Pt] += 5.0 * R[15];  /* R15: OH(S) */
    rr[IDX_C_Pt] += 3.0 * R[15];  /* R15: C(S) */
    rr[IDX_C3H6_Pt] += -1.0 * R[16];  /* R16: C3H6(S) */
    rr[IDX_CC2H5_Pt] += 1.0 * R[16];  /* R16: CC2H5(S) */
    rr[IDX_H_Pt] += 1.0 * R[16];  /* R16: H(S) */
    rr[IDX_CC2H5_Pt] += -1.0 * R[17];  /* R17: CC2H5(S) */
    rr[IDX_H_Pt] += -1.0 * R[17];  /* R17: H(S) */
    rr[IDX_C3H6_Pt] += 1.0 * R[17];  /* R17: C3H6(S) */
    rr[IDX_CC2H5_Pt] += -1.0 * R[18];  /* R18: CC2H5(S) */
    rr[IDX_Pt_Vac] += -1.0 * R[18];  /* R18: Pt(S) */
    rr[IDX_C2H3_Pt] += 1.0 * R[18];  /* R18: C2H3(S) */
    rr[IDX_CH2_Pt] += 1.0 * R[18];  /* R18: CH2(S) */
    rr[IDX_C2H3_Pt] += -1.0 * R[19];  /* R19: C2H3(S) */
    rr[IDX_CH2_Pt] += -1.0 * R[19];  /* R19: CH2(S) */
    rr[IDX_CC2H5_Pt] += 1.0 * R[19];  /* R19: CC2H5(S) */
    rr[IDX_Pt_Vac] += 1.0 * R[19];  /* R19: Pt(S) */
    rr[IDX_C2H3_Pt] += -1.0 * R[20];  /* R20: C2H3(S) */
    rr[IDX_Pt_Vac] += -1.0 * R[20];  /* R20: Pt(S) */
    rr[IDX_CH3_Pt] += 1.0 * R[20];  /* R20: CH3(S) */
    rr[IDX_C_Pt] += 1.0 * R[20];  /* R20: C(S) */
    rr[IDX_CH3_Pt] += -1.0 * R[21];  /* R21: CH3(S) */
    rr[IDX_C_Pt] += -1.0 * R[21];  /* R21: C(S) */
    rr[IDX_C2H3_Pt] += 1.0 * R[21];  /* R21: C2H3(S) */
    rr[IDX_Pt_Vac] += 1.0 * R[21];  /* R21: Pt(S) */
    rr[IDX_CH3_Pt] += -1.0 * R[22];  /* R22: CH3(S) */
    rr[IDX_Pt_Vac] += -1.0 * R[22];  /* R22: Pt(S) */
    rr[IDX_CH2_Pt] += 1.0 * R[22];  /* R22: CH2(S) */
    rr[IDX_H_Pt] += 1.0 * R[22];  /* R22: H(S) */
    rr[IDX_CH2_Pt] += -1.0 * R[23];  /* R23: CH2(S) */
    rr[IDX_H_Pt] += -1.0 * R[23];  /* R23: H(S) */
    rr[IDX_CH3_Pt] += 1.0 * R[23];  /* R23: CH3(S) */
    rr[IDX_Pt_Vac] += 1.0 * R[23];  /* R23: Pt(S) */
    rr[IDX_CH2_Pt] += -1.0 * R[24];  /* R24: CH2(S) */
    rr[IDX_Pt_Vac] += -1.0 * R[24];  /* R24: Pt(S) */
    rr[IDX_CH_Pt] += 1.0 * R[24];  /* R24: CH(S) */
    rr[IDX_H_Pt] += 1.0 * R[24];  /* R24: H(S) */
    rr[IDX_CH_Pt] += -1.0 * R[25];  /* R25: CH(S) */
    rr[IDX_H_Pt] += -1.0 * R[25];  /* R25: H(S) */
    rr[IDX_CH2_Pt] += 1.0 * R[25];  /* R25: CH2(S) */
    rr[IDX_Pt_Vac] += 1.0 * R[25];  /* R25: Pt(S) */
    rr[IDX_CH_Pt] += -1.0 * R[26];  /* R26: CH(S) */
    rr[IDX_Pt_Vac] += -1.0 * R[26];  /* R26: Pt(S) */
    rr[IDX_C_Pt] += 1.0 * R[26];  /* R26: C(S) */
    rr[IDX_H_Pt] += 1.0 * R[26];  /* R26: H(S) */
    rr[IDX_C_Pt] += -1.0 * R[27];  /* R27: C(S) */
    rr[IDX_H_Pt] += -1.0 * R[27];  /* R27: H(S) */
    rr[IDX_CH_Pt] += 1.0 * R[27];  /* R27: CH(S) */
    rr[IDX_Pt_Vac] += 1.0 * R[27];  /* R27: Pt(S) */
    rr[IDX_C2H3_Pt] += -1.0 * R[28];  /* R28: C2H3(S) */
    rr[IDX_O_Pt] += -1.0 * R[28];  /* R28: O(S) */
    rr[IDX_Pt_Vac] += 1.0 * R[28];  /* R28: Pt(S) */
    rr[IDX_CH3CO_Pt] += 1.0 * R[28];  /* R28: CH3CO(S) */
    rr[IDX_CH3CO_Pt] += -1.0 * R[29];  /* R29: CH3CO(S) */
    rr[IDX_Pt_Vac] += -1.0 * R[29];  /* R29: Pt(S) */
    rr[IDX_C2H3_Pt] += 1.0 * R[29];  /* R29: C2H3(S) */
    rr[IDX_O_Pt] += 1.0 * R[29];  /* R29: O(S) */
    rr[IDX_CH3_Pt] += -1.0 * R[30];  /* R30: CH3(S) */
    rr[IDX_CO_Pt] += -1.0 * R[30];  /* R30: CO(S) */
    rr[IDX_Pt_Vac] += 1.0 * R[30];  /* R30: Pt(S) */
    rr[IDX_CH3CO_Pt] += 1.0 * R[30];  /* R30: CH3CO(S) */
    rr[IDX_CH3CO_Pt] += -1.0 * R[31];  /* R31: CH3CO(S) */
    rr[IDX_Pt_Vac] += -1.0 * R[31];  /* R31: Pt(S) */
    rr[IDX_CH3_Pt] += 1.0 * R[31];  /* R31: CH3(S) */
    rr[IDX_CO_Pt] += 1.0 * R[31];  /* R31: CO(S) */
    rr[IDX_CH3_Pt] += -1.0 * R[32];  /* R32: CH3(S) */
    rr[IDX_O_Pt] += -1.0 * R[32];  /* R32: O(S) */
    rr[IDX_CH2_Pt] += 1.0 * R[32];  /* R32: CH2(S) */
    rr[IDX_OH_Pt] += 1.0 * R[32];  /* R32: OH(S) */
    rr[IDX_CH2_Pt] += -1.0 * R[33];  /* R33: CH2(S) */
    rr[IDX_OH_Pt] += -1.0 * R[33];  /* R33: OH(S) */
    rr[IDX_CH3_Pt] += 1.0 * R[33];  /* R33: CH3(S) */
    rr[IDX_O_Pt] += 1.0 * R[33];  /* R33: O(S) */
    rr[IDX_CH2_Pt] += -1.0 * R[34];  /* R34: CH2(S) */
    rr[IDX_O_Pt] += -1.0 * R[34];  /* R34: O(S) */
    rr[IDX_CH_Pt] += 1.0 * R[34];  /* R34: CH(S) */
    rr[IDX_OH_Pt] += 1.0 * R[34];  /* R34: OH(S) */
    rr[IDX_CH_Pt] += -1.0 * R[35];  /* R35: CH(S) */
    rr[IDX_OH_Pt] += -1.0 * R[35];  /* R35: OH(S) */
    rr[IDX_CH2_Pt] += 1.0 * R[35];  /* R35: CH2(S) */
    rr[IDX_O_Pt] += 1.0 * R[35];  /* R35: O(S) */
    rr[IDX_CH_Pt] += -1.0 * R[36];  /* R36: CH(S) */
    rr[IDX_O_Pt] += -1.0 * R[36];  /* R36: O(S) */
    rr[IDX_C_Pt] += 1.0 * R[36];  /* R36: C(S) */
    rr[IDX_OH_Pt] += 1.0 * R[36];  /* R36: OH(S) */
    rr[IDX_C_Pt] += -1.0 * R[37];  /* R37: C(S) */
    rr[IDX_OH_Pt] += -1.0 * R[37];  /* R37: OH(S) */
    rr[IDX_CH_Pt] += 1.0 * R[37];  /* R37: CH(S) */
    rr[IDX_O_Pt] += 1.0 * R[37];  /* R37: O(S) */
    rr[IDX_O_Pt] += -1.0 * R[38];  /* R38: O(S) */
    rr[IDX_H_Pt] += -1.0 * R[38];  /* R38: H(S) */
    rr[IDX_OH_Pt] += 1.0 * R[38];  /* R38: OH(S) */
    rr[IDX_Pt_Vac] += 1.0 * R[38];  /* R38: Pt(S) */
    rr[IDX_OH_Pt] += -1.0 * R[39];  /* R39: OH(S) */
    rr[IDX_Pt_Vac] += -1.0 * R[39];  /* R39: Pt(S) */
    rr[IDX_O_Pt] += 1.0 * R[39];  /* R39: O(S) */
    rr[IDX_H_Pt] += 1.0 * R[39];  /* R39: H(S) */
    rr[IDX_H_Pt] += -1.0 * R[40];  /* R40: H(S) */
    rr[IDX_OH_Pt] += -1.0 * R[40];  /* R40: OH(S) */
    rr[IDX_H2O_Pt] += 1.0 * R[40];  /* R40: H2O(S) */
    rr[IDX_Pt_Vac] += 1.0 * R[40];  /* R40: Pt(S) */
    rr[IDX_H2O_Pt] += -1.0 * R[41];  /* R41: H2O(S) */
    rr[IDX_Pt_Vac] += -1.0 * R[41];  /* R41: Pt(S) */
    rr[IDX_H_Pt] += 1.0 * R[41];  /* R41: H(S) */
    rr[IDX_OH_Pt] += 1.0 * R[41];  /* R41: OH(S) */
    rr[IDX_OH_Pt] += -2.0 * R[42];  /* R42: OH(S) */
    rr[IDX_H2O_Pt] += 1.0 * R[42];  /* R42: H2O(S) */
    rr[IDX_O_Pt] += 1.0 * R[42];  /* R42: O(S) */
    rr[IDX_H2O_Pt] += -1.0 * R[43];  /* R43: H2O(S) */
    rr[IDX_O_Pt] += -1.0 * R[43];  /* R43: O(S) */
    rr[IDX_OH_Pt] += 2.0 * R[43];  /* R43: OH(S) */
    rr[IDX_CO_Pt] += -1.0 * R[44];  /* R44: CO(S) */
    rr[IDX_O_Pt] += -1.0 * R[44];  /* R44: O(S) */
    rr[IDX_CO2_Pt] += 1.0 * R[44];  /* R44: CO2(S) */
    rr[IDX_Pt_Vac] += 1.0 * R[44];  /* R44: Pt(S) */
    rr[IDX_CO2_Pt] += -1.0 * R[45];  /* R45: CO2(S) */
    rr[IDX_Pt_Vac] += -1.0 * R[45];  /* R45: Pt(S) */
    rr[IDX_CO_Pt] += 1.0 * R[45];  /* R45: CO(S) */
    rr[IDX_O_Pt] += 1.0 * R[45];  /* R45: O(S) */
    rr[IDX_C_Pt] += -1.0 * R[46];  /* R46: C(S) */
    rr[IDX_O_Pt] += -1.0 * R[46];  /* R46: O(S) */
    rr[IDX_CO_Pt] += 1.0 * R[46];  /* R46: CO(S) */
    rr[IDX_Pt_Vac] += 1.0 * R[46];  /* R46: Pt(S) */
    rr[IDX_CO_Pt] += -1.0 * R[47];  /* R47: CO(S) */
    rr[IDX_Pt_Vac] += -1.0 * R[47];  /* R47: Pt(S) */
    rr[IDX_C_Pt] += 1.0 * R[47];  /* R47: C(S) */
    rr[IDX_O_Pt] += 1.0 * R[47];  /* R47: O(S) */
    rr[IDX_NO] += -1.0 * R[48];  /* R48: NO */
    rr[IDX_Pt_Vac] += -1.0 * R[48];  /* R48: Pt(S) */
    rr[IDX_NO_Pt] += 1.0 * R[48];  /* R48: NO(S) */
    rr[IDX_NO_Pt] += -1.0 * R[49];  /* R49: NO(S) */
    rr[IDX_NO] += 1.0 * R[49];  /* R49: NO */
    rr[IDX_Pt_Vac] += 1.0 * R[49];  /* R49: Pt(S) */
    rr[IDX_N_Pt] += -2.0 * R[50];  /* R50: N(S) */
    rr[IDX_N2] += 1.0 * R[50];  /* R50: N2 */
    rr[IDX_Pt_Vac] += 2.0 * R[50];  /* R50: Pt(S) */
    rr[IDX_NO_Pt] += -1.0 * R[51];  /* R51: NO(S) */
    rr[IDX_Pt_Vac] += -1.0 * R[51];  /* R51: Pt(S) */
    rr[IDX_N_Pt] += 1.0 * R[51];  /* R51: N(S) */
    rr[IDX_O_Pt] += 1.0 * R[51];  /* R51: O(S) */
    rr[IDX_N_Pt] += -1.0 * R[52];  /* R52: N(S) */
    rr[IDX_O_Pt] += -1.0 * R[52];  /* R52: O(S) */
    rr[IDX_NO_Pt] += 1.0 * R[52];  /* R52: NO(S) */
    rr[IDX_Pt_Vac] += 1.0 * R[52];  /* R52: Pt(S) */
    rr[IDX_O2] += -1.0 * R[53];  /* R53: O2 */
    rr[IDX_Rh_Vac] += -2.0 * R[53];  /* R53: Rh(S1) */
    rr[IDX_O_Rh] += 2.0 * R[53];  /* R53: O(S1) */
    rr[IDX_CO] += -1.0 * R[54];  /* R54: CO */
    rr[IDX_Rh_Vac] += -1.0 * R[54];  /* R54: Rh(S1) */
    rr[IDX_CO_Rh] += 1.0 * R[54];  /* R54: CO(S1) */
    rr[IDX_NO] += -1.0 * R[55];  /* R55: NO */
    rr[IDX_Rh_Vac] += -1.0 * R[55];  /* R55: Rh(S1) */
    rr[IDX_NO_Rh] += 1.0 * R[55];  /* R55: NO(S1) */
    rr[IDX_O_Rh] += -2.0 * R[56];  /* R56: O(S1) */
    rr[IDX_O2] += 1.0 * R[56];  /* R56: O2 */
    rr[IDX_Rh_Vac] += 2.0 * R[56];  /* R56: Rh(S1) */
    rr[IDX_CO_Rh] += -1.0 * R[57];  /* R57: CO(S1) */
    rr[IDX_CO] += 1.0 * R[57];  /* R57: CO */
    rr[IDX_Rh_Vac] += 1.0 * R[57];  /* R57: Rh(S1) */
    rr[IDX_NO_Rh] += -1.0 * R[58];  /* R58: NO(S1) */
    rr[IDX_NO] += 1.0 * R[58];  /* R58: NO */
    rr[IDX_Rh_Vac] += 1.0 * R[58];  /* R58: Rh(S1) */
    rr[IDX_N_Rh] += -2.0 * R[59];  /* R59: N(S1) */
    rr[IDX_N2] += 1.0 * R[59];  /* R59: N2 */
    rr[IDX_Rh_Vac] += 2.0 * R[59];  /* R59: Rh(S1) */
    rr[IDX_CO_Rh] += -1.0 * R[60];  /* R60: CO(S1) */
    rr[IDX_O_Rh] += -1.0 * R[60];  /* R60: O(S1) */
    rr[IDX_CO2] += 1.0 * R[60];  /* R60: CO2 */
    rr[IDX_Rh_Vac] += 2.0 * R[60];  /* R60: Rh(S1) */
    rr[IDX_NO_Rh] += -1.0 * R[61];  /* R61: NO(S1) */
    rr[IDX_Rh_Vac] += -1.0 * R[61];  /* R61: Rh(S1) */
    rr[IDX_N_Rh] += 1.0 * R[61];  /* R61: N(S1) */
    rr[IDX_O_Rh] += 1.0 * R[61];  /* R61: O(S1) */
    }
}
