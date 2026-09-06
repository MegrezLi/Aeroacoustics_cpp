#pragma once
#ifdef _WIN32
#define AERO_EXPORT __declspec(dllexport)
#else
#define AERO_EXPORT
#endif
#ifdef __cplusplus
extern "C" {
#endif
/* op: 1 LBLVS,2 TBLTE,3 TIPNOIS,4 InflowNoise,5 BLUNT,6 Guidati,7 TNO.
 * x[26]: alpha,chord,U,theta,phi,span,R,stall,d99P,dstarS,dstarP,TI,h,psi,
 * thick1%,thick10%,sound_speed,nu,rho,Lturb,alprat,CfS,CfP,d99S,UeS,UeP.
 * Alpha is radians ONLY for op 4, degrees otherwise.
 * flags[3]: BL method, trip, rounded tip. Outputs concatenated by mechanism:
 * op2 [pressure,suction,separation], op7 [pressure,suction]; each n frequencies.
 * Returns 0 on success; caller allocates 3*n doubles. Errors never cross the ABI.
 */
AERO_EXPORT int aeroacoustics_kernel(int op,int n,const double* freq,const double* x,const int* flags,double* output);
AERO_EXPORT const char* aeroacoustics_last_error(void);
AERO_EXPORT const char* aeroacoustics_backend(void);
AERO_EXPORT unsigned long long aeroacoustics_call_count(void);
#ifdef __cplusplus
}
#endif
