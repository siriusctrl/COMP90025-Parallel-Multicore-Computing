#ifndef __OPENMP_H__
#define __OPENMP_H__

double cal_pi();
double cal_pi_non_for();
double cal_pi_atom(int n_thread);
double cal_pi_reduction(int n_thread);

#endif // __OPENMP_H__
