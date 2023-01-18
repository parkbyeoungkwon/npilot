#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_3276638675524358659);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_6624635078069616892);
void car_H_mod_fun(double *state, double *out_5697227785510417705);
void car_f_fun(double *state, double dt, double *out_2314309569023769267);
void car_F_fun(double *state, double dt, double *out_4231482060143634001);
void car_h_25(double *state, double *unused, double *out_3800024309291946346);
void car_H_25(double *state, double *unused, double *out_484264312352185458);
void car_h_24(double *state, double *unused, double *out_8354558351955910234);
void car_H_24(double *state, double *unused, double *out_1688385286653314108);
void car_h_30(double *state, double *unused, double *out_199560818977883780);
void car_H_30(double *state, double *unused, double *out_4043432017775422740);
void car_h_26(double *state, double *unused, double *out_8605217638622110433);
void car_H_26(double *state, double *unused, double *out_3257239006521870766);
void car_h_27(double *state, double *unused, double *out_4782108402627244458);
void car_H_27(double *state, double *unused, double *out_6218195329575847651);
void car_h_29(double *state, double *unused, double *out_1988726823723106478);
void car_H_29(double *state, double *unused, double *out_3533200673461030556);
void car_h_28(double *state, double *unused, double *out_8021905800001350323);
void car_H_28(double *state, double *unused, double *out_1569570401895704305);
void car_h_31(double *state, double *unused, double *out_1231066527595281980);
void car_H_31(double *state, double *unused, double *out_3883447108755222242);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}