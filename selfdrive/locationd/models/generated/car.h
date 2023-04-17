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
void car_err_fun(double *nom_x, double *delta_x, double *out_3212818796643484509);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_5510234247823456276);
void car_H_mod_fun(double *state, double *out_3235254973669982056);
void car_f_fun(double *state, double dt, double *out_146651213496092769);
void car_F_fun(double *state, double dt, double *out_7094746895558734709);
void car_h_25(double *state, double *unused, double *out_212152291061883490);
void car_H_25(double *state, double *unused, double *out_2670607375307301223);
void car_h_24(double *state, double *unused, double *out_6372502206104763617);
void car_H_24(double *state, double *unused, double *out_4083635623345710060);
void car_h_30(double *state, double *unused, double *out_6379266331033017793);
void car_H_30(double *state, double *unused, double *out_1857088954820306975);
void car_h_26(double *state, double *unused, double *out_7065804475015688178);
void car_H_26(double *state, double *unused, double *out_1070895943566755001);
void car_h_27(double *state, double *unused, double *out_1298878939418347817);
void car_H_27(double *state, double *unused, double *out_4031852266620731886);
void car_h_29(double *state, double *unused, double *out_8651204453027463807);
void car_H_29(double *state, double *unused, double *out_1346857610505914791);
void car_h_28(double *state, double *unused, double *out_3846011123697299720);
void car_H_28(double *state, double *unused, double *out_616772661059411460);
void car_h_31(double *state, double *unused, double *out_3849601819315905840);
void car_H_31(double *state, double *unused, double *out_1697104045800106477);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}