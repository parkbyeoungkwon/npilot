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
void car_err_fun(double *nom_x, double *delta_x, double *out_4527979494488094157);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_8670707127363216169);
void car_H_mod_fun(double *state, double *out_6553232040040776566);
void car_f_fun(double *state, double dt, double *out_8513490376515808541);
void car_F_fun(double *state, double dt, double *out_529957493507700796);
void car_h_25(double *state, double *unused, double *out_2278518569996568210);
void car_H_25(double *state, double *unused, double *out_2705346990763697114);
void car_h_24(double *state, double *unused, double *out_6144112832497924479);
void car_H_24(double *state, double *unused, double *out_532697391758197548);
void car_h_30(double *state, double *unused, double *out_5226244560352053171);
void car_H_30(double *state, double *unused, double *out_1822349339363911084);
void car_h_26(double *state, double *unused, double *out_4238530341350657680);
void car_H_26(double *state, double *unused, double *out_1036156328110359110);
void car_h_27(double *state, double *unused, double *out_7280191837657413855);
void car_H_27(double *state, double *unused, double *out_3997112651164335995);
void car_h_29(double *state, double *unused, double *out_5829560445028530688);
void car_H_29(double *state, double *unused, double *out_1312117995049518900);
void car_h_28(double *state, double *unused, double *out_4243405786213081170);
void car_H_28(double *state, double *unused, double *out_651512276515807351);
void car_h_31(double *state, double *unused, double *out_9049353796346919146);
void car_H_31(double *state, double *unused, double *out_1662364430343710586);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}