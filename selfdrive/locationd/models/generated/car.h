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
void car_err_fun(double *nom_x, double *delta_x, double *out_5238830884312407088);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_8143899348751214684);
void car_H_mod_fun(double *state, double *out_3308815358614772559);
void car_f_fun(double *state, double dt, double *out_219863977268319485);
void car_F_fun(double *state, double dt, double *out_8110082327665247600);
void car_h_25(double *state, double *unused, double *out_6658966305743246042);
void car_H_25(double *state, double *unused, double *out_3520617132260708968);
void car_h_24(double *state, double *unused, double *out_5997382673890137841);
void car_H_24(double *state, double *unused, double *out_8267126459227218002);
void car_h_30(double *state, double *unused, double *out_5016942541615954638);
void car_H_30(double *state, double *unused, double *out_1002284173753460341);
void car_h_26(double *state, double *unused, double *out_2309032993367267470);
void car_H_26(double *state, double *unused, double *out_7262120451134765192);
void car_h_27(double *state, double *unused, double *out_6221676583905795709);
void car_H_27(double *state, double *unused, double *out_1221309897430482876);
void car_h_29(double *state, double *unused, double *out_5033518376948244111);
void car_H_29(double *state, double *unused, double *out_492052829439068157);
void car_h_28(double *state, double *unused, double *out_6999160273968430960);
void car_H_28(double *state, double *unused, double *out_5574451846508598731);
void car_h_31(double *state, double *unused, double *out_1980911170975295842);
void car_H_31(double *state, double *unused, double *out_7888328553368116668);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}