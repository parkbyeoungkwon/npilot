#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_3856071707786988958);
void live_err_fun(double *nom_x, double *delta_x, double *out_6083673053396515662);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_5989358798141820053);
void live_H_mod_fun(double *state, double *out_6752933429230893614);
void live_f_fun(double *state, double dt, double *out_4349171889326912352);
void live_F_fun(double *state, double dt, double *out_6295655534664862528);
void live_h_4(double *state, double *unused, double *out_7589561656465924740);
void live_H_4(double *state, double *unused, double *out_4332055988748560851);
void live_h_9(double *state, double *unused, double *out_963519030010363412);
void live_H_9(double *state, double *unused, double *out_4090866342118970206);
void live_h_10(double *state, double *unused, double *out_1873821913890153121);
void live_H_10(double *state, double *unused, double *out_3372944201445050667);
void live_h_12(double *state, double *unused, double *out_2527492314541715057);
void live_H_12(double *state, double *unused, double *out_687400419283400944);
void live_h_31(double *state, double *unused, double *out_1525015336238717412);
void live_H_31(double *state, double *unused, double *out_3432963451608414653);
void live_h_32(double *state, double *unused, double *out_2026364883032920571);
void live_H_32(double *state, double *unused, double *out_6137987964510164590);
void live_h_13(double *state, double *unused, double *out_4717903246551368742);
void live_H_13(double *state, double *unused, double *out_1475130968955094840);
void live_h_14(double *state, double *unused, double *out_963519030010363412);
void live_H_14(double *state, double *unused, double *out_4090866342118970206);
void live_h_33(double *state, double *unused, double *out_6774399766974439182);
void live_H_33(double *state, double *unused, double *out_6583520456247272257);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}