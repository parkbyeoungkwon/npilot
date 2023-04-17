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
void live_H(double *in_vec, double *out_1526956633722536676);
void live_err_fun(double *nom_x, double *delta_x, double *out_3230416160365452320);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_3644280506359970566);
void live_H_mod_fun(double *state, double *out_5597765590530589087);
void live_f_fun(double *state, double dt, double *out_1109019082588435211);
void live_F_fun(double *state, double dt, double *out_4209121349994159134);
void live_h_4(double *state, double *unused, double *out_127005917296737428);
void live_H_4(double *state, double *unused, double *out_1751574332541175876);
void live_h_9(double *state, double *unused, double *out_295549267304732882);
void live_H_9(double *state, double *unused, double *out_5908742068895953359);
void live_h_10(double *state, double *unused, double *out_8994351018960577770);
void live_H_10(double *state, double *unused, double *out_881387217823812750);
void live_h_12(double *state, double *unused, double *out_8105328428021366339);
void live_H_12(double *state, double *unused, double *out_1130475307493582209);
void live_h_31(double *state, double *unused, double *out_5160887999776252386);
void live_H_31(double *state, double *unused, double *out_1615087724831431500);
void live_h_32(double *state, double *unused, double *out_4726594278692124627);
void live_H_32(double *state, double *unused, double *out_1746252509772171969);
void live_h_13(double *state, double *unused, double *out_1833653487685933304);
void live_H_13(double *state, double *unused, double *out_5941475568526639255);
void live_h_14(double *state, double *unused, double *out_295549267304732882);
void live_H_14(double *state, double *unused, double *out_5908742068895953359);
void live_h_33(double *state, double *unused, double *out_7478846720422514676);
void live_H_33(double *state, double *unused, double *out_4765644729470289104);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}