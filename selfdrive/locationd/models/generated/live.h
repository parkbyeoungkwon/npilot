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
void live_H(double *in_vec, double *out_7776051101100143741);
void live_err_fun(double *nom_x, double *delta_x, double *out_1100551690712499759);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_147841421639408624);
void live_H_mod_fun(double *state, double *out_7407040773614503798);
void live_f_fun(double *state, double dt, double *out_562507417187797478);
void live_F_fun(double *state, double dt, double *out_3286422164419981033);
void live_h_4(double *state, double *unused, double *out_7149850188331743031);
void live_H_4(double *state, double *unused, double *out_3214534449346219015);
void live_h_9(double *state, double *unused, double *out_2669286062041613761);
void live_H_9(double *state, double *unused, double *out_6103396001626298357);
void live_h_10(double *state, double *unused, double *out_8886796444612020855);
void live_H_10(double *state, double *unused, double *out_2138968565216532455);
void live_h_12(double *state, double *unused, double *out_7318216388374232033);
void live_H_12(double *state, double *unused, double *out_3835633474393812682);
void live_h_31(double *state, double *unused, double *out_5660530859326060859);
void live_H_31(double *state, double *unused, double *out_6581196506718826391);
void live_h_32(double *state, double *unused, double *out_317394649047073963);
void live_H_32(double *state, double *unused, double *out_2035860113937450807);
void live_h_13(double *state, double *unused, double *out_818184101308438788);
void live_H_13(double *state, double *unused, double *out_6083956319532643888);
void live_h_14(double *state, double *unused, double *out_2669286062041613761);
void live_H_14(double *state, double *unused, double *out_6103396001626298357);
void live_h_33(double *state, double *unused, double *out_3058552237343974245);
void live_H_33(double *state, double *unused, double *out_8714990562351867621);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}