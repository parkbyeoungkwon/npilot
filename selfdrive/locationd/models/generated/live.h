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
void live_H(double *in_vec, double *out_7402070304101322236);
void live_err_fun(double *nom_x, double *delta_x, double *out_460665167769025615);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_6625240037991151673);
void live_H_mod_fun(double *state, double *out_6172559208002130840);
void live_f_fun(double *state, double dt, double *out_3016561442420806540);
void live_F_fun(double *state, double dt, double *out_6833536681640660446);
void live_h_4(double *state, double *unused, double *out_1804053058314886230);
void live_H_4(double *state, double *unused, double *out_1970134765174672164);
void live_h_9(double *state, double *unused, double *out_1266080828111930518);
void live_H_9(double *state, double *unused, double *out_9189390373270431982);
void live_h_10(double *state, double *unused, double *out_2426365442037283271);
void live_H_10(double *state, double *unused, double *out_6739971248024717518);
void live_h_12(double *state, double *unused, double *out_178634462741645431);
void live_H_12(double *state, double *unused, double *out_4411123611868060832);
void live_h_31(double *state, double *unused, double *out_8342889678123320940);
void live_H_31(double *state, double *unused, double *out_1665560579543047123);
void live_h_32(double *state, double *unused, double *out_8151265232605344810);
void live_H_32(double *state, double *unused, double *out_3691957312836379938);
void live_h_13(double *state, double *unused, double *out_6316801585668289047);
void live_H_13(double *state, double *unused, double *out_5625168891511364287);
void live_h_14(double *state, double *unused, double *out_1266080828111930518);
void live_H_14(double *state, double *unused, double *out_9189390373270431982);
void live_h_33(double *state, double *unused, double *out_6178890621864734831);
void live_H_33(double *state, double *unused, double *out_2913360957888557647);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}