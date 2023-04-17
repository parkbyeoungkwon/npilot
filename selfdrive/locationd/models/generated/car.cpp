#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with sympy 1.9                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_3212818796643484509) {
   out_3212818796643484509[0] = delta_x[0] + nom_x[0];
   out_3212818796643484509[1] = delta_x[1] + nom_x[1];
   out_3212818796643484509[2] = delta_x[2] + nom_x[2];
   out_3212818796643484509[3] = delta_x[3] + nom_x[3];
   out_3212818796643484509[4] = delta_x[4] + nom_x[4];
   out_3212818796643484509[5] = delta_x[5] + nom_x[5];
   out_3212818796643484509[6] = delta_x[6] + nom_x[6];
   out_3212818796643484509[7] = delta_x[7] + nom_x[7];
   out_3212818796643484509[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_5510234247823456276) {
   out_5510234247823456276[0] = -nom_x[0] + true_x[0];
   out_5510234247823456276[1] = -nom_x[1] + true_x[1];
   out_5510234247823456276[2] = -nom_x[2] + true_x[2];
   out_5510234247823456276[3] = -nom_x[3] + true_x[3];
   out_5510234247823456276[4] = -nom_x[4] + true_x[4];
   out_5510234247823456276[5] = -nom_x[5] + true_x[5];
   out_5510234247823456276[6] = -nom_x[6] + true_x[6];
   out_5510234247823456276[7] = -nom_x[7] + true_x[7];
   out_5510234247823456276[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_3235254973669982056) {
   out_3235254973669982056[0] = 1.0;
   out_3235254973669982056[1] = 0;
   out_3235254973669982056[2] = 0;
   out_3235254973669982056[3] = 0;
   out_3235254973669982056[4] = 0;
   out_3235254973669982056[5] = 0;
   out_3235254973669982056[6] = 0;
   out_3235254973669982056[7] = 0;
   out_3235254973669982056[8] = 0;
   out_3235254973669982056[9] = 0;
   out_3235254973669982056[10] = 1.0;
   out_3235254973669982056[11] = 0;
   out_3235254973669982056[12] = 0;
   out_3235254973669982056[13] = 0;
   out_3235254973669982056[14] = 0;
   out_3235254973669982056[15] = 0;
   out_3235254973669982056[16] = 0;
   out_3235254973669982056[17] = 0;
   out_3235254973669982056[18] = 0;
   out_3235254973669982056[19] = 0;
   out_3235254973669982056[20] = 1.0;
   out_3235254973669982056[21] = 0;
   out_3235254973669982056[22] = 0;
   out_3235254973669982056[23] = 0;
   out_3235254973669982056[24] = 0;
   out_3235254973669982056[25] = 0;
   out_3235254973669982056[26] = 0;
   out_3235254973669982056[27] = 0;
   out_3235254973669982056[28] = 0;
   out_3235254973669982056[29] = 0;
   out_3235254973669982056[30] = 1.0;
   out_3235254973669982056[31] = 0;
   out_3235254973669982056[32] = 0;
   out_3235254973669982056[33] = 0;
   out_3235254973669982056[34] = 0;
   out_3235254973669982056[35] = 0;
   out_3235254973669982056[36] = 0;
   out_3235254973669982056[37] = 0;
   out_3235254973669982056[38] = 0;
   out_3235254973669982056[39] = 0;
   out_3235254973669982056[40] = 1.0;
   out_3235254973669982056[41] = 0;
   out_3235254973669982056[42] = 0;
   out_3235254973669982056[43] = 0;
   out_3235254973669982056[44] = 0;
   out_3235254973669982056[45] = 0;
   out_3235254973669982056[46] = 0;
   out_3235254973669982056[47] = 0;
   out_3235254973669982056[48] = 0;
   out_3235254973669982056[49] = 0;
   out_3235254973669982056[50] = 1.0;
   out_3235254973669982056[51] = 0;
   out_3235254973669982056[52] = 0;
   out_3235254973669982056[53] = 0;
   out_3235254973669982056[54] = 0;
   out_3235254973669982056[55] = 0;
   out_3235254973669982056[56] = 0;
   out_3235254973669982056[57] = 0;
   out_3235254973669982056[58] = 0;
   out_3235254973669982056[59] = 0;
   out_3235254973669982056[60] = 1.0;
   out_3235254973669982056[61] = 0;
   out_3235254973669982056[62] = 0;
   out_3235254973669982056[63] = 0;
   out_3235254973669982056[64] = 0;
   out_3235254973669982056[65] = 0;
   out_3235254973669982056[66] = 0;
   out_3235254973669982056[67] = 0;
   out_3235254973669982056[68] = 0;
   out_3235254973669982056[69] = 0;
   out_3235254973669982056[70] = 1.0;
   out_3235254973669982056[71] = 0;
   out_3235254973669982056[72] = 0;
   out_3235254973669982056[73] = 0;
   out_3235254973669982056[74] = 0;
   out_3235254973669982056[75] = 0;
   out_3235254973669982056[76] = 0;
   out_3235254973669982056[77] = 0;
   out_3235254973669982056[78] = 0;
   out_3235254973669982056[79] = 0;
   out_3235254973669982056[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_146651213496092769) {
   out_146651213496092769[0] = state[0];
   out_146651213496092769[1] = state[1];
   out_146651213496092769[2] = state[2];
   out_146651213496092769[3] = state[3];
   out_146651213496092769[4] = state[4];
   out_146651213496092769[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_146651213496092769[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_146651213496092769[7] = state[7];
   out_146651213496092769[8] = state[8];
}
void F_fun(double *state, double dt, double *out_7094746895558734709) {
   out_7094746895558734709[0] = 1;
   out_7094746895558734709[1] = 0;
   out_7094746895558734709[2] = 0;
   out_7094746895558734709[3] = 0;
   out_7094746895558734709[4] = 0;
   out_7094746895558734709[5] = 0;
   out_7094746895558734709[6] = 0;
   out_7094746895558734709[7] = 0;
   out_7094746895558734709[8] = 0;
   out_7094746895558734709[9] = 0;
   out_7094746895558734709[10] = 1;
   out_7094746895558734709[11] = 0;
   out_7094746895558734709[12] = 0;
   out_7094746895558734709[13] = 0;
   out_7094746895558734709[14] = 0;
   out_7094746895558734709[15] = 0;
   out_7094746895558734709[16] = 0;
   out_7094746895558734709[17] = 0;
   out_7094746895558734709[18] = 0;
   out_7094746895558734709[19] = 0;
   out_7094746895558734709[20] = 1;
   out_7094746895558734709[21] = 0;
   out_7094746895558734709[22] = 0;
   out_7094746895558734709[23] = 0;
   out_7094746895558734709[24] = 0;
   out_7094746895558734709[25] = 0;
   out_7094746895558734709[26] = 0;
   out_7094746895558734709[27] = 0;
   out_7094746895558734709[28] = 0;
   out_7094746895558734709[29] = 0;
   out_7094746895558734709[30] = 1;
   out_7094746895558734709[31] = 0;
   out_7094746895558734709[32] = 0;
   out_7094746895558734709[33] = 0;
   out_7094746895558734709[34] = 0;
   out_7094746895558734709[35] = 0;
   out_7094746895558734709[36] = 0;
   out_7094746895558734709[37] = 0;
   out_7094746895558734709[38] = 0;
   out_7094746895558734709[39] = 0;
   out_7094746895558734709[40] = 1;
   out_7094746895558734709[41] = 0;
   out_7094746895558734709[42] = 0;
   out_7094746895558734709[43] = 0;
   out_7094746895558734709[44] = 0;
   out_7094746895558734709[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_7094746895558734709[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_7094746895558734709[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7094746895558734709[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_7094746895558734709[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_7094746895558734709[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_7094746895558734709[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_7094746895558734709[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_7094746895558734709[53] = -9.8000000000000007*dt;
   out_7094746895558734709[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_7094746895558734709[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_7094746895558734709[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7094746895558734709[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7094746895558734709[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_7094746895558734709[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_7094746895558734709[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_7094746895558734709[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_7094746895558734709[62] = 0;
   out_7094746895558734709[63] = 0;
   out_7094746895558734709[64] = 0;
   out_7094746895558734709[65] = 0;
   out_7094746895558734709[66] = 0;
   out_7094746895558734709[67] = 0;
   out_7094746895558734709[68] = 0;
   out_7094746895558734709[69] = 0;
   out_7094746895558734709[70] = 1;
   out_7094746895558734709[71] = 0;
   out_7094746895558734709[72] = 0;
   out_7094746895558734709[73] = 0;
   out_7094746895558734709[74] = 0;
   out_7094746895558734709[75] = 0;
   out_7094746895558734709[76] = 0;
   out_7094746895558734709[77] = 0;
   out_7094746895558734709[78] = 0;
   out_7094746895558734709[79] = 0;
   out_7094746895558734709[80] = 1;
}
void h_25(double *state, double *unused, double *out_212152291061883490) {
   out_212152291061883490[0] = state[6];
}
void H_25(double *state, double *unused, double *out_2670607375307301223) {
   out_2670607375307301223[0] = 0;
   out_2670607375307301223[1] = 0;
   out_2670607375307301223[2] = 0;
   out_2670607375307301223[3] = 0;
   out_2670607375307301223[4] = 0;
   out_2670607375307301223[5] = 0;
   out_2670607375307301223[6] = 1;
   out_2670607375307301223[7] = 0;
   out_2670607375307301223[8] = 0;
}
void h_24(double *state, double *unused, double *out_6372502206104763617) {
   out_6372502206104763617[0] = state[4];
   out_6372502206104763617[1] = state[5];
}
void H_24(double *state, double *unused, double *out_4083635623345710060) {
   out_4083635623345710060[0] = 0;
   out_4083635623345710060[1] = 0;
   out_4083635623345710060[2] = 0;
   out_4083635623345710060[3] = 0;
   out_4083635623345710060[4] = 1;
   out_4083635623345710060[5] = 0;
   out_4083635623345710060[6] = 0;
   out_4083635623345710060[7] = 0;
   out_4083635623345710060[8] = 0;
   out_4083635623345710060[9] = 0;
   out_4083635623345710060[10] = 0;
   out_4083635623345710060[11] = 0;
   out_4083635623345710060[12] = 0;
   out_4083635623345710060[13] = 0;
   out_4083635623345710060[14] = 1;
   out_4083635623345710060[15] = 0;
   out_4083635623345710060[16] = 0;
   out_4083635623345710060[17] = 0;
}
void h_30(double *state, double *unused, double *out_6379266331033017793) {
   out_6379266331033017793[0] = state[4];
}
void H_30(double *state, double *unused, double *out_1857088954820306975) {
   out_1857088954820306975[0] = 0;
   out_1857088954820306975[1] = 0;
   out_1857088954820306975[2] = 0;
   out_1857088954820306975[3] = 0;
   out_1857088954820306975[4] = 1;
   out_1857088954820306975[5] = 0;
   out_1857088954820306975[6] = 0;
   out_1857088954820306975[7] = 0;
   out_1857088954820306975[8] = 0;
}
void h_26(double *state, double *unused, double *out_7065804475015688178) {
   out_7065804475015688178[0] = state[7];
}
void H_26(double *state, double *unused, double *out_1070895943566755001) {
   out_1070895943566755001[0] = 0;
   out_1070895943566755001[1] = 0;
   out_1070895943566755001[2] = 0;
   out_1070895943566755001[3] = 0;
   out_1070895943566755001[4] = 0;
   out_1070895943566755001[5] = 0;
   out_1070895943566755001[6] = 0;
   out_1070895943566755001[7] = 1;
   out_1070895943566755001[8] = 0;
}
void h_27(double *state, double *unused, double *out_1298878939418347817) {
   out_1298878939418347817[0] = state[3];
}
void H_27(double *state, double *unused, double *out_4031852266620731886) {
   out_4031852266620731886[0] = 0;
   out_4031852266620731886[1] = 0;
   out_4031852266620731886[2] = 0;
   out_4031852266620731886[3] = 1;
   out_4031852266620731886[4] = 0;
   out_4031852266620731886[5] = 0;
   out_4031852266620731886[6] = 0;
   out_4031852266620731886[7] = 0;
   out_4031852266620731886[8] = 0;
}
void h_29(double *state, double *unused, double *out_8651204453027463807) {
   out_8651204453027463807[0] = state[1];
}
void H_29(double *state, double *unused, double *out_1346857610505914791) {
   out_1346857610505914791[0] = 0;
   out_1346857610505914791[1] = 1;
   out_1346857610505914791[2] = 0;
   out_1346857610505914791[3] = 0;
   out_1346857610505914791[4] = 0;
   out_1346857610505914791[5] = 0;
   out_1346857610505914791[6] = 0;
   out_1346857610505914791[7] = 0;
   out_1346857610505914791[8] = 0;
}
void h_28(double *state, double *unused, double *out_3846011123697299720) {
   out_3846011123697299720[0] = state[0];
}
void H_28(double *state, double *unused, double *out_616772661059411460) {
   out_616772661059411460[0] = 1;
   out_616772661059411460[1] = 0;
   out_616772661059411460[2] = 0;
   out_616772661059411460[3] = 0;
   out_616772661059411460[4] = 0;
   out_616772661059411460[5] = 0;
   out_616772661059411460[6] = 0;
   out_616772661059411460[7] = 0;
   out_616772661059411460[8] = 0;
}
void h_31(double *state, double *unused, double *out_3849601819315905840) {
   out_3849601819315905840[0] = state[8];
}
void H_31(double *state, double *unused, double *out_1697104045800106477) {
   out_1697104045800106477[0] = 0;
   out_1697104045800106477[1] = 0;
   out_1697104045800106477[2] = 0;
   out_1697104045800106477[3] = 0;
   out_1697104045800106477[4] = 0;
   out_1697104045800106477[5] = 0;
   out_1697104045800106477[6] = 0;
   out_1697104045800106477[7] = 0;
   out_1697104045800106477[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_3212818796643484509) {
  err_fun(nom_x, delta_x, out_3212818796643484509);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_5510234247823456276) {
  inv_err_fun(nom_x, true_x, out_5510234247823456276);
}
void car_H_mod_fun(double *state, double *out_3235254973669982056) {
  H_mod_fun(state, out_3235254973669982056);
}
void car_f_fun(double *state, double dt, double *out_146651213496092769) {
  f_fun(state,  dt, out_146651213496092769);
}
void car_F_fun(double *state, double dt, double *out_7094746895558734709) {
  F_fun(state,  dt, out_7094746895558734709);
}
void car_h_25(double *state, double *unused, double *out_212152291061883490) {
  h_25(state, unused, out_212152291061883490);
}
void car_H_25(double *state, double *unused, double *out_2670607375307301223) {
  H_25(state, unused, out_2670607375307301223);
}
void car_h_24(double *state, double *unused, double *out_6372502206104763617) {
  h_24(state, unused, out_6372502206104763617);
}
void car_H_24(double *state, double *unused, double *out_4083635623345710060) {
  H_24(state, unused, out_4083635623345710060);
}
void car_h_30(double *state, double *unused, double *out_6379266331033017793) {
  h_30(state, unused, out_6379266331033017793);
}
void car_H_30(double *state, double *unused, double *out_1857088954820306975) {
  H_30(state, unused, out_1857088954820306975);
}
void car_h_26(double *state, double *unused, double *out_7065804475015688178) {
  h_26(state, unused, out_7065804475015688178);
}
void car_H_26(double *state, double *unused, double *out_1070895943566755001) {
  H_26(state, unused, out_1070895943566755001);
}
void car_h_27(double *state, double *unused, double *out_1298878939418347817) {
  h_27(state, unused, out_1298878939418347817);
}
void car_H_27(double *state, double *unused, double *out_4031852266620731886) {
  H_27(state, unused, out_4031852266620731886);
}
void car_h_29(double *state, double *unused, double *out_8651204453027463807) {
  h_29(state, unused, out_8651204453027463807);
}
void car_H_29(double *state, double *unused, double *out_1346857610505914791) {
  H_29(state, unused, out_1346857610505914791);
}
void car_h_28(double *state, double *unused, double *out_3846011123697299720) {
  h_28(state, unused, out_3846011123697299720);
}
void car_H_28(double *state, double *unused, double *out_616772661059411460) {
  H_28(state, unused, out_616772661059411460);
}
void car_h_31(double *state, double *unused, double *out_3849601819315905840) {
  h_31(state, unused, out_3849601819315905840);
}
void car_H_31(double *state, double *unused, double *out_1697104045800106477) {
  H_31(state, unused, out_1697104045800106477);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
