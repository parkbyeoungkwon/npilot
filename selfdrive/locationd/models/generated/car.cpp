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
void err_fun(double *nom_x, double *delta_x, double *out_3276638675524358659) {
   out_3276638675524358659[0] = delta_x[0] + nom_x[0];
   out_3276638675524358659[1] = delta_x[1] + nom_x[1];
   out_3276638675524358659[2] = delta_x[2] + nom_x[2];
   out_3276638675524358659[3] = delta_x[3] + nom_x[3];
   out_3276638675524358659[4] = delta_x[4] + nom_x[4];
   out_3276638675524358659[5] = delta_x[5] + nom_x[5];
   out_3276638675524358659[6] = delta_x[6] + nom_x[6];
   out_3276638675524358659[7] = delta_x[7] + nom_x[7];
   out_3276638675524358659[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_6624635078069616892) {
   out_6624635078069616892[0] = -nom_x[0] + true_x[0];
   out_6624635078069616892[1] = -nom_x[1] + true_x[1];
   out_6624635078069616892[2] = -nom_x[2] + true_x[2];
   out_6624635078069616892[3] = -nom_x[3] + true_x[3];
   out_6624635078069616892[4] = -nom_x[4] + true_x[4];
   out_6624635078069616892[5] = -nom_x[5] + true_x[5];
   out_6624635078069616892[6] = -nom_x[6] + true_x[6];
   out_6624635078069616892[7] = -nom_x[7] + true_x[7];
   out_6624635078069616892[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_5697227785510417705) {
   out_5697227785510417705[0] = 1.0;
   out_5697227785510417705[1] = 0;
   out_5697227785510417705[2] = 0;
   out_5697227785510417705[3] = 0;
   out_5697227785510417705[4] = 0;
   out_5697227785510417705[5] = 0;
   out_5697227785510417705[6] = 0;
   out_5697227785510417705[7] = 0;
   out_5697227785510417705[8] = 0;
   out_5697227785510417705[9] = 0;
   out_5697227785510417705[10] = 1.0;
   out_5697227785510417705[11] = 0;
   out_5697227785510417705[12] = 0;
   out_5697227785510417705[13] = 0;
   out_5697227785510417705[14] = 0;
   out_5697227785510417705[15] = 0;
   out_5697227785510417705[16] = 0;
   out_5697227785510417705[17] = 0;
   out_5697227785510417705[18] = 0;
   out_5697227785510417705[19] = 0;
   out_5697227785510417705[20] = 1.0;
   out_5697227785510417705[21] = 0;
   out_5697227785510417705[22] = 0;
   out_5697227785510417705[23] = 0;
   out_5697227785510417705[24] = 0;
   out_5697227785510417705[25] = 0;
   out_5697227785510417705[26] = 0;
   out_5697227785510417705[27] = 0;
   out_5697227785510417705[28] = 0;
   out_5697227785510417705[29] = 0;
   out_5697227785510417705[30] = 1.0;
   out_5697227785510417705[31] = 0;
   out_5697227785510417705[32] = 0;
   out_5697227785510417705[33] = 0;
   out_5697227785510417705[34] = 0;
   out_5697227785510417705[35] = 0;
   out_5697227785510417705[36] = 0;
   out_5697227785510417705[37] = 0;
   out_5697227785510417705[38] = 0;
   out_5697227785510417705[39] = 0;
   out_5697227785510417705[40] = 1.0;
   out_5697227785510417705[41] = 0;
   out_5697227785510417705[42] = 0;
   out_5697227785510417705[43] = 0;
   out_5697227785510417705[44] = 0;
   out_5697227785510417705[45] = 0;
   out_5697227785510417705[46] = 0;
   out_5697227785510417705[47] = 0;
   out_5697227785510417705[48] = 0;
   out_5697227785510417705[49] = 0;
   out_5697227785510417705[50] = 1.0;
   out_5697227785510417705[51] = 0;
   out_5697227785510417705[52] = 0;
   out_5697227785510417705[53] = 0;
   out_5697227785510417705[54] = 0;
   out_5697227785510417705[55] = 0;
   out_5697227785510417705[56] = 0;
   out_5697227785510417705[57] = 0;
   out_5697227785510417705[58] = 0;
   out_5697227785510417705[59] = 0;
   out_5697227785510417705[60] = 1.0;
   out_5697227785510417705[61] = 0;
   out_5697227785510417705[62] = 0;
   out_5697227785510417705[63] = 0;
   out_5697227785510417705[64] = 0;
   out_5697227785510417705[65] = 0;
   out_5697227785510417705[66] = 0;
   out_5697227785510417705[67] = 0;
   out_5697227785510417705[68] = 0;
   out_5697227785510417705[69] = 0;
   out_5697227785510417705[70] = 1.0;
   out_5697227785510417705[71] = 0;
   out_5697227785510417705[72] = 0;
   out_5697227785510417705[73] = 0;
   out_5697227785510417705[74] = 0;
   out_5697227785510417705[75] = 0;
   out_5697227785510417705[76] = 0;
   out_5697227785510417705[77] = 0;
   out_5697227785510417705[78] = 0;
   out_5697227785510417705[79] = 0;
   out_5697227785510417705[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_2314309569023769267) {
   out_2314309569023769267[0] = state[0];
   out_2314309569023769267[1] = state[1];
   out_2314309569023769267[2] = state[2];
   out_2314309569023769267[3] = state[3];
   out_2314309569023769267[4] = state[4];
   out_2314309569023769267[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_2314309569023769267[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_2314309569023769267[7] = state[7];
   out_2314309569023769267[8] = state[8];
}
void F_fun(double *state, double dt, double *out_4231482060143634001) {
   out_4231482060143634001[0] = 1;
   out_4231482060143634001[1] = 0;
   out_4231482060143634001[2] = 0;
   out_4231482060143634001[3] = 0;
   out_4231482060143634001[4] = 0;
   out_4231482060143634001[5] = 0;
   out_4231482060143634001[6] = 0;
   out_4231482060143634001[7] = 0;
   out_4231482060143634001[8] = 0;
   out_4231482060143634001[9] = 0;
   out_4231482060143634001[10] = 1;
   out_4231482060143634001[11] = 0;
   out_4231482060143634001[12] = 0;
   out_4231482060143634001[13] = 0;
   out_4231482060143634001[14] = 0;
   out_4231482060143634001[15] = 0;
   out_4231482060143634001[16] = 0;
   out_4231482060143634001[17] = 0;
   out_4231482060143634001[18] = 0;
   out_4231482060143634001[19] = 0;
   out_4231482060143634001[20] = 1;
   out_4231482060143634001[21] = 0;
   out_4231482060143634001[22] = 0;
   out_4231482060143634001[23] = 0;
   out_4231482060143634001[24] = 0;
   out_4231482060143634001[25] = 0;
   out_4231482060143634001[26] = 0;
   out_4231482060143634001[27] = 0;
   out_4231482060143634001[28] = 0;
   out_4231482060143634001[29] = 0;
   out_4231482060143634001[30] = 1;
   out_4231482060143634001[31] = 0;
   out_4231482060143634001[32] = 0;
   out_4231482060143634001[33] = 0;
   out_4231482060143634001[34] = 0;
   out_4231482060143634001[35] = 0;
   out_4231482060143634001[36] = 0;
   out_4231482060143634001[37] = 0;
   out_4231482060143634001[38] = 0;
   out_4231482060143634001[39] = 0;
   out_4231482060143634001[40] = 1;
   out_4231482060143634001[41] = 0;
   out_4231482060143634001[42] = 0;
   out_4231482060143634001[43] = 0;
   out_4231482060143634001[44] = 0;
   out_4231482060143634001[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_4231482060143634001[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_4231482060143634001[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_4231482060143634001[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_4231482060143634001[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_4231482060143634001[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_4231482060143634001[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_4231482060143634001[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_4231482060143634001[53] = -9.8000000000000007*dt;
   out_4231482060143634001[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_4231482060143634001[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_4231482060143634001[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4231482060143634001[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4231482060143634001[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_4231482060143634001[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_4231482060143634001[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_4231482060143634001[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_4231482060143634001[62] = 0;
   out_4231482060143634001[63] = 0;
   out_4231482060143634001[64] = 0;
   out_4231482060143634001[65] = 0;
   out_4231482060143634001[66] = 0;
   out_4231482060143634001[67] = 0;
   out_4231482060143634001[68] = 0;
   out_4231482060143634001[69] = 0;
   out_4231482060143634001[70] = 1;
   out_4231482060143634001[71] = 0;
   out_4231482060143634001[72] = 0;
   out_4231482060143634001[73] = 0;
   out_4231482060143634001[74] = 0;
   out_4231482060143634001[75] = 0;
   out_4231482060143634001[76] = 0;
   out_4231482060143634001[77] = 0;
   out_4231482060143634001[78] = 0;
   out_4231482060143634001[79] = 0;
   out_4231482060143634001[80] = 1;
}
void h_25(double *state, double *unused, double *out_3800024309291946346) {
   out_3800024309291946346[0] = state[6];
}
void H_25(double *state, double *unused, double *out_484264312352185458) {
   out_484264312352185458[0] = 0;
   out_484264312352185458[1] = 0;
   out_484264312352185458[2] = 0;
   out_484264312352185458[3] = 0;
   out_484264312352185458[4] = 0;
   out_484264312352185458[5] = 0;
   out_484264312352185458[6] = 1;
   out_484264312352185458[7] = 0;
   out_484264312352185458[8] = 0;
}
void h_24(double *state, double *unused, double *out_8354558351955910234) {
   out_8354558351955910234[0] = state[4];
   out_8354558351955910234[1] = state[5];
}
void H_24(double *state, double *unused, double *out_1688385286653314108) {
   out_1688385286653314108[0] = 0;
   out_1688385286653314108[1] = 0;
   out_1688385286653314108[2] = 0;
   out_1688385286653314108[3] = 0;
   out_1688385286653314108[4] = 1;
   out_1688385286653314108[5] = 0;
   out_1688385286653314108[6] = 0;
   out_1688385286653314108[7] = 0;
   out_1688385286653314108[8] = 0;
   out_1688385286653314108[9] = 0;
   out_1688385286653314108[10] = 0;
   out_1688385286653314108[11] = 0;
   out_1688385286653314108[12] = 0;
   out_1688385286653314108[13] = 0;
   out_1688385286653314108[14] = 1;
   out_1688385286653314108[15] = 0;
   out_1688385286653314108[16] = 0;
   out_1688385286653314108[17] = 0;
}
void h_30(double *state, double *unused, double *out_199560818977883780) {
   out_199560818977883780[0] = state[4];
}
void H_30(double *state, double *unused, double *out_4043432017775422740) {
   out_4043432017775422740[0] = 0;
   out_4043432017775422740[1] = 0;
   out_4043432017775422740[2] = 0;
   out_4043432017775422740[3] = 0;
   out_4043432017775422740[4] = 1;
   out_4043432017775422740[5] = 0;
   out_4043432017775422740[6] = 0;
   out_4043432017775422740[7] = 0;
   out_4043432017775422740[8] = 0;
}
void h_26(double *state, double *unused, double *out_8605217638622110433) {
   out_8605217638622110433[0] = state[7];
}
void H_26(double *state, double *unused, double *out_3257239006521870766) {
   out_3257239006521870766[0] = 0;
   out_3257239006521870766[1] = 0;
   out_3257239006521870766[2] = 0;
   out_3257239006521870766[3] = 0;
   out_3257239006521870766[4] = 0;
   out_3257239006521870766[5] = 0;
   out_3257239006521870766[6] = 0;
   out_3257239006521870766[7] = 1;
   out_3257239006521870766[8] = 0;
}
void h_27(double *state, double *unused, double *out_4782108402627244458) {
   out_4782108402627244458[0] = state[3];
}
void H_27(double *state, double *unused, double *out_6218195329575847651) {
   out_6218195329575847651[0] = 0;
   out_6218195329575847651[1] = 0;
   out_6218195329575847651[2] = 0;
   out_6218195329575847651[3] = 1;
   out_6218195329575847651[4] = 0;
   out_6218195329575847651[5] = 0;
   out_6218195329575847651[6] = 0;
   out_6218195329575847651[7] = 0;
   out_6218195329575847651[8] = 0;
}
void h_29(double *state, double *unused, double *out_1988726823723106478) {
   out_1988726823723106478[0] = state[1];
}
void H_29(double *state, double *unused, double *out_3533200673461030556) {
   out_3533200673461030556[0] = 0;
   out_3533200673461030556[1] = 1;
   out_3533200673461030556[2] = 0;
   out_3533200673461030556[3] = 0;
   out_3533200673461030556[4] = 0;
   out_3533200673461030556[5] = 0;
   out_3533200673461030556[6] = 0;
   out_3533200673461030556[7] = 0;
   out_3533200673461030556[8] = 0;
}
void h_28(double *state, double *unused, double *out_8021905800001350323) {
   out_8021905800001350323[0] = state[0];
}
void H_28(double *state, double *unused, double *out_1569570401895704305) {
   out_1569570401895704305[0] = 1;
   out_1569570401895704305[1] = 0;
   out_1569570401895704305[2] = 0;
   out_1569570401895704305[3] = 0;
   out_1569570401895704305[4] = 0;
   out_1569570401895704305[5] = 0;
   out_1569570401895704305[6] = 0;
   out_1569570401895704305[7] = 0;
   out_1569570401895704305[8] = 0;
}
void h_31(double *state, double *unused, double *out_1231066527595281980) {
   out_1231066527595281980[0] = state[8];
}
void H_31(double *state, double *unused, double *out_3883447108755222242) {
   out_3883447108755222242[0] = 0;
   out_3883447108755222242[1] = 0;
   out_3883447108755222242[2] = 0;
   out_3883447108755222242[3] = 0;
   out_3883447108755222242[4] = 0;
   out_3883447108755222242[5] = 0;
   out_3883447108755222242[6] = 0;
   out_3883447108755222242[7] = 0;
   out_3883447108755222242[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_3276638675524358659) {
  err_fun(nom_x, delta_x, out_3276638675524358659);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_6624635078069616892) {
  inv_err_fun(nom_x, true_x, out_6624635078069616892);
}
void car_H_mod_fun(double *state, double *out_5697227785510417705) {
  H_mod_fun(state, out_5697227785510417705);
}
void car_f_fun(double *state, double dt, double *out_2314309569023769267) {
  f_fun(state,  dt, out_2314309569023769267);
}
void car_F_fun(double *state, double dt, double *out_4231482060143634001) {
  F_fun(state,  dt, out_4231482060143634001);
}
void car_h_25(double *state, double *unused, double *out_3800024309291946346) {
  h_25(state, unused, out_3800024309291946346);
}
void car_H_25(double *state, double *unused, double *out_484264312352185458) {
  H_25(state, unused, out_484264312352185458);
}
void car_h_24(double *state, double *unused, double *out_8354558351955910234) {
  h_24(state, unused, out_8354558351955910234);
}
void car_H_24(double *state, double *unused, double *out_1688385286653314108) {
  H_24(state, unused, out_1688385286653314108);
}
void car_h_30(double *state, double *unused, double *out_199560818977883780) {
  h_30(state, unused, out_199560818977883780);
}
void car_H_30(double *state, double *unused, double *out_4043432017775422740) {
  H_30(state, unused, out_4043432017775422740);
}
void car_h_26(double *state, double *unused, double *out_8605217638622110433) {
  h_26(state, unused, out_8605217638622110433);
}
void car_H_26(double *state, double *unused, double *out_3257239006521870766) {
  H_26(state, unused, out_3257239006521870766);
}
void car_h_27(double *state, double *unused, double *out_4782108402627244458) {
  h_27(state, unused, out_4782108402627244458);
}
void car_H_27(double *state, double *unused, double *out_6218195329575847651) {
  H_27(state, unused, out_6218195329575847651);
}
void car_h_29(double *state, double *unused, double *out_1988726823723106478) {
  h_29(state, unused, out_1988726823723106478);
}
void car_H_29(double *state, double *unused, double *out_3533200673461030556) {
  H_29(state, unused, out_3533200673461030556);
}
void car_h_28(double *state, double *unused, double *out_8021905800001350323) {
  h_28(state, unused, out_8021905800001350323);
}
void car_H_28(double *state, double *unused, double *out_1569570401895704305) {
  H_28(state, unused, out_1569570401895704305);
}
void car_h_31(double *state, double *unused, double *out_1231066527595281980) {
  h_31(state, unused, out_1231066527595281980);
}
void car_H_31(double *state, double *unused, double *out_3883447108755222242) {
  H_31(state, unused, out_3883447108755222242);
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
