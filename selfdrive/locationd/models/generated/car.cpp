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
void err_fun(double *nom_x, double *delta_x, double *out_5238830884312407088) {
   out_5238830884312407088[0] = delta_x[0] + nom_x[0];
   out_5238830884312407088[1] = delta_x[1] + nom_x[1];
   out_5238830884312407088[2] = delta_x[2] + nom_x[2];
   out_5238830884312407088[3] = delta_x[3] + nom_x[3];
   out_5238830884312407088[4] = delta_x[4] + nom_x[4];
   out_5238830884312407088[5] = delta_x[5] + nom_x[5];
   out_5238830884312407088[6] = delta_x[6] + nom_x[6];
   out_5238830884312407088[7] = delta_x[7] + nom_x[7];
   out_5238830884312407088[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_8143899348751214684) {
   out_8143899348751214684[0] = -nom_x[0] + true_x[0];
   out_8143899348751214684[1] = -nom_x[1] + true_x[1];
   out_8143899348751214684[2] = -nom_x[2] + true_x[2];
   out_8143899348751214684[3] = -nom_x[3] + true_x[3];
   out_8143899348751214684[4] = -nom_x[4] + true_x[4];
   out_8143899348751214684[5] = -nom_x[5] + true_x[5];
   out_8143899348751214684[6] = -nom_x[6] + true_x[6];
   out_8143899348751214684[7] = -nom_x[7] + true_x[7];
   out_8143899348751214684[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_3308815358614772559) {
   out_3308815358614772559[0] = 1.0;
   out_3308815358614772559[1] = 0;
   out_3308815358614772559[2] = 0;
   out_3308815358614772559[3] = 0;
   out_3308815358614772559[4] = 0;
   out_3308815358614772559[5] = 0;
   out_3308815358614772559[6] = 0;
   out_3308815358614772559[7] = 0;
   out_3308815358614772559[8] = 0;
   out_3308815358614772559[9] = 0;
   out_3308815358614772559[10] = 1.0;
   out_3308815358614772559[11] = 0;
   out_3308815358614772559[12] = 0;
   out_3308815358614772559[13] = 0;
   out_3308815358614772559[14] = 0;
   out_3308815358614772559[15] = 0;
   out_3308815358614772559[16] = 0;
   out_3308815358614772559[17] = 0;
   out_3308815358614772559[18] = 0;
   out_3308815358614772559[19] = 0;
   out_3308815358614772559[20] = 1.0;
   out_3308815358614772559[21] = 0;
   out_3308815358614772559[22] = 0;
   out_3308815358614772559[23] = 0;
   out_3308815358614772559[24] = 0;
   out_3308815358614772559[25] = 0;
   out_3308815358614772559[26] = 0;
   out_3308815358614772559[27] = 0;
   out_3308815358614772559[28] = 0;
   out_3308815358614772559[29] = 0;
   out_3308815358614772559[30] = 1.0;
   out_3308815358614772559[31] = 0;
   out_3308815358614772559[32] = 0;
   out_3308815358614772559[33] = 0;
   out_3308815358614772559[34] = 0;
   out_3308815358614772559[35] = 0;
   out_3308815358614772559[36] = 0;
   out_3308815358614772559[37] = 0;
   out_3308815358614772559[38] = 0;
   out_3308815358614772559[39] = 0;
   out_3308815358614772559[40] = 1.0;
   out_3308815358614772559[41] = 0;
   out_3308815358614772559[42] = 0;
   out_3308815358614772559[43] = 0;
   out_3308815358614772559[44] = 0;
   out_3308815358614772559[45] = 0;
   out_3308815358614772559[46] = 0;
   out_3308815358614772559[47] = 0;
   out_3308815358614772559[48] = 0;
   out_3308815358614772559[49] = 0;
   out_3308815358614772559[50] = 1.0;
   out_3308815358614772559[51] = 0;
   out_3308815358614772559[52] = 0;
   out_3308815358614772559[53] = 0;
   out_3308815358614772559[54] = 0;
   out_3308815358614772559[55] = 0;
   out_3308815358614772559[56] = 0;
   out_3308815358614772559[57] = 0;
   out_3308815358614772559[58] = 0;
   out_3308815358614772559[59] = 0;
   out_3308815358614772559[60] = 1.0;
   out_3308815358614772559[61] = 0;
   out_3308815358614772559[62] = 0;
   out_3308815358614772559[63] = 0;
   out_3308815358614772559[64] = 0;
   out_3308815358614772559[65] = 0;
   out_3308815358614772559[66] = 0;
   out_3308815358614772559[67] = 0;
   out_3308815358614772559[68] = 0;
   out_3308815358614772559[69] = 0;
   out_3308815358614772559[70] = 1.0;
   out_3308815358614772559[71] = 0;
   out_3308815358614772559[72] = 0;
   out_3308815358614772559[73] = 0;
   out_3308815358614772559[74] = 0;
   out_3308815358614772559[75] = 0;
   out_3308815358614772559[76] = 0;
   out_3308815358614772559[77] = 0;
   out_3308815358614772559[78] = 0;
   out_3308815358614772559[79] = 0;
   out_3308815358614772559[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_219863977268319485) {
   out_219863977268319485[0] = state[0];
   out_219863977268319485[1] = state[1];
   out_219863977268319485[2] = state[2];
   out_219863977268319485[3] = state[3];
   out_219863977268319485[4] = state[4];
   out_219863977268319485[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_219863977268319485[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_219863977268319485[7] = state[7];
   out_219863977268319485[8] = state[8];
}
void F_fun(double *state, double dt, double *out_8110082327665247600) {
   out_8110082327665247600[0] = 1;
   out_8110082327665247600[1] = 0;
   out_8110082327665247600[2] = 0;
   out_8110082327665247600[3] = 0;
   out_8110082327665247600[4] = 0;
   out_8110082327665247600[5] = 0;
   out_8110082327665247600[6] = 0;
   out_8110082327665247600[7] = 0;
   out_8110082327665247600[8] = 0;
   out_8110082327665247600[9] = 0;
   out_8110082327665247600[10] = 1;
   out_8110082327665247600[11] = 0;
   out_8110082327665247600[12] = 0;
   out_8110082327665247600[13] = 0;
   out_8110082327665247600[14] = 0;
   out_8110082327665247600[15] = 0;
   out_8110082327665247600[16] = 0;
   out_8110082327665247600[17] = 0;
   out_8110082327665247600[18] = 0;
   out_8110082327665247600[19] = 0;
   out_8110082327665247600[20] = 1;
   out_8110082327665247600[21] = 0;
   out_8110082327665247600[22] = 0;
   out_8110082327665247600[23] = 0;
   out_8110082327665247600[24] = 0;
   out_8110082327665247600[25] = 0;
   out_8110082327665247600[26] = 0;
   out_8110082327665247600[27] = 0;
   out_8110082327665247600[28] = 0;
   out_8110082327665247600[29] = 0;
   out_8110082327665247600[30] = 1;
   out_8110082327665247600[31] = 0;
   out_8110082327665247600[32] = 0;
   out_8110082327665247600[33] = 0;
   out_8110082327665247600[34] = 0;
   out_8110082327665247600[35] = 0;
   out_8110082327665247600[36] = 0;
   out_8110082327665247600[37] = 0;
   out_8110082327665247600[38] = 0;
   out_8110082327665247600[39] = 0;
   out_8110082327665247600[40] = 1;
   out_8110082327665247600[41] = 0;
   out_8110082327665247600[42] = 0;
   out_8110082327665247600[43] = 0;
   out_8110082327665247600[44] = 0;
   out_8110082327665247600[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_8110082327665247600[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_8110082327665247600[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8110082327665247600[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_8110082327665247600[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_8110082327665247600[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_8110082327665247600[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_8110082327665247600[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_8110082327665247600[53] = -9.8000000000000007*dt;
   out_8110082327665247600[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_8110082327665247600[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_8110082327665247600[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8110082327665247600[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8110082327665247600[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_8110082327665247600[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_8110082327665247600[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_8110082327665247600[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_8110082327665247600[62] = 0;
   out_8110082327665247600[63] = 0;
   out_8110082327665247600[64] = 0;
   out_8110082327665247600[65] = 0;
   out_8110082327665247600[66] = 0;
   out_8110082327665247600[67] = 0;
   out_8110082327665247600[68] = 0;
   out_8110082327665247600[69] = 0;
   out_8110082327665247600[70] = 1;
   out_8110082327665247600[71] = 0;
   out_8110082327665247600[72] = 0;
   out_8110082327665247600[73] = 0;
   out_8110082327665247600[74] = 0;
   out_8110082327665247600[75] = 0;
   out_8110082327665247600[76] = 0;
   out_8110082327665247600[77] = 0;
   out_8110082327665247600[78] = 0;
   out_8110082327665247600[79] = 0;
   out_8110082327665247600[80] = 1;
}
void h_25(double *state, double *unused, double *out_6658966305743246042) {
   out_6658966305743246042[0] = state[6];
}
void H_25(double *state, double *unused, double *out_3520617132260708968) {
   out_3520617132260708968[0] = 0;
   out_3520617132260708968[1] = 0;
   out_3520617132260708968[2] = 0;
   out_3520617132260708968[3] = 0;
   out_3520617132260708968[4] = 0;
   out_3520617132260708968[5] = 0;
   out_3520617132260708968[6] = 1;
   out_3520617132260708968[7] = 0;
   out_3520617132260708968[8] = 0;
}
void h_24(double *state, double *unused, double *out_5997382673890137841) {
   out_5997382673890137841[0] = state[4];
   out_5997382673890137841[1] = state[5];
}
void H_24(double *state, double *unused, double *out_8267126459227218002) {
   out_8267126459227218002[0] = 0;
   out_8267126459227218002[1] = 0;
   out_8267126459227218002[2] = 0;
   out_8267126459227218002[3] = 0;
   out_8267126459227218002[4] = 1;
   out_8267126459227218002[5] = 0;
   out_8267126459227218002[6] = 0;
   out_8267126459227218002[7] = 0;
   out_8267126459227218002[8] = 0;
   out_8267126459227218002[9] = 0;
   out_8267126459227218002[10] = 0;
   out_8267126459227218002[11] = 0;
   out_8267126459227218002[12] = 0;
   out_8267126459227218002[13] = 0;
   out_8267126459227218002[14] = 1;
   out_8267126459227218002[15] = 0;
   out_8267126459227218002[16] = 0;
   out_8267126459227218002[17] = 0;
}
void h_30(double *state, double *unused, double *out_5016942541615954638) {
   out_5016942541615954638[0] = state[4];
}
void H_30(double *state, double *unused, double *out_1002284173753460341) {
   out_1002284173753460341[0] = 0;
   out_1002284173753460341[1] = 0;
   out_1002284173753460341[2] = 0;
   out_1002284173753460341[3] = 0;
   out_1002284173753460341[4] = 1;
   out_1002284173753460341[5] = 0;
   out_1002284173753460341[6] = 0;
   out_1002284173753460341[7] = 0;
   out_1002284173753460341[8] = 0;
}
void h_26(double *state, double *unused, double *out_2309032993367267470) {
   out_2309032993367267470[0] = state[7];
}
void H_26(double *state, double *unused, double *out_7262120451134765192) {
   out_7262120451134765192[0] = 0;
   out_7262120451134765192[1] = 0;
   out_7262120451134765192[2] = 0;
   out_7262120451134765192[3] = 0;
   out_7262120451134765192[4] = 0;
   out_7262120451134765192[5] = 0;
   out_7262120451134765192[6] = 0;
   out_7262120451134765192[7] = 1;
   out_7262120451134765192[8] = 0;
}
void h_27(double *state, double *unused, double *out_6221676583905795709) {
   out_6221676583905795709[0] = state[3];
}
void H_27(double *state, double *unused, double *out_1221309897430482876) {
   out_1221309897430482876[0] = 0;
   out_1221309897430482876[1] = 0;
   out_1221309897430482876[2] = 0;
   out_1221309897430482876[3] = 1;
   out_1221309897430482876[4] = 0;
   out_1221309897430482876[5] = 0;
   out_1221309897430482876[6] = 0;
   out_1221309897430482876[7] = 0;
   out_1221309897430482876[8] = 0;
}
void h_29(double *state, double *unused, double *out_5033518376948244111) {
   out_5033518376948244111[0] = state[1];
}
void H_29(double *state, double *unused, double *out_492052829439068157) {
   out_492052829439068157[0] = 0;
   out_492052829439068157[1] = 1;
   out_492052829439068157[2] = 0;
   out_492052829439068157[3] = 0;
   out_492052829439068157[4] = 0;
   out_492052829439068157[5] = 0;
   out_492052829439068157[6] = 0;
   out_492052829439068157[7] = 0;
   out_492052829439068157[8] = 0;
}
void h_28(double *state, double *unused, double *out_6999160273968430960) {
   out_6999160273968430960[0] = state[0];
}
void H_28(double *state, double *unused, double *out_5574451846508598731) {
   out_5574451846508598731[0] = 1;
   out_5574451846508598731[1] = 0;
   out_5574451846508598731[2] = 0;
   out_5574451846508598731[3] = 0;
   out_5574451846508598731[4] = 0;
   out_5574451846508598731[5] = 0;
   out_5574451846508598731[6] = 0;
   out_5574451846508598731[7] = 0;
   out_5574451846508598731[8] = 0;
}
void h_31(double *state, double *unused, double *out_1980911170975295842) {
   out_1980911170975295842[0] = state[8];
}
void H_31(double *state, double *unused, double *out_7888328553368116668) {
   out_7888328553368116668[0] = 0;
   out_7888328553368116668[1] = 0;
   out_7888328553368116668[2] = 0;
   out_7888328553368116668[3] = 0;
   out_7888328553368116668[4] = 0;
   out_7888328553368116668[5] = 0;
   out_7888328553368116668[6] = 0;
   out_7888328553368116668[7] = 0;
   out_7888328553368116668[8] = 1;
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
void car_err_fun(double *nom_x, double *delta_x, double *out_5238830884312407088) {
  err_fun(nom_x, delta_x, out_5238830884312407088);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_8143899348751214684) {
  inv_err_fun(nom_x, true_x, out_8143899348751214684);
}
void car_H_mod_fun(double *state, double *out_3308815358614772559) {
  H_mod_fun(state, out_3308815358614772559);
}
void car_f_fun(double *state, double dt, double *out_219863977268319485) {
  f_fun(state,  dt, out_219863977268319485);
}
void car_F_fun(double *state, double dt, double *out_8110082327665247600) {
  F_fun(state,  dt, out_8110082327665247600);
}
void car_h_25(double *state, double *unused, double *out_6658966305743246042) {
  h_25(state, unused, out_6658966305743246042);
}
void car_H_25(double *state, double *unused, double *out_3520617132260708968) {
  H_25(state, unused, out_3520617132260708968);
}
void car_h_24(double *state, double *unused, double *out_5997382673890137841) {
  h_24(state, unused, out_5997382673890137841);
}
void car_H_24(double *state, double *unused, double *out_8267126459227218002) {
  H_24(state, unused, out_8267126459227218002);
}
void car_h_30(double *state, double *unused, double *out_5016942541615954638) {
  h_30(state, unused, out_5016942541615954638);
}
void car_H_30(double *state, double *unused, double *out_1002284173753460341) {
  H_30(state, unused, out_1002284173753460341);
}
void car_h_26(double *state, double *unused, double *out_2309032993367267470) {
  h_26(state, unused, out_2309032993367267470);
}
void car_H_26(double *state, double *unused, double *out_7262120451134765192) {
  H_26(state, unused, out_7262120451134765192);
}
void car_h_27(double *state, double *unused, double *out_6221676583905795709) {
  h_27(state, unused, out_6221676583905795709);
}
void car_H_27(double *state, double *unused, double *out_1221309897430482876) {
  H_27(state, unused, out_1221309897430482876);
}
void car_h_29(double *state, double *unused, double *out_5033518376948244111) {
  h_29(state, unused, out_5033518376948244111);
}
void car_H_29(double *state, double *unused, double *out_492052829439068157) {
  H_29(state, unused, out_492052829439068157);
}
void car_h_28(double *state, double *unused, double *out_6999160273968430960) {
  h_28(state, unused, out_6999160273968430960);
}
void car_H_28(double *state, double *unused, double *out_5574451846508598731) {
  H_28(state, unused, out_5574451846508598731);
}
void car_h_31(double *state, double *unused, double *out_1980911170975295842) {
  h_31(state, unused, out_1980911170975295842);
}
void car_H_31(double *state, double *unused, double *out_7888328553368116668) {
  H_31(state, unused, out_7888328553368116668);
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
