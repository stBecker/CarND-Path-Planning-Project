#ifndef MPC_H
#define MPC_H

#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/Dense"
#include "Eigen-3.3/Eigen/Eigen"

using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;


struct State {
  double p, v, a;
  State(double p_, double v_, double a_) {
    p = p_;
    v = v_;
    a = a_;
  };
};


class JMT {

public:
  VectorXd coeff_;

  JMT(State initial_state, State final_state, double t) {
    MatrixXd A = MatrixXd(3, 3);
    VectorXd b = VectorXd(3);
    VectorXd x = VectorXd(3);
    coeff_ = VectorXd(6);

    const double t2 = t * t;
    const double t3 = t * t2;
    const double t4 = t * t3;
    const double t5 = t * t4;

    A << 
      t3, t4, t5,
      3 * t2, 4 * t3, 5 * t4,
      6 * t, 12 * t2, 20 * t3;

    b << 
      final_state.p - (initial_state.p + initial_state.v * t + 0.5 * initial_state.a * t2),
      final_state.v - (initial_state.v + initial_state.a * t),
      final_state.a - initial_state.a;

    x = A.inverse() * b;

    coeff_ << 
      initial_state.p,
      initial_state.v,
      initial_state.a / 2.0,
      x[0],
      x[1],
      x[2];
  
  };

  double Evaluate(double t) {
    const double t2 = t * t;
    const double t3 = t * t2;
    const double t4 = t * t3;
    const double t5 = t * t4;

    Eigen::VectorXd T = VectorXd(6);
    T << 1.0, t, t2, t3, t4, t5;

    return T.transpose() * coeff_;
  };
};

#endif /* MPC_H */
