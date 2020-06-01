#pragma once

#include <Eigen/Dense>

//Using Eigen based QP solver. Can be downloaded from https://www.cs.cmu.edu/~bstephe1/eiquadprog.hpp
#include "eiquadprog.hpp" 

using namespace Eigen;


//Function to compute Kronecker product
MatrixXd kron(const MatrixXd& A, const MatrixXd& B) {

	MatrixXd C = MatrixXd::Zero(A.rows() * B.rows(), A.cols() * B.cols());

	for (int i = 0; i < A.rows(); i++) {
		for (int j = 0; j < A.cols(); j++) {
			C.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) = A(i, j) * B;
		}
	}

	return C;
}

//Structure containing the parameters needed for the MPC controller
struct controller_param {

	MatrixXd A; //State-space representation of dynamics: x_i+1 = A * x_i + B * u_i
	MatrixXd B;
	MatrixXd Q; //State penalty
	MatrixXd R; //Input penalty
	MatrixXd P; //Terminal state penalty
	int N; //Horizon length
	VectorXd x_min; //State constraints
	VectorXd x_max;
	VectorXd u_min; //Input constraints
	VectorXd u_max;

};

//Class for MPC controller. It solves the constrained finite time optimal control problem by reformulating it as
//a QP without substituting the dynamics.
class MPC_controller {

public:

	//Constructor that initializes the optimization problem
	MPC_controller(const controller_param& params) :
		A(params.A), N(params.N), d_x(params.A.rows()), d_u(params.B.cols())
	{
		//Make sure that dimensions match
		assert(params.A.rows() == params.A.cols() && "A must be square!");
		assert(params.A.rows() == params.B.rows() && "A and B must have the same number of rows!");
		assert(params.A.rows() == params.Q.rows() && params.A.cols() == params.Q.cols() &&
			"A and Q must have the same dimensions!");
		assert(params.R.cols() == params.B.cols() && params.R.cols() == params.R.rows() &&
			"R must be square and its dimensions must be consistent with B!");
		assert(params.P.rows() == params.Q.rows() && params.P.cols() == params.Q.cols() &&
			"P and Q must have the same dimensions!");
		assert(params.x_min.rows() == params.A.rows() && params.x_max.rows() == params.A.rows() &&
			"State constraints must be consistent with A!");
		assert(params.u_min.rows() == params.B.cols() && params.u_max.rows() == params.B.cols() &&
			"State constraints must be consistent with A!");

		//The optimization variables are x1,...,xN and u0,...,uN-1
		int size_x = (N - 1) * d_x; //Size of block corresponding to states
		int size_u = (N - 1) * d_u; //Size of block corresponding to inputs

		//The augmented cost matrix contains the state, input and terminal state costs
		H = MatrixXd::Zero(size_x + size_u, size_x + size_u);
		H.topLeftCorner((N - 2) * d_x, (N - 2) * d_x) = kron(MatrixXd::Identity(N - 2, N - 2), params.Q);
		H.bottomRightCorner(size_u, size_u) = kron(MatrixXd::Identity(N - 1, N - 1), params.R);
		H.block((N - 2) * d_x, (N - 2) * d_x, d_x, d_x) = params.P;

		//The matrix for the equality constraints contain the dynamics at each step
		MatrixXd aux_mat = MatrixXd::Zero(N - 1, N - 1); //Just an auxiliary matrix to create G_eq
		aux_mat.diagonal(-1) = VectorXd::Ones(N - 2);

		G_eq = MatrixXd::Zero(size_x, size_x + size_u);
		G_eq.leftCols(size_x) = MatrixXd::Identity(size_x, size_x) - kron(aux_mat, params.A);
		G_eq.rightCols(size_u) = -1.0 * kron(MatrixXd::Identity(N - 1, N - 1), params.B);
		
		//The matrix for the inequality constraints including trivially true constraints as well for now
		MatrixXd G_ineq_full = MatrixXd::Zero(2 * (size_x + size_u), size_x + size_u);
		G_ineq_full.block(0, 0, size_x, size_x) = MatrixXd::Identity(size_x, size_x);
		G_ineq_full.block(size_x, 0, size_x, size_x) = -1.0 * MatrixXd::Identity(size_x, size_x);
		G_ineq_full.block(2 * size_x, size_x, size_u, size_u) = MatrixXd::Identity(size_u, size_u);
		G_ineq_full.block(2 * size_x + size_u, size_x, size_u, size_u) = -1.0 * MatrixXd::Identity(size_u, size_u);
		
		//Maximal and minimal values
		VectorXd b_ineq_full = VectorXd::Zero(2 * (size_x + size_u));
		b_ineq_full.segment(0, size_x) = params.x_max.replicate(N - 1, 1);
		b_ineq_full.segment(size_x, size_x) = -1.0 * params.x_min.replicate(N - 1, 1);
		b_ineq_full.segment(2 * size_x, size_u) = params.u_max.replicate(N - 1, 1);
		b_ineq_full.segment(2 * size_x + size_u, size_u) = -1.0 * params.u_min.replicate(N - 1, 1);

		//Counting not trivially true constraints
		int ntt_consts = 0;
		//State maximum constraints
		for (int i = 0; i < d_x; i++) {
			if (params.x_max(i) != std::numeric_limits<double>::infinity()) { ntt_consts++; }
		}
		//State minimum constraints
		for (int i = 0; i < d_x; i++) {
			if (params.x_min(i) != -1.0 * std::numeric_limits<double>::infinity()) { ntt_consts++; }
		}
		//Input maximum constraints
		for (int i = 0; i < d_u; i++) {
			if (params.u_max(i) != std::numeric_limits<double>::infinity()) { ntt_consts++; }
		}
		//Input minimum constraints
		for (int i = 0; i < d_u; i++) {
			if (params.u_min(i) != -1.0 * std::numeric_limits<double>::infinity()) { ntt_consts++; }
		}

		//Initialize the G_ineq and b_ineq with the right size
		G_ineq = MatrixXd::Zero((N - 1) * ntt_consts, size_x + size_u);
		b_ineq = VectorXd::Zero((N - 1) * ntt_consts);

		int j = 0;

		//Fill G_ineq and b_ineq with the non-trivial constraints
		for (int i = 0; i < 2 * (size_x + size_u); i++) {
			if (b_ineq_full(i) != std::numeric_limits<double>::infinity()) {
				G_ineq.row(j) = G_ineq_full.row(i);
				b_ineq(j) = b_ineq_full(i);
				j++;
			}
		}

		//Creating a zero vector for the linear term in the QP
		f0 = VectorXd::Zero((N - 1) * (d_x + d_u));

		//Precomputing the Cholesky decomposition of H for efficiency
		H_chol.compute(H);

		//Compute the trace of H for the QP solver
		H_tr = H.trace();
	}

	//Solver method for solving the problem with given initial state
	VectorXd solve(const VectorXd& x0) {

		//Initial state is needed for the dynamics, so for the equality constraints as well
		b_eq = VectorXd::Zero((N - 1) * d_x);
		b_eq.segment(0, d_x) = A * x0;
		
		//Allocate memory for the solution
		VectorXd x((N - 1) * (d_x + d_u));

		//Solve the QP
		solve_quadprog2(H_chol, H_tr, f0, G_eq.transpose(), -b_eq, -G_ineq.transpose(), b_ineq, x);

		//Return only the first control input of the sequence
		return x.segment((N - 1) * d_x, d_u);
	}
private:

	MatrixXd A, H, G_eq, G_ineq;
	VectorXd u, b_eq, b_ineq, f0;
	int N, d_x, d_u;
	double H_tr;
	LLT<MatrixXd, Lower> H_chol; //Cholesky decomposition of the cost matrix

};

