#include <iostream>
#include <iomanip>
#include <Eigen/Dense>

#include "MPC_controller.h"

using namespace Eigen;


int main(void) {

	//Structure that contains parameters for the MPC controller
	controller_param params;

	//The physical dimensions of the inverted pendulum and the slide
	double l_rod = 1; //Length of rod in meters
	double R_rod = 0.01; //Radius of rod in meters
	double m_rod = 1; //Mass of rod in kilograms
	double m_slide = 2; //Mass of slide in kilograms
	double g = 9.81; //Gravitational acceleration m/s^2

	//Equation of motion of the system linearized around the upper equilibrium of the pendulum
	//dx/dt = A_cont * x + B_cont * u, where x contains the angle and angular velocity of the rod
	//and the horizontal position and velocity of the slide
	double C = R_rod * R_rod * m_rod / 4.0 + l_rod * l_rod * m_rod / 12.0 +
		l_rod * l_rod * m_rod * m_slide / (4.0 * (m_rod + m_slide)); //Auxiliary constant 

	MatrixXd A_cont(4, 4), B_cont(4, 1);
	A_cont << 0, 1, 0, 0,
		g* l_rod* m_rod / (2.0 * C), 0, 0, 0,
		0, 0, 0, 1,
		g* l_rod* l_rod* m_rod* m_rod / (4.0 * (m_rod + m_slide) * C), 0, 0, 0;
	B_cont << 0, l_rod* m_rod / (C * 2 * (m_rod + m_slide)), 0,
		m_rod* (3.0 * R_rod * R_rod + 4.0 * l_rod * l_rod) / (12.0 * C * (m_rod + m_slide));

	//Discretizing dynamics using Euler discretization
	double T_s = 0.1; //Sampling time
	MatrixXd A_d(4, 4), B_d(4, 1);
	A_d = MatrixXd::Identity(4, 4) + T_s * A_cont;
	B_d = T_s * B_cont;
	
	params.A = A_d;
	params.B = B_d;

	//Defining appropriate state cost matrix
	MatrixXd Q = MatrixXd::Zero(4, 4);
	Q.diagonal() << 1000, 1, 1000, 1;
	params.Q = Q;

	//Defining input penalty as a matrix
	params.R = MatrixXd::Identity(1, 1);

	//Defining terminal state cost
	params.P = 1000.0 * MatrixXd::Identity(4, 4);

	//Time horizon length
	params.N = 15; //1.5 seconds

	//State constraints only on the position of the slide
	//Set constraints to infinity at unconstrained states
	const double inf = std::numeric_limits<double>::infinity();
	params.x_min.resize(4);
	params.x_max.resize(4);
	params.x_min << -inf, -inf, -1, -inf;
	params.x_max << inf, inf, 1, inf;

	//Input constraints
	params.u_min.resize(1);
	params.u_max.resize(1);
	params.u_min = -20.0 * VectorXd::Ones(1);
	params.u_max = 20.0 * VectorXd::Ones(1);
	
	//Creating MPC controller object
	MPC_controller ctrl(params);

	//Defining a current state
	VectorXd x0(4);
	x0 << 0.1, 0, 0, 0;

	//Calculating optimal control input
	VectorXd u = ctrl.solve(x0);

	//Printing results
	std::cout << "Initial state:" << std::endl;
	std::cout << x0 << std::endl << std::endl;
	std::cout << "Optimal control input:" << std::endl;
	std::cout << u << std::endl;

}