/****************************************************************************
 *
 *   Copyright (c) 2016 PX4 Development Team. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in
 *    the documentation and/or other materials provided with the
 *    distribution.
 * 3. Neither the name PX4 nor the names of its contributors may be
 *    used to endorse or promote products derived from this software
 *    without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
 * OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 * AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 ****************************************************************************/

#include "../IEKF.hpp"

void IEKF::correctFlow(const optical_flow_s *msg)
{
	// requires terrain estimate
	if (!getTerrainValid()) {
		return;
	}

	//ROS_INFO("correct flow");
	// return if no new data
	float dt = 0;

	if (!_sensorFlow.ready(msg->timestamp, dt)) {
		return;
	}

	// compute agl
	float agl = _x(X::asl) - _x(X::terrain_asl);

	// return if too close to ground
	if (agl < 0.2f) {
		return;
	}

	// init global reference
	if (agl > 0.2f && !_origin.xyInitialized()) {
		float lat_deg = 47.397742f;
		float lon_deg = 8.545594f;
		ROS_INFO("flow origin init lat: %12.6f deg lon: %12.6f deg",
			 double(lat_deg), double(lon_deg));
		_origin.xyInitialize(lat_deg, lon_deg, msg->timestamp);
	}

	// state info
	Dcmf C_nb = Quaternion<float>(
			    _x(X::q_nb_0), _x(X::q_nb_1),
			    _x(X::q_nb_2), _x(X::q_nb_3));

	Vector3f omega_nb_b(
		_u(U::omega_nb_bX), _u(U::omega_nb_bY), _u(U::omega_nb_bZ));

	float vel_N = _x(X::vel_N);
	float vel_E = _x(X::vel_E);
	float vel_D = _x(X::vel_D);

	// expected measurement
	Vector<float, Y_flow::n> yh;
	yh(0) = - (C_nb(0, 1) * vel_N + C_nb(1, 1) * vel_E
		   + C_nb(2, 1) * vel_D) * C_nb(2, 2) / agl + omega_nb_b(1);
	yh(1) = - (C_nb(0, 0) * vel_N + C_nb(1, 0) * vel_E
		   + C_nb(2, 0) * vel_D) * C_nb(2, 2) / agl + omega_nb_b(0);

	// measurement
	Vector<float, Y_flow::n> y;
	float flow_dt = msg->integration_timespan / 1.0e6f;
	y(0) = msg->pixel_flow_x_integral / flow_dt;
	y(1) = -msg->pixel_flow_y_integral / flow_dt;

	//ROS_INFO("flow X %10.4f Y %10.4f", double(y(0)), double(y(1)));
	//ROS_INFO("vel N %10.4f E %10.4f", double(vel_N), double(vel_E));

	// residual
	Vector<float, Y_flow::n> r = y - yh;

	//ROS_INFO("flow: (%10g, %10g)", double(y(0)), double(y(1)));
	//ROS_INFO("float dt: %10.4f", double(dt));

	// define R
	Matrix<float, Y_flow::n, Y_flow::n> R;
	R(Y_flow::flowX, Y_flow::flowX) = flow_sigma_rw * flow_sigma_rw / dt;
	R(Y_flow::flowY, Y_flow::flowY) = flow_sigma_rw * flow_sigma_rw / dt;

	// define H
	// Note: this measurement is not invariant due to
	// rotation matrix
	Matrix<float, Y_flow::n, Xe::n> H;
	float x0 = 1 / agl;
	float x1 = 2 * C_nb(1, 2);
	float x2 = C_nb(0, 1) * vel_N;
	float x3 = C_nb(1, 1) * vel_E;
	float x4 = 2 * C_nb(2, 2) * vel_D;
	float x5 = C_nb(2, 1) * vel_D;
	float x6 = 2 * C_nb(2, 2) * vel_E;
	float x7 = 2 * x0;
	float x8 = C_nb(2, 2) * vel_D;
	float x9 = C_nb(2, 2) * vel_N;
	float x10 = 2 * C_nb(2, 2) * x0;
	float x11 = C_nb(2, 2) * x0;
	float x12 = C_nb(2, 2) / agl / agl;
	float x13 = x12 * (x2 + x3 + x5);
	float x14 = C_nb(0, 0) * vel_N;
	float x15 = C_nb(1, 0) * vel_E;
	float x16 = C_nb(2, 0) * vel_D;
	float x17 = x12 * (x14 + x15 + x16);
	H(Y_flow::flowX, Xe::rot_N) = x0 * (-C_nb(1, 1) * x4
					    + C_nb(2, 1) * x6 - x1 * x2 - x1 * x3 - x1 * x5);
	H(Y_flow::flowX, Xe::rot_E) = x7 * (C_nb(0, 1) * x8 + C_nb(0, 2) * x2
					    + C_nb(0, 2) * x3 + C_nb(0, 2) * x5 - C_nb(2, 1) * x9);
	H(Y_flow::flowX, Xe::rot_D) = -x10 * (C_nb(0, 1) * vel_E - C_nb(1, 1) * vel_N);
	H(Y_flow::flowX, Xe::vel_N) = -C_nb(0, 1) * x11;
	H(Y_flow::flowX, Xe::vel_E) = -C_nb(1, 1) * x11;
	H(Y_flow::flowX, Xe::vel_D) = -C_nb(2, 1) * x11;
	H(Y_flow::flowX, Xe::gyro_bias_N) = -C_nb(0, 0);
	H(Y_flow::flowX, Xe::gyro_bias_E) = -C_nb(1, 0);
	H(Y_flow::flowX, Xe::gyro_bias_D) = -C_nb(2, 0);
	H(Y_flow::flowX, Xe::asl) = x13;
	H(Y_flow::flowX, Xe::terrain_asl) = -x13;
	H(Y_flow::flowY, Xe::rot_N) = x0 * (-C_nb(1, 0) * x4
					    + C_nb(2, 0) * x6 - x1 * x14 - x1 * x15 - x1 * x16);
	H(Y_flow::flowY, Xe::rot_E) = x7 * (C_nb(0, 0) * x8
					    + C_nb(0, 2) * x14 + C_nb(0, 2) * x15
					    + C_nb(0, 2) * x16 - C_nb(2, 0) * x9);
	H(Y_flow::flowY, Xe::rot_D) = -x10 * (C_nb(0, 0) * vel_E - C_nb(1, 0) * vel_N);
	H(Y_flow::flowY, Xe::vel_N) = -C_nb(0, 0) * x11;
	H(Y_flow::flowY, Xe::vel_E) = -C_nb(1, 0) * x11;
	H(Y_flow::flowY, Xe::vel_D) = -C_nb(2, 0) * x11;
	H(Y_flow::flowY, Xe::gyro_bias_N) = C_nb(0, 1);
	H(Y_flow::flowY, Xe::gyro_bias_E) = C_nb(1, 1);
	H(Y_flow::flowY, Xe::gyro_bias_D) = C_nb(2, 1);
	H(Y_flow::flowY, Xe::asl) = x17;
	H(Y_flow::flowY, Xe::terrain_asl) = -x17;

	// kalman correction
	Vector<float, Xe::n> dxe;
	SquareMatrix<float, Xe::n> dP;
	_sensorFlow.kalmanCorrectCond(_P, H, R, r, dxe, dP);

	if (_sensorFlow.shouldCorrect()) {
		// only correct velocity
		for (int i = 0; i < Xe::n; i++) {
			if (i != Xe::vel_N && i != Xe::vel_D) {
				dxe(i) = 0;
			}
		}

		setX(applyErrorCorrection(dxe));
		setP(_P + dP);
	}
}
