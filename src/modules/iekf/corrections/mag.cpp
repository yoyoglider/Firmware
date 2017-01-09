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

void IEKF::correctMag(const sensor_combined_s *msg)
{
	//ROS_INFO("correct mag");
	// return if no new data
	float dt = 0;
	uint64_t timestamp = msg->timestamp + msg->magnetometer_timestamp_relative;

	if (!_sensorMag.ready(timestamp, dt)) {
		return;
	}

	// calculate residual
	Quatf q_nb(_x(X::q_nb_0), _x(X::q_nb_1),
		   _x(X::q_nb_2), _x(X::q_nb_3));
	Vector3f y_b = Vector3f(
			       msg->magnetometer_ga[0],
			       msg->magnetometer_ga[1],
			       msg->magnetometer_ga[2]).unit();
	Vector3f yh = Dcmf(Eulerf(
				   0,
				   -1 * deg2radf * _magInclDeg, /* pitch down convention */
				   deg2radf * _magDeclDeg)).T() * Vector3f(1, 0, 0);
	Vector3f y = q_nb.conjugate(y_b);
	Vector2f r(y(0) - yh(0), y(1) - yh(1));

	// define R
	float covNE = mag_sigma_rw * mag_sigma_rw / dt;
	Matrix<float, Y_mag::n, Y_mag::n> R;
	R(Y_mag::mag_N, Y_mag::mag_N) = covNE;
	R(Y_mag::mag_E, Y_mag::mag_E) = covNE;

	// define H
	Matrix<float, Y_mag::n, Xe::n> H;
	Matrix3f tmp = yh.hat() * 2;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			H(Y_mag::mag_N + i, Xe::rot_N + j) = tmp(i, j);
		}
	}

	// kalman correction
	Vector<float, Xe::n> dxe;
	SquareMatrix<float, Xe::n> dP;
	_sensorMag.kalmanCorrectCond(_P, H, R, r, dxe, dP);
	Vector<float, X::n> x = applyErrorCorrection(dxe);

	if (_sensorMag.shouldCorrect()) {
		// only correct heading
		for (int i = 0; i < Xe::n; i++) {
			if (i != Xe::rot_D) {
				dxe(i) = 0;
			}
		}

		setX(x);
		setP(_P + dP);
	}
}
