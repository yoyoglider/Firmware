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

void IEKF::correctAccel(const sensor_combined_s *msg)
{
	//ROS_INFO("correct accel");
	// return if no new data
	float dt = 0;
	uint64_t timestamp = msg->timestamp + msg->accelerometer_timestamp_relative;

	if (!_sensorAccel.ready(timestamp, dt)) {
		return;
	}

	// measurement
	Vector3f y_b(
		msg->accelerometer_m_s2[0],
		msg->accelerometer_m_s2[1],
		msg->accelerometer_m_s2[2]);

	// don't correct if accelerating
	float relNormError = (Vector3f(y_b / _x(X::accel_scale)).norm()
			      - _g_n.norm()) / _g_n.norm();

	// calculate residual
	Quatf q_nb(_x(X::q_nb_0), _x(X::q_nb_1),
		   _x(X::q_nb_2), _x(X::q_nb_3));
	Vector3f y_g_n = q_nb.conjugate(y_b / _x(X::accel_scale));
	Vector3f r = y_g_n - _g_n;

	// define R
	// worst accel dir change is if accel is normal to gravity,
	// assume this and calculate angle covariance based on accel norm error
	Matrix<float, Y_accel::n, Y_accel::n> R;
	float cov = (accel_sigma_rw * accel_sigma_rw + relNormError * relNormError) / dt;
	R(Y_accel::accel_bX, Y_accel::accel_bX) = cov;
	R(Y_accel::accel_bY, Y_accel::accel_bY) = cov;
	R(Y_accel::accel_bZ, Y_accel::accel_bZ) = cov;

	// define H
	Matrix<float, Y_accel::n, Xe::n> H;
	Matrix3f tmp = _g_n.unit().hat() * 2;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			H(Y_accel::accel_bX + i, Xe::rot_N + j) = tmp(i, j);
		}
	}

	// kalman correction
	Vector<float, Xe::n> dxe;
	SquareMatrix<float, Xe::n> dP;
	_sensorAccel.kalmanCorrectCond(_P, H, R, r, dxe, dP);
	Vector<float, X::n> x = applyErrorCorrection(dxe);

	if (_sensorAccel.shouldCorrect()) {
		// only correct roll/pitch
		for (int i = 0; i < Xe::n; i++) {
			if (i != Xe::rot_N && i != Xe::rot_E) {
				dxe(i) = 0;
			}
		}

		setX(applyErrorCorrection(dxe));
		setP(_P + dP);
	}
}
