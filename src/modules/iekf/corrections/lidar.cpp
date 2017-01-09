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

void IEKF::correctLidar(const distance_sensor_s *msg)
{
	//ROS_INFO("correct lidar");
	//ROS_INFO("terrain stddev %10.4f",  double(sqrtf(_P(Xe::terrain_asl, Xe::terrain_asl))));
	// return if no new data
	float dt = 0;

	if (!_sensorLidar.ready(msg->timestamp, dt)) {
		return;
	}

	// attitude info
	Dcmf C_nb = Quaternion<float>(
			    _x(X::q_nb_0), _x(X::q_nb_1),
			    _x(X::q_nb_2), _x(X::q_nb_3));

	// abort if too large of an angle
	if (C_nb(2, 2) < 1e-1f) {
		ROS_INFO("lidar correction aborted, too large of an angle");
		return;
	}

	// expected measurement
	float agl = _x(X::asl) - _x(X::terrain_asl);
	float yh = agl / C_nb(2, 2);

	// measured airspeed
	float y = msg->current_distance;

	Vector<float, 1> r;
	r(0) = y - yh;

	// define R
	Matrix<float, Y_distance_down::n, Y_distance_down::n> R;
	R(Y_distance_down::d, Y_distance_down::d) = 2.5e-4f / dt;

	// define H
	// Note: this measurement is not invariant due to
	// rotation matrix
	Matrix<float, Y_distance_down::n, Xe::n> H;
	float x0 = 2 * agl / C_nb(2, 2) / C_nb(2, 2);
	float x1 = 1 / C_nb(2, 2);
	H(Y_distance_down::d, Xe::rot_N) = -C_nb(1, 2) * x0;
	H(Y_distance_down::d, Xe::rot_E) = C_nb(0, 2) * x0;
	H(Y_distance_down::d, Xe::asl) = x1;
	H(Y_distance_down::d, Xe::terrain_asl) = -x1;

	// kalman correction
	Vector<float, Xe::n> dxe;
	SquareMatrix<float, Xe::n> dP;
	_sensorLidar.kalmanCorrectCond(_P, H, R, r, dxe, dP);

	if (_sensorLidar.shouldCorrect()) {
		// only correct altitude states
		for (int i = 0; i < Xe::n; i++) {
			if (i != Xe::vel_D && i != Xe::asl && i != Xe::terrain_asl) {
				dxe(i) = 0;
			}
		}

		setX(applyErrorCorrection(dxe));
		setP(_P + dP);
	}
}
