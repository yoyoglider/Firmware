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

void IEKF::correctGps(const vehicle_gps_position_s *msg)
{
	//ROS_INFO("correct gps");
	// return if no new data
	float dt = 0;

	if (!_sensorGps.ready(msg->timestamp, dt)) {
		return;
	}

	// check for good gps signal
	if (msg->satellites_used < 6 || msg->fix_type < 3) {
		return;
	}

	double lat_deg = msg->lat * 1e-7;
	double lon_deg = msg->lon * 1e-7;
	float alt_m = msg->alt * 1e-3;

	// init global reference
	if (!_origin.xyInitialized()) {
		ROS_INFO("gps origin init lat: %12.6f deg lon: %12.6f deg",
			 double(lat_deg), double(lon_deg));
		_origin.xyInitialize(lat_deg, lon_deg, msg->timestamp);

		// init origin alt
		ROS_INFO("gps origin init alt %12.2f m", double(alt_m));
		_origin.altInitialize(alt_m, msg->timestamp);
		float deltaAlt = alt_m - _baroAsl;
		_baroOffset -= deltaAlt;
		_x(X::asl) = alt_m;

		// if we have no terrain data, guess we are on the ground
		if (!getTerrainValid()) {
			_x(X::terrain_asl) = alt_m;

		} else {
			_x(X::terrain_asl) += deltaAlt;
		}

	}

	// calculate residual
	float gps_pos_N = 0;
	float gps_pos_E = 0;
	_origin.globalToLocalXY(lat_deg, lon_deg, gps_pos_N, gps_pos_E);

	Vector<float, Y_gps::n> y;
	y(Y_gps::pos_N) = gps_pos_N;
	y(Y_gps::pos_E) = gps_pos_E;
	y(Y_gps::asl) = alt_m;
	y(Y_gps::vel_N) = msg->vel_n_m_s;
	y(Y_gps::vel_E) = msg->vel_e_m_s;
	y(Y_gps::vel_D) = msg->vel_d_m_s;

	//ROS_INFO("gps vx: %10.4f , iekf vx: %10.4f",
	//double(msg->vel_n_m_s), double(_x(X::vel_N)));

	Vector<float, Y_gps::n> yh;
	yh(Y_gps::pos_N) = _x(X::pos_N);
	yh(Y_gps::pos_E) = _x(X::pos_E);
	yh(Y_gps::asl) = _x(X::asl);
	yh(Y_gps::vel_N) = _x(X::vel_N);
	yh(Y_gps::vel_E) = _x(X::vel_E);
	yh(Y_gps::vel_D) = _x(X::vel_D);

	Vector<float, Y_gps::n> r = y - yh;

	// define R
	Matrix<float, Y_gps::n, Y_gps::n> R;
	float gps_xy_var = gps_xy_sigma_rw * gps_xy_sigma_rw / dt;
	float gps_z_var = gps_z_sigma_rw * gps_z_sigma_rw / dt;
	float gps_vxy_var = gps_vxy_sigma_rw * gps_vxy_sigma_rw / dt;
	float gps_vz_var = gps_vz_sigma_rw * gps_vz_sigma_rw / dt;
	R(Y_gps::pos_N, Y_gps::pos_N) = gps_xy_var;
	R(Y_gps::pos_E, Y_gps::pos_E) = gps_xy_var;
	R(Y_gps::asl, Y_gps::asl) = gps_z_var;
	R(Y_gps::vel_N, Y_gps::vel_N) = gps_vxy_var;
	R(Y_gps::vel_E, Y_gps::vel_E) = gps_vxy_var;
	R(Y_gps::vel_D, Y_gps::vel_D) = gps_vz_var;

	// define H
	Matrix<float, Y_gps::n, Xe::n> H;
	H(Y_gps::pos_N, Xe::pos_N) = 1;
	H(Y_gps::pos_E, Xe::pos_E) = 1;
	H(Y_gps::asl, Xe::asl) = 1;
	H(Y_gps::vel_N, Xe::vel_N) = 1;
	H(Y_gps::vel_E, Xe::vel_E) = 1;
	H(Y_gps::vel_D, Xe::vel_D) = 1;

	// kalman correction
	Vector<float, Xe::n> dxe;
	SquareMatrix<float, Xe::n> dP;
	_sensorGps.kalmanCorrectCond(_P, H, R, r, dxe, dP);

	if (_sensorGps.shouldCorrect()) {
		setX(applyErrorCorrection(dxe));
		setP(_P + dP);
	}
}
