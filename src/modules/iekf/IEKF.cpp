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

#include "IEKF.hpp"

IEKF::IEKF() :
	_nh(), // node handlke
	// sensors
	_sensorAccel("accel", 20, 50),
	_sensorMag("mag", 20, 50),
	_sensorBaro("baro", 20, 50),
	_sensorGps("gps", 20, 50),
	_sensorAirspeed("airspeed", 20 , 50),
	_sensorFlow("flow", 20, 50),
	_sensorSonar("sonar", 20, 50),
	_sensorLidar("lidar", 20, 50),
	_sensorVision("vision", 20, 50),
	_sensorMocap("mocap", 20, 50),
	// subscriptions
	_subImu(),
	_subGps(),
	_subAirspeed(),
	_subFlow(),
	_subDistance(),
	_subVision(),
	_subMocap(),
	// publications
	_pubAttitude(),
	_pubLocalPosition(),
	_pubGlobalPosition(),
	_pubControlState(),
	_pubEstimatorStatus(),
	// data
	_x0(),
	_xMin(),
	_xMax(),
	_P0Diag(),
	_x(),
	_P(),
	_u(),
	_g_n(0, 0, -9.8),
	_origin(),
	_baroAsl(0),
	_baroOffset(0),
	_gpsUSec(0),
	_magDeclDeg(0),
	_magInclDeg(0),
	_attitudeInitialized(false)
{
	// for quaterinons we bound at 2
	// so it has a chance to
	// do a renormalization first
	_xMin(X::q_nb_0) = -1.1;
	_xMin(X::q_nb_1) = -1.1;
	_xMin(X::q_nb_2) = -1.1;
	_xMin(X::q_nb_3) = -1.1;
	_xMin(X::vel_N) = -100;
	_xMin(X::vel_E) = -100;
	_xMin(X::vel_D) = -100;
	_xMin(X::gyro_bias_bX) = -0.1;
	_xMin(X::gyro_bias_bY) = -0.1;
	_xMin(X::gyro_bias_bZ) = -0.1;
	_xMin(X::accel_scale) = 0.8;
	_xMin(X::pos_N) = -1e9;
	_xMin(X::pos_E) = -1e9;
	_xMin(X::asl) = -1e9;
	_xMin(X::terrain_asl) = -1e6;
	_xMin(X::baro_bias) = -1e6;
	_xMin(X::wind_N) = -100;
	_xMin(X::wind_E) = -100;
	_xMin(X::wind_D) = -100;

	_xMax(X::q_nb_0) = 1.1;
	_xMax(X::q_nb_1) = 1.1;
	_xMax(X::q_nb_2) = 1.1;
	_xMax(X::q_nb_3) = 1.1;
	_xMax(X::vel_N) = 100;
	_xMax(X::vel_E) = 100;
	_xMax(X::vel_D) = 100;
	_xMax(X::gyro_bias_bX) = 0.1;
	_xMax(X::gyro_bias_bY) = 0.1;
	_xMax(X::gyro_bias_bZ) = 0.1;
	_xMax(X::accel_scale) = 1.5;
	_xMax(X::pos_N) = 1e9;
	_xMax(X::pos_E) = 1e9;
	_xMax(X::asl) = 1e9;
	_xMax(X::terrain_asl) = 1e6;
	_xMax(X::baro_bias) = 1e6;
	_xMax(X::wind_N) = 100;
	_xMax(X::wind_E) = 100;
	_xMax(X::wind_D) = 100;

	// initialize state
	_x0(X::q_nb_0) = 1;
	_x0(X::q_nb_1) = 0;
	_x0(X::q_nb_2) = 0;
	_x0(X::q_nb_3) = 0;
	_x0(X::accel_scale) = 1;
	setX(_x0);

	// initialize covariance
	_P0Diag(Xe::rot_N) = 0;
	_P0Diag(Xe::rot_E) = 0;
	_P0Diag(Xe::rot_D) = 0;
	_P0Diag(Xe::vel_N) = 0;
	_P0Diag(Xe::vel_E) = 0;
	_P0Diag(Xe::vel_D) = 0;
	_P0Diag(Xe::gyro_bias_N) = 1e-6f;
	_P0Diag(Xe::gyro_bias_E) = 1e-6f;
	_P0Diag(Xe::gyro_bias_D) = 1e-6f;
	_P0Diag(Xe::accel_scale) = 1e-6f;
	_P0Diag(Xe::pos_N) = 0;
	_P0Diag(Xe::pos_E) = 0;
	_P0Diag(Xe::asl) = 0;
	_P0Diag(Xe::terrain_asl) = 0;
	_P0Diag(Xe::baro_bias) = 0;
	_P0Diag(Xe::wind_N) = 0;
	_P0Diag(Xe::wind_E) = 0;
	_P0Diag(Xe::wind_D) = 0;
	setP(diag(_P0Diag));


	// magnetic field values for SITL
	// TODO load from GPS or params
	_magDeclDeg = magDeclDeg;
	_magInclDeg = magInclDeg;
}

void IEKF::init()
{
	// subscriptions
	_subImu = _nh.subscribe("sensor_combined", 0, &IEKF::callbackImu, this);
	_subGps = _nh.subscribe("vehicle_gps_position", 0, &IEKF::correctGps, this);
	_subAirspeed = _nh.subscribe("airspeed", 0, &IEKF::correctAirspeed, this);
	_subFlow = _nh.subscribe("optical_flow", 0, &IEKF::correctFlow, this);
	_subDistance = _nh.subscribe("distance_sensor", 0, &IEKF::callbackDistance, this);
	_subVision = _nh.subscribe("vision_position_estimate", 0, &IEKF::correctVision, this);
	_subMocap = _nh.subscribe("att_pos_mocap", 0, &IEKF::correctMocap, this);
	// publications
	_pubAttitude = _nh.advertise<vehicle_attitude_s>("vehicle_attitude", 0);
	_pubLocalPosition = _nh.advertise<vehicle_local_position_s>("vehicle_local_position", 0);
	_pubGlobalPosition = _nh.advertise<vehicle_global_position_s>("vehicle_global_position", 0);
	_pubControlState = _nh.advertise<control_state_s>("control_state", 0);
	_pubEstimatorStatus = _nh.advertise<estimator_status_s>("estimator_status", 0);
}

Vector<float, X::n> IEKF::dynamics(float t, const Vector<float, X::n> &x, const Vector<float, U::n> &u)
{
	Quatf q_nb(x(X::q_nb_0), x(X::q_nb_1), x(X::q_nb_2), x(X::q_nb_3));
	Vector3f a_b(u(U::accel_bX), u(U::accel_bY), u(U::accel_bZ));
	Vector3f a_n = q_nb.conjugate(a_b / x(X::accel_scale));
	Vector3f as_n = a_n - _g_n;
	Vector3f gyro_bias_b(x(X::gyro_bias_bX), x(X::gyro_bias_bY), x(X::gyro_bias_bZ));
	Vector3f omega_nb_b(u(U::omega_nb_bX), u(U::omega_nb_bY), u(U::omega_nb_bZ));
	Vector3f omega_nb_b_corrected = omega_nb_b - gyro_bias_b;
	Quatf dq_nb = q_nb * Quatf(0, omega_nb_b_corrected(0),
				   omega_nb_b_corrected(1), omega_nb_b_corrected(2)) * 0.5f;

	Vector<float, X::n> dx;

	if (getAttitudeValid()) {
		dx(X::q_nb_0) = dq_nb(0);
		dx(X::q_nb_1) = dq_nb(1);
		dx(X::q_nb_2) = dq_nb(2);
		dx(X::q_nb_3) = dq_nb(3);
	}

	//ROS_INFO("as_n");
	//as_n.print();
	//ROS_INFO("a_b");
	//a_b.print();
	//ROS_INFO("a_n");
	//a_n.print();

	if (getVelocityValid()) {
		dx(X::vel_N) = as_n(0);
		dx(X::vel_E) = as_n(1);
		dx(X::vel_D) = as_n(2);
	}

	dx(X::gyro_bias_bX) = -_x(X::gyro_bias_bX) / gyro_correlation_time;
	dx(X::gyro_bias_bY) = -_x(X::gyro_bias_bY) / gyro_correlation_time;
	dx(X::gyro_bias_bZ) = -_x(X::gyro_bias_bZ) / gyro_correlation_time;
	dx(X::accel_scale) = 0; // TODO

	if (getPositionValid()) {
		dx(X::pos_N) = x(X::vel_N);
		dx(X::pos_E) = x(X::vel_E);
		dx(X::asl) = -x(X::vel_D);
	}

	// want terrain dynamics to be static, so when out of range it keeps
	// last estimate and doesn't decay
	dx(X::terrain_asl) = 0;
	dx(X::baro_bias) = -x(X::baro_bias) / baro_correlation_time;
	dx(X::wind_N) = -_x(X::wind_N) / wind_correlation_time;
	dx(X::wind_E) = -_x(X::wind_E) / wind_correlation_time;
	dx(X::wind_D) = -_x(X::wind_D) / wind_correlation_time;
	return dx;
}

void IEKF::callbackImu(const sensor_combined_s *msg)
{
	//ROS_INFO("imu callback");
	_u(U::omega_nb_bX) = msg->gyro_rad[0];
	_u(U::omega_nb_bY) = msg->gyro_rad[1];
	_u(U::omega_nb_bZ) = msg->gyro_rad[2];
	_u(U::accel_bX) = msg->accelerometer_m_s2[0];
	_u(U::accel_bY) = msg->accelerometer_m_s2[1];
	_u(U::accel_bZ) = msg->accelerometer_m_s2[2];

	if (_attitudeInitialized) {
		// predict driven by gyro callback
		if (msg->gyro_integral_dt > 0) {
			predict(msg->gyro_integral_dt);
		};

		// correct  if new data
		correctAccel(msg);

		correctMag(msg);

		correctBaro(msg);

	} else {
		initializeAttitude(msg);
	}
}

void IEKF::initializeAttitude(const sensor_combined_s *msg)
{
	// return if no new mag data
	float dt = 0;
	uint64_t timestamp = msg->timestamp + msg->magnetometer_timestamp_relative;

	if (!_sensorMag.ready(timestamp, dt)) {
		return;
	}

	// TODO handle freefall initialization (where accel norm is small and can't be used)
	Vector3f accel = Vector3f(msg->accelerometer_m_s2[0],
				  msg->accelerometer_m_s2[1],
				  msg->accelerometer_m_s2[2]);

	Vector3f mag = Vector3f(
			       msg->magnetometer_ga[0],
			       msg->magnetometer_ga[1],
			       msg->magnetometer_ga[2]).unit();

	Vector3f k = -accel.unit(); // then only acceleration is likely opposing gravity
	Vector3f i = (mag - mag.dot(k) * k).unit(); // project onto NE plane
	Vector3f j = k.cross(i).unit();

	Dcmf C_nb;
	C_nb.setRow(0, i);
	C_nb.setRow(1, j);
	C_nb.setRow(2, k);

	// account for magnetic declination
	//Quatf q_nb = Dcmf(Dcmf(AxisAnglef(Vector3f(0, 0, 1), deg2radf*_magDeclDeg))*C_nb);
	Quatf q_nb = C_nb;

	_x(X::q_nb_0) = q_nb(0);
	_x(X::q_nb_1) = q_nb(1);
	_x(X::q_nb_2) = q_nb(2);
	_x(X::q_nb_3) = q_nb(3);
	_attitudeInitialized = true;

	Eulerf euler = q_nb;
	ROS_INFO("initial euler angles (deg) roll: %10.4f pitch: %10.4f yaw: %10.4f",
		 double(rad2degf * euler(0)),
		 double(rad2degf * euler(1)),
		 double(rad2degf * euler(2)));
}

void IEKF::callbackDistance(const distance_sensor_s *msg)
{
	// require attitude to be initialized
	if (!_attitudeInitialized) {
		return;
	}

	// if not pointing down (roll 180 by convention), do not use
	if (msg->orientation != 8) {
		ROS_INFO("distance sensor wrong orientation %d", msg->orientation);
		return;
	}

	// if above max distance/ <= 0, out of range
	if (msg->current_distance > msg->max_distance ||
	    msg->current_distance < msg->min_distance) {
		return;
	}

	// if below 0, don't correct and warn
	if (msg->current_distance < 0) {
		ROS_WARN("distance below 0");
		return;
	}

	// call correct correction function based on type
	if (msg->type == distance_sensor_s::MAV_DISTANCE_SENSOR_ULTRASOUND) {
		correctSonar(msg);

	} else if (msg->type == distance_sensor_s::MAV_DISTANCE_SENSOR_LASER) {
		correctLidar(msg);
	}
}

void IEKF::correctMocap(const att_pos_mocap_s *msg)
{
	//ROS_INFO("correct mocap");
	// return if no new data
	float dt = 0;

	if (!_sensorMocap.ready(msg->timestamp, dt)) {
		return;
	}
}

void IEKF::predict(float dt)
{
	// normalize quaternions if needed
	Quatf q_nb(
		_x(X::q_nb_0), _x(X::q_nb_1),
		_x(X::q_nb_2), _x(X::q_nb_3));

	if (fabsf(Quatf(
			  _x(X::q_nb_0), _x(X::q_nb_1),
			  _x(X::q_nb_2), _x(X::q_nb_3)).norm() - 1.0f) > 1e-3f) {
		ROS_INFO("normalizing quaternion, norm was %6.4f\n",
			 double(q_nb.norm()));
		q_nb.normalize();
		_x(X::q_nb_0) = q_nb(0);
		_x(X::q_nb_1) = q_nb(1);
		_x(X::q_nb_2) = q_nb(2);
		_x(X::q_nb_3) = q_nb(3);
	}

	// rotation rate
	Vector3f omega_nb_b(
		_u(U::omega_nb_bX), _u(U::omega_nb_bY), _u(U::omega_nb_bZ));
	Vector3f gyro_bias_b(
		_x(X::gyro_bias_bX), _x(X::gyro_bias_bY), _x(X::gyro_bias_bZ));
	Vector3f omega_nb_b_corrected = omega_nb_b - gyro_bias_b;

	// define process noise matrix
	Matrix<float, Xe::n, Xe::n> Q;
	float gyro_Q_rrw = gyro_sigma_rrw * gyro_sigma_rrw;
	float gyro_Q_rw = gyro_sigma_rw * gyro_sigma_rw;
	float accel_Q_rw = accel_sigma_rw * accel_sigma_rw;
	float baro_Q_rrw = baro_sigma_rrw * baro_sigma_rrw;
	Q(Xe::rot_N, Xe::rot_N) = gyro_Q_rw;
	Q(Xe::rot_E, Xe::rot_E) = gyro_Q_rw;
	Q(Xe::rot_D, Xe::rot_D) = gyro_Q_rw;
	Q(Xe::vel_N, Xe::vel_N) = accel_Q_rw;
	Q(Xe::vel_E, Xe::vel_E) = accel_Q_rw;
	Q(Xe::vel_D, Xe::vel_D) = accel_Q_rw;
	Q(Xe::gyro_bias_N, Xe::gyro_bias_N) = gyro_Q_rrw;
	Q(Xe::gyro_bias_E, Xe::gyro_bias_E) = gyro_Q_rrw;
	Q(Xe::gyro_bias_D, Xe::gyro_bias_D) = gyro_Q_rrw;
	Q(Xe::accel_scale, Xe::accel_scale) = 1e-6f;
	Q(Xe::pos_N, Xe::pos_N) = 1e-1f;
	Q(Xe::pos_E, Xe::pos_E) = 1e-1f;
	Q(Xe::asl, Xe::asl) = 1e-1f;
	Q(Xe::terrain_asl, Xe::terrain_asl) = 1e-1f;
	Q(Xe::baro_bias, Xe::baro_bias) = baro_Q_rrw;
	Q(Xe::wind_N, Xe::wind_N) = 1e-2f;
	Q(Xe::wind_E, Xe::wind_E) = 1e-2f;
	Q(Xe::wind_D, Xe::wind_D) = 1e-2f;

	// define A matrix
	Matrix<float, Xe::n, Xe::n> A;

	// derivative of rotation error is -0.5 * gyro bias
	A(Xe::rot_N, Xe::Xe::gyro_bias_N) = -0.5;
	A(Xe::rot_E, Xe::Xe::gyro_bias_E) = -0.5;
	A(Xe::rot_D, Xe::Xe::gyro_bias_D) = -0.5;

	// derivative of velocity
	Vector3f a_b(_u(U::accel_bX), _u(U::accel_bY), _u(U::accel_bZ));
	Vector3f J_a_n = q_nb.conjugate(a_b / _x(X::accel_scale));
	Matrix<float, 3, 3> a_tmp = -J_a_n.hat() * 2;

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			A(Xe::vel_N + i, Xe::rot_N + j) = a_tmp(i, j);
		}

		A(Xe::vel_N + i, Xe::accel_scale) = -J_a_n(i);
	}

	// derivative of gyro bias
	Vector3f J_omega_n = q_nb.conjugate(omega_nb_b_corrected);
	Matrix<float, 3, 3> g_tmp = J_omega_n.hat();

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			A(Xe::gyro_bias_N + i, Xe::gyro_bias_N + j) = g_tmp(i, j);
		}
	}

	// derivative of position is velocity
	A(Xe::pos_N, Xe::vel_N) = 1;
	A(Xe::pos_E, Xe::vel_E) = 1;
	A(Xe::asl, Xe::vel_D) = -1;

	// derivative of terrain alt is zero

	// derivative of baro bias
	A(Xe::baro_bias, Xe::baro_bias) = -1.0f / baro_correlation_time;

	// derivative of gyro bias
	A(Xe::gyro_bias_N, Xe::gyro_bias_N) = -1 / gyro_correlation_time;
	A(Xe::gyro_bias_E, Xe::gyro_bias_E) = -1 / gyro_correlation_time;
	A(Xe::gyro_bias_D, Xe::gyro_bias_D) = -1 / gyro_correlation_time;

	//ROS_INFO("A:");
	//for (int i=0;i<Xe::n; i++) {
	//for (int j=0;j<Xe::n; j++) {
	//printf("%10.3f, ", double(A(i, j)));
	//}
	//printf("\n");
	//}

	// continuous time kalman filter prediction
	// integrate runge kutta 4th order
	// TODO move rk4 algorithm to matrixlib
	// https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods
	float h = dt;
	Vector<float, X::n> k1, k2, k3, k4;
	k1 = dynamics(0, _x, _u);
	k2 = dynamics(h / 2, _x + k1 * h / 2, _u);
	k3 = dynamics(h / 2, _x + k2 * h / 2, _u);
	k4 = dynamics(h, _x + k3 * h, _u);
	Vector<float, X::n> dx = (k1 + k2 * 2 + k3 * 2 + k4) * (h / 6);

	//ROS_INFO("dx predict \n");
	//dx.print();
	setX(_x + dx);

	// propgate covariance using euler integration
	Matrix<float, Xe::n, Xe::n> dP = (A * _P + _P * A.T() + Q) * dt;
	setP(_P + dP);

	//ROS_INFO("P:");
	//_P.diag().print();
}

Vector<float, X::n> IEKF::applyErrorCorrection(const Vector<float, Xe::n> &d_xe)
{
	Quatf q_nb(_x(X::q_nb_0), _x(X::q_nb_1), _x(X::q_nb_2), _x(X::q_nb_3));
	Quatf d_q_nb = Quatf(0.0f,
			     d_xe(Xe::rot_N), d_xe(Xe::rot_E), d_xe(Xe::rot_D)) * q_nb;

	//ROS_INFO("d_q_nb");
	//d_q_nb.print();
	Vector3f d_gyro_bias_b = q_nb.conjugate_inversed(
					 Vector3f(d_xe(Xe::gyro_bias_N),
							 d_xe(Xe::gyro_bias_E),
							 d_xe(Xe::gyro_bias_D)));

	// linear term correction is the same
	// as the error correction
	Vector<float, X::n> x = _x;
	x(X::q_nb_0) += d_q_nb(0);
	x(X::q_nb_1) += d_q_nb(1);
	x(X::q_nb_2) += d_q_nb(2);
	x(X::q_nb_3) += d_q_nb(3);
	x(X::vel_N) += d_xe(Xe::vel_N);
	x(X::vel_E) += d_xe(Xe::vel_E);
	x(X::vel_D) += d_xe(Xe::vel_D);
	x(X::gyro_bias_bX) += d_gyro_bias_b(0);
	x(X::gyro_bias_bY) += d_gyro_bias_b(1);
	x(X::gyro_bias_bZ) +=  d_gyro_bias_b(2);
	x(X::accel_scale) += _x(X::accel_scale) * d_xe(Xe::accel_scale);
	x(X::pos_N) += d_xe(Xe::pos_N);
	x(X::pos_E) += d_xe(Xe::pos_E);
	x(X::asl) += d_xe(Xe::asl);
	x(X::terrain_asl) += d_xe(Xe::terrain_asl);
	x(X::baro_bias) += d_xe(Xe::baro_bias);
	x(X::wind_N) += d_xe(Xe::wind_N);
	x(X::wind_E) += d_xe(Xe::wind_E);
	x(X::wind_D) += d_xe(Xe::wind_D);
	return x;
}

void IEKF::setP(const SquareMatrix<float, Xe::n> &P)
{
	_P = P;

	for (int i = 0; i < Xe::n; i++) {
		// only operate on upper triangle, then copy to lower

		// don't allow NaN or large numbers
		for (int j = 0; j <= i; j++) {
			if (!PX4_ISFINITE(_P(i, j))) {
				ROS_WARN("P(%d, %d) NaN, resetting", i, j);

				if (i == j) {
					_P(i, j) = _P0Diag(i);

				} else {
					_P(i, j) = 0;
				}
			}

			if (_P(i, j) > 1e6f) {
				// upper bound
				_P(i, j) = 1e6f;
			}
		}

		// force non-negative diagonal
		if (_P(i, i) < 0) {
			ROS_WARN("P(%d, %d) < 0, setting to P0 val", i, i, double(0));
			_P(i, i) = _P0Diag(i);
		}

		// force symmetry, copy uppper triangle to lower
		for (int j = 0; j < i; j++) {
			_P(j, i) = _P(i, j);
		}
	}
}

void IEKF::setX(const Vector<float, X::n> &x)
{
	// set private state
	_x = x;

	// normalize quaternion
	Quatf q_nb(
		_x(X::q_nb_0), _x(X::q_nb_1),
		_x(X::q_nb_2), _x(X::q_nb_3));

	if (fabsf(q_nb.norm() - 1.0f) > 1e-3f) {
		ROS_INFO("normalizing quaternion, norm was %6.4f\n", double(q_nb.norm()));
		q_nb.normalize();
		_x(X::q_nb_0) = q_nb(0);
		_x(X::q_nb_1) = q_nb(1);
		_x(X::q_nb_2) = q_nb(2);
		_x(X::q_nb_3) = q_nb(3);
	}

	// saturate
	for (int i = 0; i < X::n; i++) {
		if (!PX4_ISFINITE(_x(i))) {
			ROS_WARN("x(%d) NaN, setting to %10.4f", i, double(_x0(i)));
			_x(i) = _x0(i);
		}

		if (_x(i) < _xMin(i)) {
			ROS_WARN("x(%d) < lower bound, saturating", i);
			_x(i) = _xMin(i);

		} else if (_x(i) > _xMax(i)) {
			ROS_WARN("x(%d) > upper bound, saturating", i);
			_x(i) = _xMax(i);
		}
	}
}

void IEKF::publish()
{
	//ROS_INFO("x:");
	//_x.print();

	//ROS_INFO("P:");
	//_P.diag().print();

	float eph = sqrt(_P(Xe::pos_N, Xe::pos_N) + _P(Xe::pos_E, Xe::pos_E));
	float epv = _P(Xe::asl, Xe::asl);
	Quatf q_nb(
		_x(X::q_nb_0), _x(X::q_nb_1),
		_x(X::q_nb_2), _x(X::q_nb_3));
	Euler<float> euler_nb = q_nb;
	Vector3f a_b(_u(U::accel_bX), _u(U::accel_bY), _u(U::accel_bZ));
	Vector3f a_n = q_nb.conjugate(a_b / _x(X::accel_scale));
	ros::Time now = ros::Time::now();

	// predicted airspeed
	Vector3f wind_n(_x(X::wind_N), _x(X::wind_E), _x(X::wind_D));
	Vector3f vel_n(_x(X::vel_N), _x(X::vel_E), _x(X::vel_D));
	Vector3f wind_rel_b = q_nb.conjugate_inversed(wind_n - vel_n);
	float airspeed = -wind_rel_b(0); // body -x component aligned with pitot tube

	//bool attitudeValid = sqrtf(_P(Xe::rot_N, Xe::rot_N)
	//+ _P(Xe::rot_E, Xe::rot_E)
	//+ _P(Xe::rot_D, Xe::rot_D)) < 0.1f;
	bool attitudeValid = getAttitudeValid();
	bool velocityValid = getVelocityValid();
	bool positionValid = getPositionValid();
	bool terrainValid = getTerrainValid();

	// publish attitude
	if (attitudeValid) {
		vehicle_attitude_s msg = {};
		msg.timestamp = now.toNSec() / 1e3;
		msg.q[0] = _x(X::q_nb_0);
		msg.q[1] = _x(X::q_nb_1);
		msg.q[2] = _x(X::q_nb_2);
		msg.q[3] = _x(X::q_nb_3);
		msg.rollspeed = _u(U::omega_nb_bX) - _x(X::gyro_bias_bX);
		msg.pitchspeed = _u(U::omega_nb_bY) - _x(X::gyro_bias_bY);
		msg.yawspeed = _u(U::omega_nb_bZ) - _x(X::gyro_bias_bZ);
		_pubAttitude.publish(msg);
	}

	// publish local position
	if (_origin.xyInitialized() && _origin.altInitialized() && velocityValid) {
		vehicle_local_position_s msg = {};
		msg.timestamp = now.toNSec() / 1e3;
		msg.xy_valid = positionValid;
		msg.z_valid = positionValid;
		msg.v_xy_valid = velocityValid;
		msg.v_z_valid = velocityValid;
		msg.x = _x(X::pos_N);
		msg.y = _x(X::pos_E);
		// TODO make agl or relative alt to origin an option ?
		msg.z = -(_x(X::asl) - _origin.getAlt());
		msg.delta_xy[0] = 0;
		msg.delta_xy[1] = 0;
		msg.delta_z = 0;
		msg.vx = _x(X::vel_N);
		msg.vy = _x(X::vel_E);
		msg.vz = _x(X::vel_D);
		msg.delta_vxy[0] = 0;
		msg.delta_vxy[1] = 0;
		msg.delta_vz = 0;
		msg.xy_reset_counter = 0;
		msg.z_reset_counter = 0;
		msg.vxy_reset_counter = 0;
		msg.vz_reset_counter = 0;
		msg.yaw = euler_nb(2);
		msg.xy_global = _origin.xyInitialized();
		msg.z_global = _origin.altInitialized();
		msg.ref_timestamp = _origin.getXYTimestamp();
		msg.ref_lat = _origin.getLatDeg();
		msg.ref_lon = _origin.getLonDeg();
		msg.ref_alt = _origin.getAlt();
		msg.dist_bottom = _x(X::asl) - _x(X::terrain_asl);
		msg.dist_bottom_rate = -_x(X::vel_D);
		msg.surface_bottom_timestamp = now.toNSec() / 1e3;
		msg.dist_bottom_valid = terrainValid;
		msg.eph = eph;
		msg.epv = epv;
		_pubLocalPosition.publish(msg);
	}

	// publish global position
	if (_origin.xyInitialized() && _origin.altInitialized() && velocityValid) {
		double lat_deg = 0;
		double lon_deg = 0;
		_origin.northEastToLatLon(_x(X::pos_N), _x(X::pos_E), lat_deg, lon_deg);
		//ROS_INFO("alt %10.4f m", double(alt_m));
		vehicle_global_position_s msg = {};
		msg.timestamp = now.toNSec() / 1e3;
		msg.time_utc_usec = _gpsUSec;
		msg.lat = lat_deg;
		msg.lon = lon_deg;
		msg.alt = _x(X::asl);
		msg.delta_lat_lon[0] = 0;
		msg.delta_lat_lon[1] = 0;
		msg.delta_alt = 0;
		msg.lat_lon_reset_counter = 0;
		msg.alt_reset_counter = 0;
		msg.vel_n = _x(X::vel_N);
		msg.vel_e = _x(X::vel_E);
		msg.vel_d = _x(X::vel_D);
		msg.yaw = euler_nb(2);
		msg.eph = eph;
		msg.epv = epv;
		msg.terrain_alt = _x(X::terrain_asl);
		msg.terrain_alt_valid = terrainValid;
		msg.dead_reckoning = false;
		msg.pressure_alt = _baroAsl;
		_pubGlobalPosition.publish(msg);
	}

	// publish control state
	{
		// specific acceleration
		control_state_s msg = {};
		msg.timestamp = now.toNSec() / 1e3;
		msg.x_acc = a_n(0);
		msg.y_acc = a_n(1);
		msg.z_acc = a_n(2);
		msg.x_vel = _x(X::vel_N);
		msg.y_vel = _x(X::vel_E);
		msg.z_vel = _x(X::vel_D);
		msg.x_pos = _x(X::pos_N);
		msg.y_pos = _x(X::pos_E);
		msg.z_pos = -(_x(X::asl) - _origin.getAlt());
		msg.airspeed = airspeed;
		msg.airspeed_valid = true;
		msg.vel_variance[0] = _P(Xe::vel_N, Xe::vel_N);
		msg.vel_variance[1] = _P(Xe::vel_E, Xe::vel_E);
		msg.vel_variance[2] = _P(Xe::vel_D, Xe::vel_D);
		msg.pos_variance[0] = _P(Xe::pos_N, Xe::pos_N);
		msg.pos_variance[1] = _P(Xe::pos_E, Xe::pos_E);
		msg.pos_variance[2] = _P(Xe::asl, Xe::asl);
		msg.q[0] = _x(X::q_nb_0);
		msg.q[1] = _x(X::q_nb_1);
		msg.q[2] = _x(X::q_nb_2);
		msg.q[3] = _x(X::q_nb_3);
		msg.delta_q_reset[0] = 0;
		msg.delta_q_reset[1] = 0;
		msg.delta_q_reset[2] = 0;
		msg.delta_q_reset[3] = 0;
		msg.quat_reset_counter = 0;
		msg.roll_rate = _u(U::omega_nb_bX) - _x(X::gyro_bias_bX);
		msg.pitch_rate = _u(U::omega_nb_bY) - _x(X::gyro_bias_bY);
		msg.yaw_rate = _u(U::omega_nb_bZ) - _x(X::gyro_bias_bZ);
		msg.horz_acc_mag = 0;
		_pubControlState.publish(msg);
	}

	// estimator status
	{
		estimator_status_s msg = {};
		msg.timestamp = now.toNSec() / 1e3;
		msg.vibe[0] = 0; // TODO
		msg.vibe[1] = 0; // TODO
		msg.vibe[2] = 0; // TODO
		msg.n_states = X::n;

		for (int i = 0; i < X::n; i++) {
			msg.states[i] = _x(i);
		}

		for (int i = 0; i < Xe::n; i++) {
			msg.covariances[i] = _P(i, i);
		}

		// XXX
		// this isn't really general and is
		// tailored to EKF2 so just dumping
		// data in the best field names available
		msg.gps_check_fail_flags = 0; // TODO
		msg.control_mode_flags = 0; // TODO
		msg.filter_fault_flags = 0; // TODO
		msg.pos_horiz_accuracy = eph;
		msg.pos_vert_accuracy = epv;
		msg.innovation_check_flags = 0; // TODO
		msg.mag_test_ratio = _sensorMag.getBeta();
		msg.vel_test_ratio = _sensorFlow.getBeta();
		msg.pos_test_ratio = _sensorGps.getBeta();
		msg.hgt_test_ratio = _sensorAccel.getBeta();
		msg.tas_test_ratio = _sensorAirspeed.getBeta();
		msg.hagl_test_ratio = _sensorLidar.getBeta();
		msg.solution_status_flags = 0; // TODO
		_pubEstimatorStatus.publish(msg);
	}
}
