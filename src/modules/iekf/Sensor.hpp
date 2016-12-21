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

#pragma once

#include <matrix/math.hpp>

using namespace matrix;

/**
 * Class to hold data for sensors and manage
 * timestamps etc. There should be one for each
 * uniquely sent and timestamped message.
 */
class Sensor
{
public:
	Sensor(const char *name, float betaMax, float condMax) :
		_timestamp(0),
		_beta(0),
		_cond(1),
		_name(name),
		_betaMax(betaMax),
		_condMax(condMax)
	{
	}

	/**
	 * Check if ready and return dt
	 */
	bool ready(uint64_t timestampNew, float &dt)
	{
		// return if no new data
		if (timestampNew != _timestamp) {
			dt = (timestampNew - _timestamp) / 1.0e6f;

			if (dt < 0) {
				return false;
			}

			_timestamp = timestampNew;

		} else {
			return false;
		}

		return true;
	}
	inline void setBeta(float beta)
	{
		_beta = beta;
	}
	inline float getBeta()
	{
		return _beta;
	}
	inline bool faulted()
	{
		return _beta > _betaMax;
	}
	inline bool poorlyConditioned()
	{
		return _cond > _condMax;
	}
	inline bool shouldCorrect()
	{
		return !faulted() && !poorlyConditioned();
	}
	void setCorrectionInfo(float beta, float cond)
	{
		_beta = beta;
		_cond = cond;

		if (beta > _betaMax) {
			ROS_WARN("%s fault: beta %10.4f", _name, double(beta));
		}

		if (cond > _condMax) {
			ROS_WARN("%s poorly conditioned %10.4f", _name, double(cond));
		}
	}

	template<class Type, size_t n_x, size_t n_y>
	void kalmanCorrectCond(
		const SquareMatrix<Type, n_x> &P,
		const Matrix<Type, n_y, n_x> &H,
		const Matrix<Type, n_y, n_y> &R,
		const Matrix<Type, n_y, 1> &r,
		Vector<Type, n_x> &dx,
		SquareMatrix<Type, n_x> &dP
	)
	{
		dx.setZero();
		dP.setZero();
		SquareMatrix<Type, n_y> S = H * P * H.T() + R;
		Vector<Type, n_y> d = cholesky(S).diag();
		float cond = d.max() / d.min();
		SquareMatrix<Type, n_y> S_I = inv(S);
		float beta = (r.T() * (S_I * r))(0, 0) / BETA_TABLE[n_y];
		Matrix<Type, n_x, n_y> K = P * H.T() * S_I;
		setCorrectionInfo(beta, cond);
		dP = -K * H * P;
		dx = K * r;
	}

private:
	uint64_t _timestamp; // time of last read
	float _beta; // if > 1, indicates fault
	float _cond; // condition number, ratio of
	// eigen values of innovation matrix
	// S = H*P*H^T + R
	const char *_name;
	float _betaMax;
	float _condMax;
};

