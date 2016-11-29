#ifndef _BRDF_MERL_
#define _BRDF_MERL_

#include <iostream>
#include <cassert>
#include <cstdlib>
#include <math.h>
#include <vector>

#include "glm/glm.hpp"
using namespace glm;

#include "brdf.h"
#include "MERL.h"

class BrdfMERL : public Brdf
{
public:
	float eval(const vec3& V, const vec3& L, const float alpha, float& pdf) const
	{
		if(V.z <= 0)
		{
			pdf = 0;
			return 0;
		}

		// Undo formatting performed in `sample` function.
		float l_theta = atanf(L.y / L.x);
		float l_phi = acosf(L.z);
		unsigned int channel = static_cast<unsigned int>(alpha);
		unsigned int v_theta_index = static_cast<unsigned int>(round((acosf(V.z) * float(_N - 1)) / (0.5 * M_PI)));
		unsigned int l_theta_index = static_cast<unsigned int>(round(l_theta * float(_NumSamples - 1) / (0.5 * M_PI)));
		unsigned int l_phi_index = static_cast<unsigned int>(round(l_phi * float(_NumSamples - 1) / (2.0 * M_PI)));

		pdf = _pdf[channel][v_theta_index][l_phi_index][l_theta_index];
		return _brdf[channel][v_theta_index][l_phi_index][l_theta_index];
	}

	vec3 sample(const vec3& V, const float alpha, const float U1, const float U2) const
	{
		unsigned int v_theta_index = static_cast<unsigned int>(round((acosf(V.z) * float(_N - 1)) / (0.5 * M_PI)));
		unsigned int l_phi_index = static_cast<unsigned int>(round(float(_NumSamples - 1) * U1)); // assume phi_light is uniformly distributed
		unsigned int l_theta_index = 0;

		int channel = static_cast<int>(alpha);
		for (int i = 0; i < _NumSamples; ++i) // argmin (cdf[i] == rand)
			if (abs(_cdf[channel][v_theta_index][l_phi_index][i] - U2) < 0.00001)
			{
				l_theta_index = i;
				break;
			}

		// Mainly transforming to this format to maintain the I/O expected by the original program.
		float o_theta = l_theta_index * 0.5 * M_PI / static_cast<float>(_NumSamples - 1);
		float o_phi = l_phi_index * 2.0 * M_PI / static_cast<float>(_NumSamples - 1);
		vec3 L = vec3(cosf(o_theta) * sinf(o_phi), sinf(o_theta) * sinf(o_phi), cosf(o_phi));
		return L;
	}

	void load(const char* filename, const int NumSamples, const int N)
	{
		bool result = read_brdf(filename, _raw_brdf_data);
		if (!result)
		{
			cout << "Error reading " << filename << endl;
			exit(-1);
		}

		_NumSamples = NumSamples; _N = N;

		_brdf.resize(3, vector<vector<vector<double> > >(N, vector<vector<double> >(NumSamples, vector<double>(NumSamples, 0))));
		_pdf.resize(3, vector<vector<vector<double> > >(N, vector<vector<double> >(NumSamples, vector<double>(NumSamples, 0))));
		_cdf.resize(3, vector<vector<vector<double> > >(N, vector<vector<double> >(NumSamples, vector<double>(NumSamples, 0))));

		// Load BRDF data per channel.
		for (int i = 0; i < N; ++i)
		{
			double theta_in = i * 0.5 * M_PI / static_cast<double>(N - 1);
			for (int j = 0; j < N; ++j)
			{
				double phi_in = j * 2.0 * M_PI / static_cast<double>(N - 1);
				for (int k = 0; k < NumSamples; ++k)
				{
					double phi_out = k * 2.0 * M_PI / static_cast<double>(NumSamples - 1);
					for (int l = 0; l < NumSamples; ++l)
					{
						double theta_out = l * 0.5 * M_PI / static_cast<double>(NumSamples - 1);

						double red, green, blue;
						lookup_brdf_val(_raw_brdf_data, theta_in, 0.0, theta_out, phi_out, red, green, blue);
						_brdf[0][i][k][l] += blue;
						_brdf[1][i][k][l] += green;
						_brdf[2][i][k][l] += red;
						//cout << "brdf: " << _brdf[0][i][k][l] << endl;
					}
					//assert(false);
				}
			}

			// average view phi contributions
			for (int k = 0; k < NumSamples; ++k)
				for (int l = 0; l < NumSamples; ++l)
				{
					_brdf[0][i][k][l] /= static_cast<double>(N);
					_brdf[1][i][k][l] /= static_cast<double>(N);
					_brdf[2][i][k][l] /= static_cast<double>(N);
				}
		}

		// Construct PDFs per channel.
		for (int c = 0; c < 3; ++c) // channel
		for (int i = 0; i < N; ++i)
		{
			// Copy brdf table to pdf table for normalization.
			_pdf[c][i] = _brdf[c][i];

			// normalize => integral(pdf) = 1
			for (int k = 0; k < NumSamples; ++k) // phi
			{
				double mag = 0.0;
				double test = 0.0;
				for (int l = 0; l < NumSamples; ++l) // theta
					mag += _pdf[c][i][k][l] * _pdf[c][i][k][l];

				cout << "mag: " << mag << endl;
				for (int l = 0; l < NumSamples; ++l) // theta
				{
					cout << "l: " << l << " pdf: " << _pdf[c][i][k][l] << " | ";
					_pdf[c][i][k][l] /= sqrt(mag);
					test += _pdf[c][i][k][l];
					cout << "l: " << l << " pdf: " << _pdf[c][i][k][l] << endl;
				}

				cout << "mag: " << sqrt(mag) << endl;
				//cout << "test: " << test << " k: " << k << " | ";
				assert(abs(1.0 - test) < 0.1);
				//assert(false);
			}
		}

		// Construct CDFs per channel.
		for (int c = 0; c < 3; ++c) // channel
			for (int i = 0; i < N; ++i) // view theta
				for (int j = 0; j < NumSamples; ++j)
				{
					_cdf[c][i][j][0] = _pdf[c][i][j][0];

					for (int k = 1; k < NumSamples; ++k)
						_cdf[c][i][j][k] = _cdf[c][i][j][k - 1] + _pdf[c][i][j][k];

					double test = _cdf[c][i][j][NumSamples - 1];
					assert(abs(1.0 - test) < 0.00001); // should span from 0 to 1
				}

		free(_raw_brdf_data);
	}

private:
	double *_raw_brdf_data;
	vector<vector<vector<vector<double> > > > _brdf; // channel x view_theta x light_phi x light_theta
	vector<vector<vector<vector<double> > > > _pdf;
	vector<vector<vector<vector<double> > > > _cdf;

	int _NumSamples;
	int _N;
};

#endif
