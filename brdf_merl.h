#ifndef _BRDF_MERL_
#define _BRDF_MERL_

#include <cassert>
#include <cstdlib>
#include <math.h>

#include "glm/glm.hpp"
using namespace glm;

//#include "brdf.h"
#include "MERL.h"

class BrdfMERL
{
public:
	float eval(const vec3& V, const vec3& L, const float alpha, float& pdf) const
	{
		if(V.z <= 0)
		{
			pdf = 0;
			return 0;
		}
	}

	vec3 sample(const vec3& V, const float alpha, const float U1, const float U2) const
	{
		const unsigned int view_theta = static_cast<unsigned int>(round(atanf(V.z / V.x) * _N));
		const unsigned int phi = static_cast<unsigned int>(round(_NumSamples * U1)); // assume phi_light is uniformly distributed

		unsigned int theta = 0;

		for (int i = 0; i < _NumSamples; ++i) // argmin (cdf[i] == rand)
			if (abs(_cdf_r[view_theta][phi][i] - U2) < 0.00001)
			{
				theta = i;
				break;
			}

		int channel = static_cast<int>(alpha);
		switch (channel)
		{
			case 0: theta = 
		}

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

		// allocate N x NumSamples x NumSamples
		_brdf_r.resize(N, vector<vector<double> >(NumSamples, vector<double>(NumSamples, 0)));
		_brdf_g.resize(N, vector<vector<double> >(NumSamples, vector<double>(NumSamples, 0)));
		_brdf_b.resize(N, vector<vector<double> >(NumSamples, vector<double>(NumSamples, 0)));

		_cdf_r.resize(N, vector<vector<double> >(NumSamples, vector<double>(NumSamples, 0)));
		_cdf_g.resize(N, vector<vector<double> >(NumSamples, vector<double>(NumSamples, 0)));
		_cdf_b.resize(N, vector<vector<double> >(NumSamples, vector<double>(NumSamples, 0)));

		// construct PDFs for each channel
		for (int i = 0; i < N; ++i) // view
		{
			double theta_in = i * 0.5 * M_PI / static_cast<double>(N);
			for (int j = 0; j < N; ++j)
			{
				// sum
				for (int k = 0; k < NumSamples; ++k) // light 
				{
					double phi_out = k * 2.0 * M_PI / static_cast<double>(NumSamples);
					for (int l = 0; l < NumSamples; ++l)
					{
						double theta_out = l * 0.5 * M_PI / static_cast<double>(NumSamples);
						double red, green, blue;
						lookup_brdf_val(_raw_brdf_data, theta_in, 0.0, theta_out, phi_out, red, green, blue);

						_brdf_r[i][k][l] += red;
						_brdf_g[i][k][l] += green;
						_brdf_b[i][k][l] += blue;
					}
				}
			}

			// Copy to temp buffer cuz we only want to save the brdf values. not the pdfs
			auto pdf_r = _brdf_r[i];
			auto pdf_g = _brdf_g[i];
			auto pdf_b = _brdf_b[i];

			// average phi contributions
			for (int k = 0; k < NumSamples; ++k)
			for (int l = 0; l < NumSamples; ++l)
			{
				pdf_r[i][k][l] /= static_cast<double>(N);
				pdf_g[i][k][l] /= static_cast<double>(N);
				pdf_b[i][k][l] /= static_cast<double>(N);
			}

			// normalize => integral(pdf) = 1
			double test = 0.0;
			for (int k = 0; k < NumSamples; ++k) //phi
			{
				double mag_r = 0.0, mag_g = 0.0, mag_b = 0.0;
				for (int l = 0; l < NumSamples; ++l) //theta
				{
					mag_r += pdf_r[i][k][l] * pdf_r[i][k][l];
					mag_g += pdf_g[i][k][l] * pdf_g[i][k][l];
					mag_b += pdf_b[i][k][l] * pdf_b[i][k][l];
				}

				for (int l = 0; l < NumSamples; ++l) //theta
				{
					phi_r[i][k][l] /= sqrt(mag_r);
					phi_g[i][k][l] /= sqrt(mag_g);
					phi_b[i][k][l] /= sqrt(mag_b);

					test += phi_r[i][k][l];
				}
			}

			assert(abs(1.0 - test) < 0.00001);
		}

		// construct CDFs
		for (int i = 0; i < N; ++i)
		{
			for (int j = 0; j < NumSamples; ++j)
			{
				_cdf_r[i][j][0] = pdf_r[i][j][0];
				_cdf_g[i][j][0] = pdf_g[i][j][0];
				_cdf_b[i][j][0] = pdf_b[i][j][0];
				for (int k = 1; k < NumSamples; ++k)
				{
					_cdf_r[i][j][k] = _cdf_r[i][j][k - 1] + pdf_r[i][j][k];
					_cdf_g[i][j][k] = _cdf_g[i][j][k - 1] + pdf_g[i][j][k];
					_cdf_b[i][j][k] = _cdf_b[i][j][k - 1] + pdf_b[i][j][k];
				}

				double test = _cdf_r[i][j][NumSamples - 1];

				assert(abs(1.0 - test) < 0.00001); // should span from 0 to 1
			}
		}

		free(_raw_brdf_data);
	}

private:
	double *_raw_brdf_data;
	vector<vector<vector<double> > > _brdf_r;
	vector<vector<vector<double> > > _brdf_g;
	vector<vector<vector<double> > > _brdf_b;

	vector<vector<vector<double> > > _cdf_r;
	vector<vector<vector<double> > > _cdf_g;
	vector<vector<vector<double> > > _cdf_b;

	int _NumSamples;
	int _N
};


#endif
