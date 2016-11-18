#ifndef _EXPORT_
#define _EXPORT_

// export data in C
void writeTabC(mat3 * tab, vec2 * tabAmplitude, int N)
{
	ofstream file("results/ltc.inc");

	file << std::fixed;
	file << std::setprecision(6);

	file << "static const int size = " << N  << ";" << endl << endl;

	file << "static const mat33 tabM[size*size] = {" << endl;
	for(int t = 0 ; t < N ; ++t)
	for(int a = 0 ; a < N ; ++a)
	{
		file << "{";
		file << tab[a + t*N][0][0] << ", " << tab[a + t*N][0][1] << ", " << tab[a + t*N][0][2] << ", ";
		file << tab[a + t*N][1][0] << ", " << tab[a + t*N][1][1] << ", " << tab[a + t*N][1][2] << ", ";
        file << tab[a + t*N][2][0] << ", " << tab[a + t*N][2][1] << ", " << tab[a + t*N][2][2] << "}";
		if(a != N-1 || t != N-1)
			file << ", ";
		file << endl;
	}
	file << "};" << endl << endl;

	file << "static const mat33 tabMinv[size*size] = {" << endl;
	for(int t = 0 ; t < N ; ++t)
	for(int a = 0 ; a < N ; ++a)
	{
		mat3 Minv = glm::inverse(tab[a + t*N]);

		file << "{";
		file << Minv[0][0] << ", " << Minv[0][1] << ", " << Minv[0][2] << ", ";
		file << Minv[1][0] << ", " << Minv[1][1] << ", " << Minv[1][2] << ", ";
        file << Minv[2][0] << ", " << Minv[2][1] << ", " << Minv[2][2] << "}";
		if(a != N-1 || t != N-1)
			file << ", ";
		file << endl;
	}
	file << "};" << endl << endl;

	file << "static const float tabAmplitude[size*size] = {" << endl;
	for(int t = 0 ; t < N ; ++t)
	for(int a = 0 ; a < N ; ++a)
	{
		file << tabAmplitude[a + t*N][0] << "f";
		if(a != N-1 || t != N-1)
			file << ", ";
		file << endl;
	}
	file << "};" << endl;

	file.close();
}

// export data in matlab
void writeTabMatlab(mat3 * tab, vec2 * tabAmplitude, int N)
{
	ofstream file("results/ltc.mat");

	file << "# name: tabAmplitude" << endl;
	file << "# type: matrix" << endl;
	file << "# ndims: 2" << endl;
	file << " " << N << " " << N << endl;

	for(int t = 0 ; t < N ; ++t)
	{
		for(int a = 0 ; a < N ; ++a)
		{
			file << tabAmplitude[a + t*N][0] << " " ;
		}
		file << endl;
	}

	for(int row = 0 ; row<3 ; ++row)
	for(int column = 0 ; column<3 ; ++column)
	{

		file << "# name: tab" << column << row << endl;
		file << "# type: matrix" << endl;
		file << "# ndims: 2" << endl;
		file << " " << N << " " << N << endl;

		for(int t = 0 ; t < N ; ++t)
		{
			for(int a = 0 ; a < N ; ++a)
			{
				file << tab[a + t*N][column][row] << " " ;
			}
			file << endl;
		}

		file << endl;
	}

	file.close();
}

#include <iostream>

// export data in JSON
// only prints out M^-1 and mag
void writeTabJSON(mat3 * tab, vec2 * tabAmplitude, int N)
{
	cout << endl << endl;

	ofstream file("results/ltc.js");

	file << std::fixed;
	file << std::setprecision(6);
	file << "var g_ltc_Minv = [" << endl;
	for(int t = 0 ; t < N ; ++t)
	for(int a = 0 ; a < N ; ++a)
	{
		float a0 = tab[a + t*N][0][0];
		float b0 = tab[a + t*N][0][2];
		float c0 = tab[a + t*N][1][1];
		float d0 = tab[a + t*N][2][0];

		float a1 = a0;
		float b1 = -b0;
		float c1 = (a0 - b0*d0) / c0;
		float d1 = -d0;

		cout << "a: " << a0 << " b: " << b0 << " c: " << c0 << " d: " << d0 << endl << endl;

		file << a1 << ", " << b1 << ", " << c1 << ", " << d1;
		if(a != N-1 || t != N-1)
			file << ", ";
	}
	file << "];" << endl;

	file << "var g_ltc_Mag = [" << endl;
	for(int t = 0 ; t < N ; ++t)
	for(int a = 0 ; a < N ; ++a)
	{
		file << tabAmplitude[a + t*N][0];
		if(a != N-1 || t != N-1)
			file << ", ";
	}
	file << "];" << endl;

	file.close();
}

// export data in dds
#include "dds.h"

void writeDDS(mat3 * tab, vec2 * tabAmplitude, int N)
{
	float * data = new float[N*N*4];

	int n = 0;
	for (int i = 0; i < N*N; ++i, n += 4)
	{
		const mat3& m = tab[i];

		float a = m[0][0];
		float b = m[0][2];
		float c = m[1][1];
		float d = m[2][0];

		// Rescaled inverse of m:
		// a 0 b   inverse   1      0      -b
		// 0 c 0     ==>     0 (a - b*d)/c  0
		// d 0 1            -d      0       a

		// Store the variable terms
		data[n + 0] =  a;
		data[n + 1] = -b;
		data[n + 2] = (a - b*d) / c;
		data[n + 3] = -d;
	}

	SaveDDS("results/ltc_mat.dds", DDS_FORMAT_R32G32B32A32_FLOAT, sizeof(float)*4, N, N, data);
	SaveDDS("results/ltc_amp.dds", DDS_FORMAT_R32G32_FLOAT,       sizeof(float)*2, N, N, tabAmplitude);

	delete [] data;
}

#include "Cimg.h"
using namespace cimg_library;

void writePNG(mat3 * tab, vec2 * tabAmplitude, int N)
{
	CImg<float> image_M(N, N, 1, 4);
	CImg<float> image_Minv(N, N, 1, 4);
	CImg<float> image_Amp(N, N, 1, 1);

	for (int i = 0; i < N; ++i)
	for (int j = 0; j < N; ++j)
	{
		const mat3& m = tab[j + i*N];
		
		float a = m[0][0];
		float b = m[0][2];
		float c = m[1][1];
		float d = m[2][0];

		image_M(i, j, 0, 0) = a;
		image_M(i, j, 0, 1) = b;
		image_M(i, j, 0, 2) = c;
		image_M(i, j, 0, 3) = d;
	
		image_Minv(i, j, 0, 0) = a;
		image_Minv(i, j, 0, 1) = -b;
		image_Minv(i, j, 0, 2) = (a - b*d) / c;
		image_Minv(i, j, 0, 3) = -d;

		image_Amp(i, j, 0, 0) = tabAmplitude[j + i*N][1];
		
		std::cout << "num: " << j + i*N << "\n";
		std::cout << a << " " << b << " " << c << " " << d << "\n\n";
	}

	image_M.save("ltc_M.png");	
	image_Minv.save("ltc_Minv.png");	
	image_Amp.save("ltc_Amp.png");	
}

#endif
