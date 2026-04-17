#pragma once

#include <cmath>
#include "blitz/array.h"

using namespace blitz;

enum MAS { // mass assignment scheme
	NGP,
	CIC,
	TSC,
	PCS
};

static inline int wrapIndex(int idx, int n) {
	int w = idx % n;
	return (w < 0) ? (w + n) : w;
}

static inline const char* masName(MAS scheme) {
	switch (scheme) {
		case MAS::NGP: return "NGP";
		case MAS::CIC: return "CIC";
		case MAS::TSC: return "TSC";
		case MAS::PCS: return "PCS";
	}
	return "UNKNOWN";
}

template<MAS S>
struct Stencil1D;

template<>
struct Stencil1D<MAS::NGP> {
	static constexpr int P = 1;

	static inline void eval(float r, int idx[P], float w[P]) {
		idx[0] = static_cast<int>(std::floor(r));
		w[0]   = 1.0f;
	}
};

template<>
struct Stencil1D<MAS::CIC> {
	static constexpr int P = 2;

	static inline void eval(float r, int idx[P], float w[P]) {
		const int i0 = static_cast<int>(std::floor(r - 0.5f));
		const float t = r - (static_cast<float>(i0) + 0.5f); // in [0,1)

		idx[0] = i0;
		idx[1] = i0 + 1;

		w[0] = 1.0f - t;
		w[1] = t;
	}
};

template<>
struct Stencil1D<MAS::TSC> {
	static constexpr int P = 3;

	static inline void eval(float r, int idx[P], float w[P]) {
		const int i0 = static_cast<int>(std::floor(r));
		const float t = r - static_cast<float>(i0); // in [0,1)

		idx[0] = i0 - 1;
		idx[1] = i0;
		idx[2] = i0 + 1;

		w[0] = 0.5f * (1.0f - t) * (1.0f - t);
		w[1] = 0.5f + t - t * t;
		w[2] = 0.5f * t * t;
	}
};

template<>
struct Stencil1D<MAS::PCS> {
	static constexpr int P = 4;

	static inline void eval(float r, int idx[P], float w[P]) {
		const float oneSixth = 1.0f / 6.0f;

		const int i0 = static_cast<int>(std::floor(r - 0.5f));
		const float s = r - (static_cast<float>(i0) + 0.5f); // in [0,1)

		const float s2 = s * s;
		const float s3 = s2 * s;

		const float omt  = 1.0f - s;
		const float omt2 = omt * omt;
		const float omt3 = omt2 * omt;

		idx[0] = i0 - 1;
		idx[1] = i0;
		idx[2] = i0 + 1;
		idx[3] = i0 + 2;

		w[0] = oneSixth * omt3;
		w[1] = oneSixth * (4.0f - 6.0f * s2 + 3.0f * s3);
		w[2] = oneSixth * (1.0f + 3.0f * s + 3.0f * s2 - 3.0f * s3);
		w[3] = oneSixth * s3;
	}
};

template<MAS S>
void assignMassToGridImpl(Array<float,3>& grid,
                          Array<float,2>& r,
                          Array<float,1>& m,
                          int N,
                          int nGrid,
                          float x_min, float x_max,
                          float y_min, float y_max,
                          float z_min, float z_max) {
	grid = 0.0f;

	using Stencil = Stencil1D<S>;
	constexpr int P = Stencil::P;

	#pragma omp parallel for
	for (int pn = r.lbound(0); pn <= r.ubound(0); ++pn) {
		const float rx = (r(pn, 0) - x_min) / (x_max - x_min) * nGrid;
		const float ry = (r(pn, 1) - y_min) / (y_max - y_min) * nGrid;
		const float rz = (r(pn, 2) - z_min) / (z_max - z_min) * nGrid;

		int ix[P], iy[P], iz[P];
		float wx[P], wy[P], wz[P];

		Stencil::eval(rx, ix, wx);
		Stencil::eval(ry, iy, wy);
		Stencil::eval(rz, iz, wz);

		const float mp = m(pn);

		for (int a = 0; a < P; ++a) {
			const int ia = wrapIndex(ix[a], nGrid);
			const float wa = wx[a];

			for (int b = 0; b < P; ++b) {
				const int jb = wrapIndex(iy[b], nGrid);
				const float wab = wa * wy[b];

				for (int c = 0; c < P; ++c) {
					const int kc = wrapIndex(iz[c], nGrid);
					const float W = wab * wz[c];

					#pragma omp atomic
					grid(ia, jb, kc) += mp * W;
				}
			}
		}
	}
}

void assignMassToGrid(Array<float,3>& grid,
                      Array<float,2>& r,
                      Array<float,1>& m,
                      int N,
                      int nGrid,
                      float x_min, float x_max,
                      float y_min, float y_max,
                      float z_min, float z_max,
                      MAS scheme);

void runAllTests();

