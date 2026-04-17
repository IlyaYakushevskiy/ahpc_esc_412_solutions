#include <vector>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include "omp.h"
#include "assign.h"

/* ===========================
   Test helpers and reference
   =========================== */

namespace {
	constexpr float TEST_EPS = 1e-5f;

	template<MAS S>
	float stencilWeight1D(float r, int i) {
		using Stencil = Stencil1D<S>;
		constexpr int P = Stencil::P;

		int idx[P];
		float w[P];

		Stencil::eval(r, idx, w);

		for (int k = 0; k < P; ++k) {
			if (idx[k] == i)
				return w[k];
		}

		return 0.0f;
	}

	void requireTrue(bool cond, const std::string& msg) {
		if (!cond) {
			throw std::runtime_error(msg);
		}
	}

	void requireNear(float actual, float expected, float eps, const std::string& msg) {
		if (std::fabs(actual - expected) > eps) {
			std::ostringstream os;
			os << msg << " | expected=" << expected << " actual=" << actual
			   << " diff=" << std::fabs(actual - expected);
			throw std::runtime_error(os.str());
		}
	}

	float referenceWeight1D(float r, int i, MAS scheme) {
		const float center = static_cast<float>(i) + 0.5f;
		const float s = std::fabs(center - r);

		switch (scheme) {
			case MAS::NGP:
				return (static_cast<float>(i) <= r && r < static_cast<float>(i + 1)) ? 1.0f : 0.0f;

			case MAS::CIC:
				return (s < 1.0f) ? (1.0f - s) : 0.0f;

			case MAS::TSC:
				return (s < 0.5f) ? (0.75f - s * s)
				     : (s < 1.5f) ? (0.5f * (1.5f - s) * (1.5f - s))
				     : 0.0f;

			case MAS::PCS: {
				const float oneSixth = 1.0f / 6.0f;
				return (s < 1.0f) ? oneSixth * (4.0f - 6.0f * s * s + 3.0f * s * s * s)
				     : (s < 2.0f) ? oneSixth * (2.0f - s) * (2.0f - s) * (2.0f - s)
				     : 0.0f;
			}
		}

		return 0.0f;
	}

	void buildReferenceGrid(Array<float,3>& ref,
	                        Array<float,2>& r,
	                        Array<float,1>& m,
	                        int N,
	                        int nGrid,
	                        float x_min, float x_max,
	                        float y_min, float y_max,
	                        float z_min, float z_max,
	                        MAS scheme) {
		ref = 0.0f;

		for (int pn = 0; pn < N; ++pn) {
			const float rx = (r(pn, 0) - x_min) / (x_max - x_min) * nGrid;
			const float ry = (r(pn, 1) - y_min) / (y_max - y_min) * nGrid;
			const float rz = (r(pn, 2) - z_min) / (z_max - z_min) * nGrid;

			// Wide enough range for every kernel support
			for (int i = static_cast<int>(std::floor(rx)) - 4; i <= static_cast<int>(std::floor(rx)) + 4; ++i) {
				for (int j = static_cast<int>(std::floor(ry)) - 4; j <= static_cast<int>(std::floor(ry)) + 4; ++j) {
					for (int k = static_cast<int>(std::floor(rz)) - 4; k <= static_cast<int>(std::floor(rz)) + 4; ++k) {
						const float W =
							referenceWeight1D(rx, i, scheme) *
							referenceWeight1D(ry, j, scheme) *
							referenceWeight1D(rz, k, scheme);

						if (W != 0.0f) {
							ref(wrapIndex(i, nGrid), wrapIndex(j, nGrid), wrapIndex(k, nGrid)) += m(pn) * W;
						}
					}
				}
			}
		}
	}

	void requireGridNear(Array<float,3>& actual,
	                     Array<float,3>& expected,
	                     int nGrid,
	                     float eps,
	                     const std::string& msg) {
		float maxDiff = 0.0f;
		int maxI = 0, maxJ = 0, maxK = 0;

		for (int i = 0; i < nGrid; ++i) {
			for (int j = 0; j < nGrid; ++j) {
				for (int k = 0; k < nGrid; ++k) {
					const float diff = std::fabs(actual(i,j,k) - expected(i,j,k));
					if (diff > maxDiff) {
						maxDiff = diff;
						maxI = i;
						maxJ = j;
						maxK = k;
					}
				}
			}
		}

		if (maxDiff > eps) {
			std::ostringstream os;
			os << msg
			   << " | maxDiff=" << maxDiff
			   << " at (" << maxI << "," << maxJ << "," << maxK << ")"
			   << " expected=" << expected(maxI,maxJ,maxK)
			   << " actual=" << actual(maxI,maxJ,maxK);
			throw std::runtime_error(os.str());
		}
	}

	float sumGrid(Array<float,3>& grid, int nGrid) {
		float sum = 0.0f;
		for (int i = 0; i < nGrid; ++i) {
			for (int j = 0; j < nGrid; ++j) {
				for (int k = 0; k < nGrid; ++k) {
					sum += grid(i,j,k);
				}
			}
		}
		return sum;
	}

	void testKernelSamples() {
		const std::vector<float> samples = {
			0.0f,
			0.25f,
			0.5f,
			1.0f,
			1.25f,
			2.5f,
			3.999f
		};

		for (int s = NGP; s <= PCS; ++s) {
			const MAS scheme = static_cast<MAS>(s);

			for (float r : samples) {
				for (int i = static_cast<int>(std::floor(r)) - 4;
					i <= static_cast<int>(std::floor(r)) + 4;
					++i) {

					float w;

					switch (scheme) {
						case MAS::NGP: w = stencilWeight1D<MAS::NGP>(r,i); break;
						case MAS::CIC: w = stencilWeight1D<MAS::CIC>(r,i); break;
						case MAS::TSC: w = stencilWeight1D<MAS::TSC>(r,i); break;
						case MAS::PCS: w = stencilWeight1D<MAS::PCS>(r,i); break;
					}

					std::ostringstream os;
					os << "Kernel sample mismatch for "
					<< masName(scheme)
					<< " r=" << r
					<< " i=" << i;

					requireNear(w, referenceWeight1D(r,i,scheme), TEST_EPS, os.str());
				}
			}
		}
	}

	void testPartitionOfUnityAndNonNegativity() {
		for (int s = NGP; s <= PCS; ++s) {
			const MAS scheme = static_cast<MAS>(s);

			for (int step = 0; step <= 400; ++step) {
				const float r = 4.0f * static_cast<float>(step) / 400.0f;

				float sum = 0.0f;

				for (int i = static_cast<int>(std::floor(r)) - 4;
					i <= static_cast<int>(std::floor(r)) + 4;
					++i) {

					float w = 0.0f;

					switch (scheme) {
						case MAS::NGP: w = stencilWeight1D<MAS::NGP>(r,i); break;
						case MAS::CIC: w = stencilWeight1D<MAS::CIC>(r,i); break;
						case MAS::TSC: w = stencilWeight1D<MAS::TSC>(r,i); break;
						case MAS::PCS: w = stencilWeight1D<MAS::PCS>(r,i); break;
					}

					std::ostringstream nonNegMsg;
					nonNegMsg << "Negative weight for "
							<< masName(scheme)
							<< " at r=" << r
							<< " i=" << i;

					requireTrue(w >= -TEST_EPS, nonNegMsg.str());

					sum += w;
				}

				std::ostringstream sumMsg;
				sumMsg << "Partition-of-unity failed for "
					<< masName(scheme)
					<< " at r=" << r;

				requireNear(sum, 1.0f, 5e-5f, sumMsg.str());
			}
		}
	}

	void testSingleParticleStencilAgainstReference() {
		const int nGrid = 8;
		Array<float,2> r(1,3);
		Array<float,1> m(1);
		Array<float,3> actual(nGrid,nGrid,nGrid);
		Array<float,3> expected(nGrid,nGrid,nGrid);

		r(0,0) = 3.25f;
		r(0,1) = 4.50f;
		r(0,2) = 5.75f;
		m(0)   = 2.5f;

		for (int s = NGP; s <= PCS; ++s) {
			const MAS scheme = static_cast<MAS>(s);

			assignMassToGrid(actual, r, m, 1, nGrid,
			                 0.0f, static_cast<float>(nGrid),
			                 0.0f, static_cast<float>(nGrid),
			                 0.0f, static_cast<float>(nGrid),
			                 scheme);

			buildReferenceGrid(expected, r, m, 1, nGrid,
			                   0.0f, static_cast<float>(nGrid),
			                   0.0f, static_cast<float>(nGrid),
			                   0.0f, static_cast<float>(nGrid),
			                   scheme);

			std::ostringstream os;
			os << "Single-particle stencil mismatch for " << masName(scheme);
			requireGridNear(actual, expected, nGrid, 5e-5f, os.str());
		}
	}

	void testPeriodicWrappingAgainstReference() {
		const int nGrid = 4;
		Array<float,2> r(1,3);
		Array<float,1> m(1);
		Array<float,3> actual(nGrid,nGrid,nGrid);
		Array<float,3> expected(nGrid,nGrid,nGrid);

		r(0,0) = 3.9f;  // near upper x-boundary, should wrap for CIC/TSC/PCS
		r(0,1) = 0.5f;
		r(0,2) = 0.5f;
		m(0)   = 1.0f;

		for (int s = NGP; s <= PCS; ++s) {
			const MAS scheme = static_cast<MAS>(s);

			assignMassToGrid(actual, r, m, 1, nGrid,
			                 0.0f, static_cast<float>(nGrid),
			                 0.0f, static_cast<float>(nGrid),
			                 0.0f, static_cast<float>(nGrid),
			                 scheme);

			buildReferenceGrid(expected, r, m, 1, nGrid,
			                   0.0f, static_cast<float>(nGrid),
			                   0.0f, static_cast<float>(nGrid),
			                   0.0f, static_cast<float>(nGrid),
			                   scheme);

			std::ostringstream os;
			os << "Periodic wrapping mismatch for " << masName(scheme);
			requireGridNear(actual, expected, nGrid, 5e-5f, os.str());
		}
	}

	void testMassConservation() {
		const int nGrid = 8;
		const int N = 4;

		Array<float,2> r(N,3);
		Array<float,1> m(N);
		Array<float,3> grid(nGrid,nGrid,nGrid);

		// Includes exact lower-boundary point to catch the old NGP bug.
		r(0,0) = 0.0f;   r(0,1) = 0.5f;   r(0,2) = 0.5f;   m(0) = 1.0f;
		r(1,0) = 1.0f;   r(1,1) = 1.0f;   r(1,2) = 1.0f;   m(1) = 2.0f;
		r(2,0) = 3.25f;  r(2,1) = 4.75f;  r(2,2) = 6.125f; m(2) = 3.5f;
		r(3,0) = 7.9f;   r(3,1) = 0.1f;   r(3,2) = 2.2f;   m(3) = 0.5f;

		const float totalMass = 1.0f + 2.0f + 3.5f + 0.5f;

		for (int s = NGP; s <= PCS; ++s) {
			const MAS scheme = static_cast<MAS>(s);

			assignMassToGrid(grid, r, m, N, nGrid,
			                 0.0f, static_cast<float>(nGrid),
			                 0.0f, static_cast<float>(nGrid),
			                 0.0f, static_cast<float>(nGrid),
			                 scheme);

			std::ostringstream os;
			os << "Mass conservation failed for " << masName(scheme);
			requireNear(sumGrid(grid, nGrid), totalMass, 1e-4f, os.str());
		}
	}
}

void runAllTests() {
	const int originalThreads = omp_get_max_threads();
	omp_set_num_threads(1);

	int passed = 0;
	int failed = 0;

	auto run = [&](const std::string& name, void (*fn)()) {
		try {
			fn();
			std::cout << "[PASS] " << name << std::endl;
			++passed;
		} catch (const std::exception& e) {
			std::cerr << "[FAIL] " << name << " -> " << e.what() << std::endl;
			++failed;
		}
	};

	run("kernel samples", testKernelSamples);
	run("partition of unity + non-negativity", testPartitionOfUnityAndNonNegativity);
	run("single-particle stencil vs reference", testSingleParticleStencilAgainstReference);
	run("periodic wrapping vs reference", testPeriodicWrappingAgainstReference);
	run("mass conservation", testMassConservation);

	omp_set_num_threads(originalThreads);

	if (failed != 0) {
		std::ostringstream os;
		os << failed << " test(s) failed";
		throw std::runtime_error(os.str());
	}

	std::cout << "All " << passed << " tests passed." << std::endl;
}