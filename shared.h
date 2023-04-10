#pragma once

#include "sobol.h"

#include <array>
#include <vector>

typedef std::array<float, 1> float1;
typedef std::array<float, 2> float2;
typedef std::array<float2, 2> float2x2; // indexed by [row][column]

static const float c_pi = 3.14159265359f;
static const float c_goldenRatioConjugate = 0.61803398875f;

template <size_t N>
struct Noise
{
	typedef std::array<float, N> TVec;

	const char* label;
	std::vector<TVec>(*Generate)(pcg32_random_t& rng, int numSamples, std::vector<TVec>& lastSamples) = nullptr;
	std::vector<float> error;
	std::vector<float> avgErrorAtSamples;
	std::vector<float> avgErrorSqAtSamples;
};

inline pcg32_random_t GetRNG(int streamIndex)
{
	pcg32_random_t rng;

#if DETERMINISTIC()
	pcg32_srandom_r(&rng, 0x1337b337, streamIndex);
#else
	std::random_device rd;
	pcg32_srandom_r(&rng, rd(), rd());
#endif

	return rng;
}

inline float RandomFloat01(pcg32_random_t& rng)
{
	return ldexpf((float)pcg32_random_r(&rng), -32);
}

float Lerp(float A, float B, float t)
{
	return A * (1.0f - t) + B * t;
}

inline float RandomFloatRange(pcg32_random_t& rng, float min, float max)
{
	return Lerp(min, max, RandomFloat01(rng));
}

// https://blog.demofox.org/2017/10/01/calculating-the-distance-between-points-in-wrap-around-toroidal-space/
template <size_t N>
float DistanceWrap(const std::array<float, N>& A, const std::array<float, N>& B)
{
	float sumDSquared = 0.0f;
	for (size_t i = 0; i < N; ++i)
	{
		float d = std::abs(A[i] - B[i]);
		d = std::min(d, 1.0f - d);
		sumDSquared += d * d;
	}
	return std::sqrt(sumDSquared);
}

template <size_t N>
float Distance(const std::array<float, N>& A, const std::array<float, N>& B)
{
	float sumDSquared = 0.0f;
	for (size_t i = 0; i < N; ++i)
	{
		float d = A[i] - B[i];
		sumDSquared += d * d;
	}
	return std::sqrt(sumDSquared);
}

template <size_t N>
std::vector<std::array<float, N>> Generate_White(pcg32_random_t& rng, int numSamples, std::vector<std::array<float, N>>& lastSamples)
{
	std::vector<std::array<float, N>> ret(numSamples);
	for (int i = 0; i < numSamples; ++i)
	{
		for (int j = 0; j < N; ++j)
			ret[i][j] = RandomFloat01(rng);
	}
	return ret;
}

template <typename T, size_t N, typename LAMBDA>
void DoTests(const char* label, T(&noiseTypes)[N], int testCount, int pointCount, const char* resultsFileName, const LAMBDA& lambda)
{
	// allocate space for the results of our test.
	// store them all out so we can work multithreadedly, then deterministically combine the results together.
	for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
	{
		T& noise = noiseTypes[noiseIndex];
		noise.error.resize(testCount * pointCount);
		noise.avgErrorAtSamples.resize(pointCount);
		noise.avgErrorSqAtSamples.resize(pointCount);
	}

	// for each test
	std::atomic<int> testsDone = 0;
	int lastPercent = -1;
	#pragma omp parallel for
	for (int testIndex = 0; testIndex < testCount; ++testIndex)
	{
		int threadId = omp_get_thread_num();
		if (threadId == 0)
		{
			int percent = int(100.0f * float(testsDone.load()) / float(testCount));
			if (lastPercent != percent)
			{
				lastPercent = percent;
				printf("\r%s: %i%%", label, percent);
			}
		}

		pcg32_random_t rng = GetRNG(testIndex);

		lambda(testIndex, rng);

		testsDone.fetch_add(1);
	}

	printf("\r%s: 100%%\n", label);

	{
		#pragma omp parallel for
		for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
		{
			T& noise = noiseTypes[noiseIndex];
			for (int testIndex = 0; testIndex < testCount; ++testIndex)
			{
				for (int pointIndex = 0; pointIndex < pointCount; ++pointIndex)
				{
					float error = noise.error[testIndex * pointCount + pointIndex];
					noise.avgErrorAtSamples[pointIndex] = Lerp(noise.avgErrorAtSamples[pointIndex], error, 1.0f / float(testIndex + 1));
					noise.avgErrorSqAtSamples[pointIndex] = Lerp(noise.avgErrorSqAtSamples[pointIndex], error * error, 1.0f / float(testIndex + 1));
				}
			}
		}

		FILE* file = nullptr;
		fopen_s(&file, resultsFileName, "wb");

		// for each sample count
		for (int pointIndex = 0; pointIndex < pointCount; ++pointIndex)
		{
			// write the csv header
			if (pointIndex == 0)
			{
				fprintf(file, "\"samples\"");
				for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
					fprintf(file, ",\"%s\"", noiseTypes[noiseIndex].label);
				fprintf(file, "\n");
			}

			// write the rmse
			fprintf(file, "\"%i\"", pointIndex + 1);
			for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
			{
				T& noise = noiseTypes[noiseIndex];
				float rmse = std::sqrt(noise.avgErrorSqAtSamples[pointIndex] - noise.avgErrorAtSamples[pointIndex] * noise.avgErrorAtSamples[pointIndex]);
				fprintf(file, ",\"%f\"", rmse);
			}
			fprintf(file, "\n");
		}

		fclose(file);
	}
}

inline float Halton(int index, float base)
{
	float result = 0.0f;
	float f = 1.0f / base;
	float i = float(index);

	for (int x = 0; x < 8; x++)
	{
		if (i <= 0.0)
			break;

		result += f * fmod(i, base);
		i = floor(i / base);
		f = f / base;
	}

	return result;
}

inline float2 operator + (const float2& A, const float2& B)
{
	return float2{ A[0] + B[0], A[1] + B[1] };
}

inline float2 operator - (const float2& A, const float2& B)
{
	return float2{ A[0] - B[0], A[1] - B[1] };
}

inline float2 operator * (const float2& A, const float& B)
{
	return float2{ A[0] * B, A[1] * B };
}

inline float2 operator *= (float2& A, const float& B)
{
	A[0] *= B;
	A[1] *= B;
	return A;
}

inline float2 operator /= (float2& A, const float& B)
{
	A[0] /= B;
	A[1] /= B;
	return A;
}

inline float2 operator += (float2& A, const float2& B)
{
	A[0] += B[0];
	A[1] += B[1];
	return A;
}

inline float2 operator -= (float2& A, const float2& B)
{
	A[0] -= B[0];
	A[1] -= B[1];
	return A;
}

inline float2 Normalize(const float2& V)
{
	float len = std::sqrt(V[0] * V[0] + V[1] * V[1]);
	return float2{ V[0] / len, V[1] / len };
}

// Gauss jordan elimination
inline float2x2 Inverse(float2x2 M)
{
	float2x2 ret = { float2{1.0f, 0.0f}, float2{0.0f, 1.0f} };

	// for better results, make sure row 0 col 0 has the largest magnitude of all rows at col 0
	if (std::abs(M[0][0]) < std::abs(M[1][0]))
	{
		std::swap(M[0], M[1]);
		std::swap(ret[0], ret[1]);
	}

	// divide [0][*] by [0][0] to make [0][0] be a 1
	{
		float scale = 1.0f / M[0][0];
		M[0] *= scale;
		ret[0] *= scale;
	}

	// subtract [0][*] from [1][*] multiplied by [1][0], to make it have a 0 in column 0
	{
		float scale = M[1][0];
		M[1] -= M[0] * scale;
		ret[1] -= ret[0] * scale;
	}

	// divide [1][*] by [1][1] to make [1][1] be a 1
	{
		float scale = 1.0f / M[1][1];
		M[1] *= scale;
		ret[1] *= scale;
	}

	// subtract [1][*] from [0][*] multiplied by [0][1], to make it have a 0 in column 1
	{
		float scale = M[0][1];
		M[0] -= M[1] * scale;
		ret[0] -= ret[1] * scale;
	}

	return ret;
}

float2 Transform(const float2& V, const float2x2& M)
{
	return float2{
		V[0] * M[0][0] + V[1] * M[1][0],
		V[0] * M[0][1] + V[1] * M[1][1]
	};
}

float2 Min(const float2& A, const float2& B)
{
	return float2{
		std::min(A[0], B[0]),
		std::min(A[1], B[1])
	};
}

float2 Max(const float2& A, const float2& B)
{
	return float2{
		std::max(A[0], B[0]),
		std::max(A[1], B[1])
	};
}
