#pragma once

#include <array>
#include <vector>

typedef std::array<float, 1> float1;
typedef std::array<float, 2> float2;

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
