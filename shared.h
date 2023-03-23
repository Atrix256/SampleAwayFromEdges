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

template <typename LAMBDA>
void DoTests(const char* label, int testCount, const LAMBDA& lambda)
{

	// for each test
	std::atomic<int> testsDone = 0;
	int lastPercent = -1;
	#pragma omp parallel for
	for (int testIndex = 0; testIndex < c_1DTestCount; ++testIndex)
	{
		int threadId = omp_get_thread_num();
		if (threadId == 0)
		{
			int percent = int(100.0f * float(testsDone.load()) / float(c_1DTestCount));
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
}

template <typename T, size_t N>
void WriteNoiseResults(T(&noiseTypes)[N], const char* fileName)
{
	#pragma omp parallel for
	for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
	{
		T& noise = noiseTypes[noiseIndex];
		for (int testIndex = 0; testIndex < c_1DTestCount; ++testIndex)
		{
			for (int pointIndex = 0; pointIndex < c_1DTestPointCount; ++pointIndex)
			{
				float error = noise.error[testIndex * c_1DTestPointCount + pointIndex];
				noise.avgErrorAtSamples[pointIndex] = Lerp(noise.avgErrorAtSamples[pointIndex], error, 1.0f / float(testIndex + 1));
				noise.avgErrorSqAtSamples[pointIndex] = Lerp(noise.avgErrorSqAtSamples[pointIndex], error * error, 1.0f / float(testIndex + 1));
			}
		}
	}

	FILE* file = nullptr;
	fopen_s(&file, fileName, "wb");

	// for each sample count
	for (int pointIndex = 0; pointIndex < c_1DTestPointCount; ++pointIndex)
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
			Noise<1>& noise = noiseTypes[noiseIndex];
			float rmse = std::sqrt(noise.avgErrorSqAtSamples[pointIndex] - noise.avgErrorAtSamples[pointIndex] * noise.avgErrorAtSamples[pointIndex]);
			fprintf(file, ",\"%f\"", rmse);
		}
		fprintf(file, "\n");
	}

	fclose(file);
}
