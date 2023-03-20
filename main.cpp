#define _CRT_SECURE_NO_WARNINGS

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

#include "pcg/pcg_basic.h"

#include <array>
#include <vector>
#include <string>
#include <random>
#include <direct.h>

#define DETERMINISTIC() false

// TODO: higher point count? and OMP to make it happen

static const int c_1DTestCount = 100;
static const int c_1DTestPointCount = 100;
static const float c_1DTestControlPointMin = -10.0f;
static const float c_1DTestControlPointMax = 10.0f;
static const int c_1DNumPointsReported = 5;

typedef std::array<float, 1> float1;
typedef std::array<float, 2> float2;

float Lerp(float A, float B, float t)
{
	return A * (1.0f - t) + B * t;
}

pcg32_random_t GetRNG()
{
	pcg32_random_t rng;

#if DETERMINISTIC()
	pcg32_srandom_r(&rng, 0x1337b337, 0xbeefcafe);
#else
	std::random_device rd;
	pcg32_srandom_r(&rng, rd(), rd());
#endif

	return rng;
}

float RandomFloat01(pcg32_random_t& rng)
{
	return ldexpf((float)pcg32_random_r(&rng), -32);
}

float RandomFloatRange(pcg32_random_t& rng, float min, float max)
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

// A,B,C,D are control points.
// t is a value between 0 and 1.
float Evaluate1DCubicBezier(float A, float B, float C, float D, float t)
{
	float s = 1.0f - t;
	return
		A * 1.0f * s * s * s +
		B * 3.0f * s * s * t +
		C * 3.0f * s * t * t +
		D * 1.0f * t * t * t;
}

// TODO: can we combine things to make it less computation?
float IndefiniteIntegral1DCubicBezier(float A, float B, float C, float D, float t)
{
	float s = 1.0f - t;
	return
		A * 1.0f * s * s * s * s / -4.0f +
		B * 3.0f * (t * t * t * t / 4.0f - t * t * t * 2.0f / 3.0f + t * t / 2.0f) +
		C * 3.0f * (t * t * t / 3.0f - t * t * t * t / 4.0f) +
		D * 1.0f * t * t * t * t / 4.0f
		;
}

// A,B,C,D are control points.
// returns the definite integral between 0 and 1
float Integral1DCubicBezier(float A, float B, float C, float D)
{
	return IndefiniteIntegral1DCubicBezier(A, B, C, D, 1.0f) - IndefiniteIntegral1DCubicBezier(A, B, C, D, 0.0f);
}

std::vector<float> Generate1D_Regular_Ends(pcg32_random_t& rng, int numSamples)
{
	std::vector<float> ret(numSamples, 0.0f);
	if (numSamples > 1)
	{
		for (int i = 0; i < numSamples; ++i)
			ret[i] = float(i) / float(numSamples - 1);
	}
	return ret;
}

std::vector<float> Generate1D_Regular_Left(pcg32_random_t& rng, int numSamples)
{
	std::vector<float> ret(numSamples);
	for (int i = 0; i < numSamples; ++i)
		ret[i] = float(i) / float(numSamples);
	return ret;
}

std::vector<float> Generate1D_Regular_Center(pcg32_random_t& rng, int numSamples)
{
	std::vector<float> ret(numSamples);
	for (int i = 0; i < numSamples; ++i)
		ret[i] = (float(i) + 0.5f) / float(numSamples);
	return ret;
}

std::vector<float> Generate1D_Regular_Center_Equal(pcg32_random_t& rng, int numSamples)
{
	std::vector<float> ret(numSamples);
	for (int i = 0; i < numSamples; ++i)
		ret[i] = (float(i) + 1.0f) / float(numSamples + 1);
	return ret;
}

std::vector<float> Generate1D_White(pcg32_random_t& rng, int numSamples)
{
	std::vector<float> ret(numSamples);
	for (int i = 0; i < numSamples; ++i)
		ret[i] = RandomFloat01(rng);
	return ret;
}

std::vector<float> Generate1D_Stratified(pcg32_random_t& rng, int numSamples)
{
	std::vector<float> ret = Generate1D_Regular_Left(rng, numSamples);
	for (float& f : ret)
		f += RandomFloatRange(rng, 0.0f, 1.0f / float(numSamples));
	return ret;
}

std::vector<float> Generate1D_Blue_Wrap(pcg32_random_t& rng, int numSamples)
{
	std::vector<float> ret(numSamples);
	for (int sampleIndex = 0; sampleIndex < numSamples; ++sampleIndex)
	{
		float bestCandidateScore = 0.0f;
		float bestCandidate = 0.0f;
		int candidateCount = sampleIndex + 1;
		for (int candidateIndex = 0; candidateIndex < candidateCount; ++candidateIndex)
		{
			float candidate = RandomFloat01(rng);
			float candidateScore = FLT_MAX;
			for (int pointIndex = 0; pointIndex < sampleIndex; ++pointIndex)
				candidateScore = std::min(candidateScore, DistanceWrap(float1{ candidate }, float1{ ret[pointIndex] }));

			if (candidateScore > bestCandidateScore)
			{
				bestCandidate = candidate;
				bestCandidateScore = candidateScore;
			}
		}
		ret[sampleIndex] = bestCandidate;
	}
	return ret;
}

std::vector<float> Generate1D_Blue_NoWrap(pcg32_random_t& rng, int numSamples)
{
	std::vector<float> ret(numSamples);
	for (int sampleIndex = 0; sampleIndex < numSamples; ++sampleIndex)
	{
		float bestCandidateScore = 0.0f;
		float bestCandidate = 0.0f;
		int candidateCount = sampleIndex + 1;
		for (int candidateIndex = 0; candidateIndex < candidateCount; ++candidateIndex)
		{
			float candidate = RandomFloat01(rng);
			float candidateScore = FLT_MAX;
			for (int pointIndex = 0; pointIndex < sampleIndex; ++pointIndex)
				candidateScore = std::min(candidateScore, Distance(float1{ candidate }, float1{ ret[pointIndex] }));

			if (candidateScore > bestCandidateScore)
			{
				bestCandidate = candidate;
				bestCandidateScore = candidateScore;
			}
		}
		ret[sampleIndex] = bestCandidate;
	}
	return ret;
}

std::vector<float> Generate1D_Blue_NoWrap_Edge(pcg32_random_t& rng, int numSamples)
{
	std::vector<float> ret(numSamples);
	for (int sampleIndex = 0; sampleIndex < numSamples; ++sampleIndex)
	{
		float bestCandidateScore = 0.0f;
		float bestCandidate = 0.0f;
		int candidateCount = sampleIndex + 1;
		for (int candidateIndex = 0; candidateIndex < candidateCount; ++candidateIndex)
		{
			float candidate = RandomFloat01(rng);
			float candidateScore = std::min(candidate, 1.0f - candidate); // initialize the score to be the distance to the edge
			for (int pointIndex = 0; pointIndex < sampleIndex; ++pointIndex)
				candidateScore = std::min(candidateScore, Distance(float1{ candidate }, float1{ ret[pointIndex] }));

			if (candidateScore > bestCandidateScore)
			{
				bestCandidate = candidate;
				bestCandidateScore = candidateScore;
			}
		}
		ret[sampleIndex] = bestCandidate;
	}
	return ret;
}

std::vector<float> Generate1D_Blue_NoWrap_HalfEdge(pcg32_random_t& rng, int numSamples)
{
	std::vector<float> ret(numSamples);
	for (int sampleIndex = 0; sampleIndex < numSamples; ++sampleIndex)
	{
		float bestCandidateScore = 0.0f;
		float bestCandidate = 0.0f;
		int candidateCount = sampleIndex + 1;
		for (int candidateIndex = 0; candidateIndex < candidateCount; ++candidateIndex)
		{
			float candidate = RandomFloat01(rng);
			float candidateScore = std::min(candidate, 1.0f - candidate) / 2.0f; // initialize the score to be half the distance to the edge
			for (int pointIndex = 0; pointIndex < sampleIndex; ++pointIndex)
				candidateScore = std::min(candidateScore, Distance(float1{ candidate }, float1{ ret[pointIndex] }));

			if (candidateScore > bestCandidateScore)
			{
				bestCandidate = candidate;
				bestCandidateScore = candidateScore;
			}
		}
		ret[sampleIndex] = bestCandidate;
	}
	return ret;
}

void Do1DTests()
{
	pcg32_random_t rng = GetRNG();

	struct Noise
	{
		const char* label;
		std::vector<float> (*Generate)(pcg32_random_t& rng, int numSamples) = nullptr;
		std::vector<float> avgErrorAtSamples;
		std::vector<float> avgErrorSqAtSamples;
	};

	Noise noiseTypes[] =
	{
		{ "Regular - Ends", Generate1D_Regular_Ends },
		{ "Regular - Left", Generate1D_Regular_Left },
		{ "Regular - Center", Generate1D_Regular_Center },
		{ "Regular - Center Equal", Generate1D_Regular_Center_Equal },
		{ "Stratified", Generate1D_Stratified },
		{ "White", Generate1D_White },
		{ "Blue - Wrap", Generate1D_Blue_Wrap },
		{ "Blue - No Wrap", Generate1D_Blue_NoWrap },
		{ "Blue - No Wrap Edge", Generate1D_Blue_NoWrap_Edge },
		{ "Blue - No Wrap Half Edge", Generate1D_Blue_NoWrap_HalfEdge},
	};

	for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
	{
		Noise& noise = noiseTypes[noiseIndex];
		noise.avgErrorAtSamples.resize(c_1DTestPointCount);
		noise.avgErrorSqAtSamples.resize(c_1DTestPointCount);
	}

	// for each test
	int lastPercent = -1;
	for (int testIndex = 0; testIndex < c_1DTestCount; ++testIndex)
	{
		int percent = int(100.0f * float(testIndex) / float(c_1DTestCount));
		if (lastPercent != percent)
		{
			lastPercent = percent;
			printf("\r%i%%", percent);
		}

		// Generate a random 1d Bezier curve
		float A = RandomFloatRange(rng, c_1DTestControlPointMin, c_1DTestControlPointMax);
		float B = RandomFloatRange(rng, c_1DTestControlPointMin, c_1DTestControlPointMax);
		float C = RandomFloatRange(rng, c_1DTestControlPointMin, c_1DTestControlPointMax);
		float D = RandomFloatRange(rng, c_1DTestControlPointMin, c_1DTestControlPointMax);

		// Calculate the definite integral of the bezier curve, between 0 and 1
		const float c_actualValue = Integral1DCubicBezier(A, B, C, D);

		// for each type of noise
		for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
		{
			Noise& noise = noiseTypes[noiseIndex];

			// for each sample count in the test
			for (int pointIndex = 0; pointIndex < c_1DTestPointCount; ++pointIndex)
			{
				// generate the samples
				std::vector<float> samples = noise.Generate(rng, pointIndex + 1);

				// integrate!
				float yAvg = 0.0f;
				for (int sampleIndex = 0; sampleIndex < pointIndex + 1; ++sampleIndex)
				{
					float y = Evaluate1DCubicBezier(A, B, C, D, samples[sampleIndex]);
					yAvg = Lerp(yAvg, y, 1.0f / float(sampleIndex + 1));
				}

				// track average error and error squared for this sample count, so we can calculate RMSE
				float error = std::abs(yAvg - c_actualValue);
				noise.avgErrorAtSamples[pointIndex] = Lerp(noise.avgErrorAtSamples[pointIndex], error, 1.0f / float(testIndex + 1));
				noise.avgErrorSqAtSamples[pointIndex] = Lerp(noise.avgErrorSqAtSamples[pointIndex], error * error, 1.0f / float(testIndex + 1));
			}
		}
	}
	printf("\r100%%\n");

	// Write out the 1d results CSV
	{
		FILE* file = nullptr;
		fopen_s(&file, "out/1d.csv", "wb");

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
				Noise& noise = noiseTypes[noiseIndex];
				float rmse = std::sqrt(noise.avgErrorSqAtSamples[pointIndex] - noise.avgErrorAtSamples[pointIndex] * noise.avgErrorAtSamples[pointIndex]);
				fprintf(file, ",\"%f\"", rmse);
			}
			fprintf(file, "\n");
		}

		fclose(file);
	}

	// write out example sample points
	{
		// generate the noise types
		std::vector<std::vector<float>> noiseSamplePoints(_countof(noiseTypes));
		for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
			noiseSamplePoints[noiseIndex] = noiseTypes[noiseIndex].Generate(rng, c_1DNumPointsReported);

		// create the file
		FILE* file = nullptr;
		fopen_s(&file, "out/1dpoints.csv", "wb");

		// write the header
		fprintf(file, "\"samples\"");
		for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
			fprintf(file, ",\"%s\"", noiseTypes[noiseIndex].label);
		fprintf(file, "\n");

		// write the sample points
		for (int i = 0; i < c_1DNumPointsReported; ++i)
		{
			fprintf(file, "\"%i\"", i + 1);
			for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
				fprintf(file, ",\"%f\"", noiseSamplePoints[noiseIndex][i]);
			fprintf(file, "\n");
		}

		fclose(file);
	}
}

int main(int argc, char** argv)
{
	_mkdir("out");

	Do1DTests();

	return 0;
}

/*
TODO:

! use omp and boost the test or sample count up!
! do a random non smooth function. random piecewise linear with 4 pieces or something?

----- 1D sampling -----

* Sample Types:
 * 0.0, 1.0 style
 * 0.0, 0.5 style
 * 0.5, 1.0 style
 * 0.25, 0.75 style
 * 0.3, 0.6 style
 * white noise?
 * stratified?
 * golden ratio?
 * MBC? (EBC?) with and without wrap around distance

* Integrate random Bezier functions
* show RMSE at each sample

* make CSV
* draw some low point count numberline for each type of noise, to show it on the blog post.

----- 2D sampling -----

* Sample Types:
 * TODO: fill out

* integrate random bezier functions in a square
* integrate a shape (random functions on it? maybe bezier triangle? or circle? dunno)
 * blue noise with rejection sampling (no edge detection)
 * blue noise with rejection sampling, but also distance to edge being considered in the candidate score
 * same, but half distance to edge
 ! for distance to edge could just test each line segment. could note that an SDF would also work.


Notes:
* this started it: https://mastodon.gamedev.place/@jkaniarz/110032776950329500
* "ok so if you do 2 samples on the unit numberline, you can do it at 0, 0.5.  or 0.5 and 1. or you can center them, which gives you the 1/4, 3/4 setup."

*/