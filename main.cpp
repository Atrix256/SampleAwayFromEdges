#define _CRT_SECURE_NO_WARNINGS

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

#include "pcg/pcg_basic.h"

#include <array>
#include <vector>
#include <string>
#include <random>
#include <direct.h>
#include <algorithm>
#include <omp.h>
#include <atomic>

#define DETERMINISTIC() false

static const int	c_1DTestCount = 1000;
static const int	c_1DTestPointCount = 200;
static const float	c_1DTestControlPointMin = -10.0f;
static const float	c_1DTestControlPointMax = 10.0f;
static const int	c_1DNumPointsReported = 20;

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

// Not the most efficient, just did an integration of each term
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

template <size_t N>
float EvaluateLinearPiecewise(const float(&pieces)[N], float t)
{
	int numPieces = N / 2;
	float pieceIndexF = t * float(numPieces);
	int pieceIndex = std::min(int(pieceIndexF), numPieces - 1);
	float pieceT = pieceIndexF - float(pieceIndex);
	return Lerp(pieces[pieceIndex * 2], pieces[pieceIndex * 2 + 1], pieceT);
}

template <size_t N>
float IntegralLinearPiecewise(const float(&pieces)[N])
{
	int numPieces = N / 2;
	float ret = 0.0f;
	for (int i = 0; i < numPieces; ++i)
	{
		float averageY = (pieces[i * 2] + pieces[i * 2 + 1]) / 2.0f; // get the average height of the piece
		ret += averageY / float(numPieces); // scale by the width of the piece
	}
	return ret;
}

std::vector<float> Generate1D_Regular_Ends(pcg32_random_t& rng, int numSamples, std::vector<float>& lastSamples)
{
	std::vector<float> ret(numSamples, 0.0f);
	if (numSamples > 1)
	{
		for (int i = 0; i < numSamples; ++i)
			ret[i] = float(i) / float(numSamples - 1);
	}
	return ret;
}

std::vector<float> Generate1D_Regular_Left(pcg32_random_t& rng, int numSamples, std::vector<float>& lastSamples)
{
	std::vector<float> ret(numSamples);
	for (int i = 0; i < numSamples; ++i)
		ret[i] = float(i) / float(numSamples);
	return ret;
}

std::vector<float> Generate1D_Regular_Center(pcg32_random_t& rng, int numSamples, std::vector<float>& lastSamples)
{
	std::vector<float> ret(numSamples);
	for (int i = 0; i < numSamples; ++i)
		ret[i] = (float(i) + 0.5f) / float(numSamples);
	return ret;
}

std::vector<float> Generate1D_Regular_Center_Equal(pcg32_random_t& rng, int numSamples, std::vector<float>& lastSamples)
{
	std::vector<float> ret(numSamples);
	for (int i = 0; i < numSamples; ++i)
		ret[i] = (float(i) + 1.0f) / float(numSamples + 1);
	return ret;
}

std::vector<float> Generate1D_White(pcg32_random_t& rng, int numSamples, std::vector<float>& lastSamples)
{
	std::vector<float> ret(numSamples);
	for (int i = 0; i < numSamples; ++i)
		ret[i] = RandomFloat01(rng);
	return ret;
}

std::vector<float> Generate1D_GoldenRatio(pcg32_random_t& rng, int numSamples, std::vector<float>& lastSamples)
{
	static const float c_goldenRatioConjugate = 0.61803398875f;
	float lastValue = 0.0f;
	std::vector<float> ret(numSamples);
	for (float& f : ret)
	{
		lastValue = std::fmodf(lastValue + c_goldenRatioConjugate, 1.0f);
		f = lastValue;
	}
	return ret;
}


std::vector<float> Generate1D_Stratified(pcg32_random_t& rng, int numSamples, std::vector<float>& lastSamples)
{
	std::vector<float> ret = Generate1D_Regular_Left(rng, numSamples, lastSamples);
	for (float& f : ret)
		f += RandomFloatRange(rng, 0.0f, 1.0f / float(numSamples));
	return ret;
}

std::vector<float> Generate1D_Blue_Wrap(pcg32_random_t& rng, int numSamples, std::vector<float>& lastSamples)
{
	std::vector<float> ret;
	int startSampleIndex = 0;
	if (lastSamples.size() < numSamples)
	{
		startSampleIndex = (int)lastSamples.size();
		ret = lastSamples;
	}
	ret.resize(numSamples);
	for (int sampleIndex = startSampleIndex; sampleIndex < numSamples; ++sampleIndex)
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

std::vector<float> Generate1D_Blue_NoWrap(pcg32_random_t& rng, int numSamples, std::vector<float>& lastSamples)
{
	std::vector<float> ret;
	int startSampleIndex = 0;
	if (lastSamples.size() < numSamples)
	{
		startSampleIndex = (int)lastSamples.size();
		ret = lastSamples;
	}
	ret.resize(numSamples);
	for (int sampleIndex = startSampleIndex; sampleIndex < numSamples; ++sampleIndex)
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

std::vector<float> Generate1D_Blue_NoWrap_Edge(pcg32_random_t& rng, int numSamples, std::vector<float>& lastSamples)
{
	std::vector<float> ret;
	int startSampleIndex = 0;
	if (lastSamples.size() < numSamples)
	{
		startSampleIndex = (int)lastSamples.size();
		ret = lastSamples;
	}
	ret.resize(numSamples);
	for (int sampleIndex = startSampleIndex; sampleIndex < numSamples; ++sampleIndex)
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

std::vector<float> Generate1D_Blue_NoWrap_HalfEdge(pcg32_random_t& rng, int numSamples, std::vector<float>& lastSamples)
{
	std::vector<float> ret;
	int startSampleIndex = 0;
	if (lastSamples.size() < numSamples)
	{
		startSampleIndex = (int)lastSamples.size();
		ret = lastSamples;
	}
	ret.resize(numSamples);
	for (int sampleIndex = startSampleIndex; sampleIndex < numSamples; ++sampleIndex)
	{
		float bestCandidateScore = 0.0f;
		float bestCandidate = 0.0f;
		int candidateCount = sampleIndex + 1;
		for (int candidateIndex = 0; candidateIndex < candidateCount; ++candidateIndex)
		{
			float candidate = RandomFloat01(rng);
			float candidateScore = std::min(candidate, 1.0f - candidate) * 2.0f; // initialize the score to be twice the distance to the edge (make distance to edge count half as much)
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
	printf("==================== 1D ====================\n");

	pcg32_random_t rng = GetRNG();

	struct Noise
	{
		const char* label;
		std::vector<float>(*Generate)(pcg32_random_t& rng, int numSamples, std::vector<float>& lastSamples) = nullptr;
		std::vector<float> error;
		std::vector<float> avgErrorAtSamples;
		std::vector<float> avgErrorSqAtSamples;
	};

	Noise noiseTypes[] =
	{
		{ "Regular - Ends", Generate1D_Regular_Ends },
		{ "Regular - Left", Generate1D_Regular_Left },
		{ "Regular - Center", Generate1D_Regular_Center },
		{ "Regular - Center Equal", Generate1D_Regular_Center_Equal },
		{ "Golden Ratio", Generate1D_GoldenRatio },
		{ "Stratified", Generate1D_Stratified },
		{ "White", Generate1D_White },
		{ "Blue - Wrap", Generate1D_Blue_Wrap },
		{ "Blue - No Wrap", Generate1D_Blue_NoWrap },
		{ "Blue - No Wrap Edge", Generate1D_Blue_NoWrap_Edge },
		{ "Blue - No Wrap Half Edge", Generate1D_Blue_NoWrap_HalfEdge},
	};

	// allocate space for the results of our test.
	// store them all out so we can work multithreadedly, then deterministically combine the results together.
	for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
	{
		Noise& noise = noiseTypes[noiseIndex];
		noise.error.resize(c_1DTestCount * c_1DTestPointCount);
		noise.avgErrorAtSamples.resize(c_1DTestPointCount);
		noise.avgErrorSqAtSamples.resize(c_1DTestPointCount);
	}

	// smooth tests
	{
		printf("Smooth Functions...\n");

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
					printf("\r%i%%", percent);
				}
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
				std::vector<float> samples;
				for (int pointIndex = 0; pointIndex < c_1DTestPointCount; ++pointIndex)
				{
					// generate the samples
					samples = noise.Generate(rng, pointIndex + 1, samples);

					// integrate!
					float yAvg = 0.0f;
					for (int sampleIndex = 0; sampleIndex < pointIndex + 1; ++sampleIndex)
					{
						float y = Evaluate1DCubicBezier(A, B, C, D, samples[sampleIndex]);
						yAvg = Lerp(yAvg, y, 1.0f / float(sampleIndex + 1));
					}

					// store the error
					noise.error[testIndex * c_1DTestPointCount + pointIndex] = std::abs(yAvg - c_actualValue);
				}
			}
			testsDone.fetch_add(1);
		}
		printf("\r100%%\n");

		// Gather the results
		{
			#pragma omp parallel for
			for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
			{
				Noise& noise = noiseTypes[noiseIndex];
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
		}

		// Write out the 1d results CSV
		{
			FILE* file = nullptr;
			fopen_s(&file, "out/1DResultsSmooth.csv", "wb");

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
	}

	// non smooth tests
	{
		printf("Non Smooth Functions...\n");

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
					printf("\r%i%%", percent);
				}
			}

			// Generate a random 1d Bezier curve
			float points[8] = {
				RandomFloatRange(rng, c_1DTestControlPointMin, c_1DTestControlPointMax),
				RandomFloatRange(rng, c_1DTestControlPointMin, c_1DTestControlPointMax),
				RandomFloatRange(rng, c_1DTestControlPointMin, c_1DTestControlPointMax),
				RandomFloatRange(rng, c_1DTestControlPointMin, c_1DTestControlPointMax),
				RandomFloatRange(rng, c_1DTestControlPointMin, c_1DTestControlPointMax),
				RandomFloatRange(rng, c_1DTestControlPointMin, c_1DTestControlPointMax),
				RandomFloatRange(rng, c_1DTestControlPointMin, c_1DTestControlPointMax),
				RandomFloatRange(rng, c_1DTestControlPointMin, c_1DTestControlPointMax)
			};

			// Calculate the definite integral of the bezier curve, between 0 and 1
			const float c_actualValue = IntegralLinearPiecewise(points);

			// for each type of noise
			for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
			{
				Noise& noise = noiseTypes[noiseIndex];

				// for each sample count in the test
				std::vector<float> samples;
				for (int pointIndex = 0; pointIndex < c_1DTestPointCount; ++pointIndex)
				{
					// generate the samples
					samples = noise.Generate(rng, pointIndex + 1, samples);

					// integrate!
					float yAvg = 0.0f;
					for (int sampleIndex = 0; sampleIndex < pointIndex + 1; ++sampleIndex)
					{
						float y = EvaluateLinearPiecewise(points, samples[sampleIndex]);
						yAvg = Lerp(yAvg, y, 1.0f / float(sampleIndex + 1));
					}

					// store the error
					noise.error[testIndex * c_1DTestPointCount + pointIndex] = std::abs(yAvg - c_actualValue);
				}
			}
			testsDone.fetch_add(1);
		}
		printf("\r100%%\n");

		// Gather the results
		{
			#pragma omp parallel for
			for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
			{
				Noise& noise = noiseTypes[noiseIndex];
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
		}

		// Write out the 1d results CSV
		{
			FILE* file = nullptr;
			fopen_s(&file, "out/1DResultsNonSmooth.csv", "wb");

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
	}

	// write out example sample points
	{
		// generate the noise types
		std::vector<std::vector<float>> noiseSamplePoints(_countof(noiseTypes));
		for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
		{
			std::vector<float> samples;
			noiseSamplePoints[noiseIndex] = noiseTypes[noiseIndex].Generate(rng, c_1DNumPointsReported, samples);
		}

		// create the file
		FILE* file = nullptr;
		fopen_s(&file, "out/1DPoints.csv", "wb");

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