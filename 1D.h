#pragma once

#include "shared.h"

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

std::vector<float1> Generate1D_Regular_Ends(pcg32_random_t& rng, int numSamples, std::vector<float1>& lastSamples)
{
	std::vector<float1> ret(numSamples);
	if (numSamples > 1)
	{
		for (int i = 0; i < numSamples; ++i)
			ret[i][0] = float(i) / float(numSamples - 1);
	}
	return ret;
}

std::vector<float1> Generate1D_Regular_Left(pcg32_random_t& rng, int numSamples, std::vector<float1>& lastSamples)
{
	std::vector<float1> ret(numSamples);
	for (int i = 0; i < numSamples; ++i)
		ret[i][0] = float(i) / float(numSamples);
	return ret;
}

std::vector<float1> Generate1D_Regular_Center(pcg32_random_t& rng, int numSamples, std::vector<float1>& lastSamples)
{
	std::vector<float1> ret(numSamples);
	for (int i = 0; i < numSamples; ++i)
		ret[i][0] = (float(i) + 0.5f) / float(numSamples);
	return ret;
}

std::vector<float1> Generate1D_Regular_Center_Equal(pcg32_random_t& rng, int numSamples, std::vector<float1>& lastSamples)
{
	std::vector<float1> ret(numSamples);
	for (int i = 0; i < numSamples; ++i)
		ret[i][0] = (float(i) + 1.0f) / float(numSamples + 1);
	return ret;
}

std::vector<float1> Generate1D_GoldenRatio(pcg32_random_t& rng, int numSamples, std::vector<float1>& lastSamples)
{
	static const float c_goldenRatioConjugate = 0.61803398875f;
	float lastValue = 0.0f;
	std::vector<float1> ret(numSamples);
	for (float1& f : ret)
	{
		lastValue = std::fmodf(lastValue + c_goldenRatioConjugate, 1.0f);
		f[0] = lastValue;
	}
	return ret;
}


std::vector<float1> Generate1D_Stratified(pcg32_random_t& rng, int numSamples, std::vector<float1>& lastSamples)
{
	std::vector<float1> ret = Generate1D_Regular_Left(rng, numSamples, lastSamples);
	for (float1& f : ret)
		f[0] += RandomFloatRange(rng, 0.0f, 1.0f / float(numSamples));
	return ret;
}

std::vector<float1> Generate1D_Blue_Wrap(pcg32_random_t& rng, int numSamples, std::vector<float1>& lastSamples)
{
	std::vector<float1> ret;
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
		ret[sampleIndex][0] = bestCandidate;
	}
	return ret;
}

std::vector<float1> Generate1D_Blue_NoWrap(pcg32_random_t& rng, int numSamples, std::vector<float1>& lastSamples)
{
	std::vector<float1> ret;
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
		ret[sampleIndex][0] = bestCandidate;
	}
	return ret;
}

std::vector<float1> Generate1D_Blue_NoWrap_Edge(pcg32_random_t& rng, int numSamples, std::vector<float1>& lastSamples)
{
	std::vector<float1> ret;
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
		ret[sampleIndex][0] = bestCandidate;
	}
	return ret;
}

std::vector<float1> Generate1D_Blue_NoWrap_HalfEdge(pcg32_random_t& rng, int numSamples, std::vector<float1>& lastSamples)
{
	std::vector<float1> ret;
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
		ret[sampleIndex][0] = bestCandidate;
	}
	return ret;
}

void Do1DTests()
{
	printf("==================== 1D ====================\n");

	Noise<1> noiseTypes[] =
	{
		{ "Regular - Ends", Generate1D_Regular_Ends },
		{ "Regular - Left", Generate1D_Regular_Left },
		{ "Regular - Center", Generate1D_Regular_Center },
		{ "Regular - Center Equal", Generate1D_Regular_Center_Equal },
		{ "Golden Ratio", Generate1D_GoldenRatio },
		{ "Stratified", Generate1D_Stratified },
		{ "White", Generate_White<1> },
		{ "Blue - Wrap", Generate1D_Blue_Wrap },
		{ "Blue - No Wrap", Generate1D_Blue_NoWrap },
		{ "Blue - No Wrap Edge", Generate1D_Blue_NoWrap_Edge },
		{ "Blue - No Wrap Half Edge", Generate1D_Blue_NoWrap_HalfEdge},
	};

	// smooth tests
	DoTests(
		"Smooth",
		noiseTypes,
		c_1DTestCount,
		c_1DTestPointCount,
		"out/1DResultsSmooth.csv",
		[&] (int testIndex, pcg32_random_t& rng)
		{
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
				Noise<1>& noise = noiseTypes[noiseIndex];

				// for each sample count in the test
				std::vector<float1> samples;
				for (int pointIndex = 0; pointIndex < c_1DTestPointCount; ++pointIndex)
				{
					// generate the samples
					samples = noise.Generate(rng, pointIndex + 1, samples);

					// integrate!
					float yAvg = 0.0f;
					for (int sampleIndex = 0; sampleIndex < pointIndex + 1; ++sampleIndex)
					{
						float y = Evaluate1DCubicBezier(A, B, C, D, samples[sampleIndex][0]);
						yAvg = Lerp(yAvg, y, 1.0f / float(sampleIndex + 1));
					}

					// store the error
					noise.error[testIndex * c_1DTestPointCount + pointIndex] = std::abs(yAvg - c_actualValue);
				}
			}
		}
	);

	// non smooth tests
	DoTests(
		"Non Smooth",
		noiseTypes,
		c_1DTestCount,
		c_1DTestPointCount,
		"out/1DResultsNonSmooth.csv",
		[&](int testIndex, pcg32_random_t& rng)
		{
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
				Noise<1>& noise = noiseTypes[noiseIndex];

				// for each sample count in the test
				std::vector<float1> samples;
				for (int pointIndex = 0; pointIndex < c_1DTestPointCount; ++pointIndex)
				{
					// generate the samples
					samples = noise.Generate(rng, pointIndex + 1, samples);

					// integrate!
					float yAvg = 0.0f;
					for (int sampleIndex = 0; sampleIndex < pointIndex + 1; ++sampleIndex)
					{
						float y = EvaluateLinearPiecewise(points, samples[sampleIndex][0]);
						yAvg = Lerp(yAvg, y, 1.0f / float(sampleIndex + 1));
					}

					// store the error
					noise.error[testIndex * c_1DTestPointCount + pointIndex] = std::abs(yAvg - c_actualValue);
				}
			}
		}
	);

	// write out example sample points
	{
		pcg32_random_t rng = GetRNG(0);

		// generate the noise types
		std::vector<std::vector<float1>> noiseSamplePoints(_countof(noiseTypes));
		for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
		{
			std::vector<float1> samples;
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
				fprintf(file, ",\"%f\"", noiseSamplePoints[noiseIndex][i][0]);
			fprintf(file, "\n");
		}

		fclose(file);
	}
}
