#pragma once

#include "shared.h"
#include "1D.h"

std::vector<float2> Generate2D_RegularGrid(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	std::vector<float2> ret(numSamples);

	// Do the full nearly square grid that is <= the number of points we want to sample.
	int sideX = int(std::sqrt(numSamples));
	int sideY = numSamples / sideX;
	int gridSamples = sideX * sideY;
	for (int i = 0; i < gridSamples; ++i)
	{
		ret[i][0] = (float(i % sideX) + 0.5f) / float(sideX);
		ret[i][1] = (float(i / sideX) + 0.5f) / float(sideY);
	}

	// Do the remainder of the samples as white noise.
	for (int i = gridSamples; i < numSamples; ++i)
	{
		ret[i][0] = RandomFloat01(rng);
		ret[i][1] = RandomFloat01(rng);
	}

	return ret;
}

std::vector<float2> Generate2D_HexGrid(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	std::vector<float2> ret(numSamples);

	// Do the full nearly square grid that is <= the number of points we want to sample.
	int sideX = int(std::sqrt(numSamples));
	int sideY = numSamples / sideX;
	int gridSamples = sideX * sideY;
	for (int i = 0; i < gridSamples; ++i)
	{
		int y = i / sideX;

		float xOffset = 0.25f + ((y % 2) ? 0.5f : 0.0f);
		ret[i][0] = (float(i % sideX) + xOffset) / float(sideX);
		ret[i][1] = (float(i / sideX) + 0.5f) / float(sideY);
	}

	// Do the remainder of the samples as white noise.
	for (int i = gridSamples; i < numSamples; ++i)
	{
		ret[i][0] = RandomFloat01(rng);
		ret[i][1] = RandomFloat01(rng);
	}

	return ret;
}

std::vector<float2> Generate2D_Stratified(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	std::vector<float2> ret(numSamples);

	// Do the full nearly square grid that is <= the number of points we want to sample.
	int sideX = int(std::sqrt(numSamples));
	int sideY = numSamples / sideX;
	int gridSamples = sideX * sideY;
	for (int i = 0; i < gridSamples; ++i)
	{
		ret[i][0] = (float(i % sideX) + RandomFloat01(rng)) / float(sideX);
		ret[i][1] = (float(i / sideX) + RandomFloat01(rng)) / float(sideY);
	}

	// Do the remainder of the samples as white noise.
	for (int i = gridSamples; i < numSamples; ++i)
	{
		ret[i][0] = RandomFloat01(rng);
		ret[i][1] = RandomFloat01(rng);
	}

	return ret;
}

std::vector<float2> Generate2D_Fibonacci(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	std::vector<float2> ret(numSamples);
	for (int i = 0; i < numSamples; ++i)
	{
		ret[i] = {
			std::fmodf(0.5f + float(i) * c_goldenRatioConjugate, 1.0f),
			(float(i) + 0.5f) / float(numSamples)
		};
	}
	return ret;
}

std::vector<float2> Generate2D_Blue_Wrap(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	std::vector<float2> ret;
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
		float2 bestCandidate = { 0.0f, 0.0f };
		int candidateCount = sampleIndex + 1;
		for (int candidateIndex = 0; candidateIndex < candidateCount; ++candidateIndex)
		{
			float2 candidate = { RandomFloat01(rng), RandomFloat01(rng) };
			float candidateScore = FLT_MAX;
			for (int pointIndex = 0; pointIndex < sampleIndex; ++pointIndex)
				candidateScore = std::min(candidateScore, DistanceWrap(candidate, ret[pointIndex]));

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

std::vector<float2> Generate2D_Blue_NoWrap(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	std::vector<float2> ret;
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
		float2 bestCandidate = { 0.0f, 0.0f };
		int candidateCount = sampleIndex + 1;
		for (int candidateIndex = 0; candidateIndex < candidateCount; ++candidateIndex)
		{
			float2 candidate = { RandomFloat01(rng), RandomFloat01(rng) };
			float candidateScore = FLT_MAX;
			for (int pointIndex = 0; pointIndex < sampleIndex; ++pointIndex)
				candidateScore = std::min(candidateScore, Distance(candidate, ret[pointIndex]));

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

std::vector<float2> Generate2D_Blue_NoWrap_Edge(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	std::vector<float2> ret;
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
		float2 bestCandidate = { 0.0f, 0.0f };
		int candidateCount = sampleIndex + 1;
		for (int candidateIndex = 0; candidateIndex < candidateCount; ++candidateIndex)
		{
			float2 candidate = { RandomFloat01(rng), RandomFloat01(rng) };

			// initialize the score to be the distance to the edge
			float distEdgeX = std::min(candidate[0], 1.0f - candidate[0]);
			float distEdgeY = std::min(candidate[1], 1.0f - candidate[1]);
			float candidateScore = std::min(distEdgeX, distEdgeY);

			for (int pointIndex = 0; pointIndex < sampleIndex; ++pointIndex)
				candidateScore = std::min(candidateScore, Distance(candidate, ret[pointIndex]));

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

std::vector<float2> Generate2D_Blue_NoWrap_HalfEdge(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	std::vector<float2> ret;
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
		float2 bestCandidate = { 0.0f, 0.0f };
		int candidateCount = sampleIndex + 1;
		for (int candidateIndex = 0; candidateIndex < candidateCount; ++candidateIndex)
		{
			float2 candidate = { RandomFloat01(rng), RandomFloat01(rng) };

			// initialize the score to be twice the distance to the edge (make distance to edge count half as much)
			float distEdgeX = std::min(candidate[0], 1.0f - candidate[0]);
			float distEdgeY = std::min(candidate[1], 1.0f - candidate[1]);
			float candidateScore = 2.0f * std::min(distEdgeX, distEdgeY);

			for (int pointIndex = 0; pointIndex < sampleIndex; ++pointIndex)
				candidateScore = std::min(candidateScore, Distance(candidate, ret[pointIndex]));

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

std::vector<float2> Generate2D_R2(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	static const float g = 1.32471795724474602596f;
	static const float a1 = 1.0f / g;
	static const float a2 = 1.0f / (g * g);

	float lastX = 0.5f;
	float lastY = 0.5f;

	std::vector<float2> ret(numSamples);

	for (float2& p : ret)
	{
		p = float2{ lastX, lastY };
		lastX = std::fmodf(lastX + a1, 1.0f);
		lastY = std::fmodf(lastY + a2, 1.0f);
	}

	return ret;
}

std::vector<float2> Generate2D_Halton23(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	std::vector<float2> ret(numSamples);
	for (int i = 0; i < numSamples; ++i)
	{
		ret[i][0] = Halton(i + 1, 2);
		ret[i][1] = Halton(i + 1, 3);
	}
	return ret;
}

std::vector<float2> Generate2D_Sobol(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	std::vector<float2> ret(numSamples);
	for (int i = 0; i < numSamples; ++i)
	{
		uint4 s = Sobol(i);
		ret[i][0] = float(s[0]) / pow(2.0f, 32.0f);
		ret[i][1] = float(s[1]) / pow(2.0f, 32.0f);
	}
	return ret;
}

std::vector<float2> Generate2D_BurleySobol(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	uint32_t seed = pcg32_random_r(&rng);

	std::vector<float2> ret(numSamples);
	for (int i = 0; i < numSamples; ++i)
	{
		uint4 s = ShuffledScrambledSobol(i, seed);
		ret[i][0] = float(s[0]) / pow(2.0f, 32.0f);
		ret[i][1] = float(s[1]) / pow(2.0f, 32.0f);
	}
	return ret;
}

void Do2DSquareTests()
{
	printf("==================== 2D Square ====================\n");

	Noise<2> noiseTypes[] =
	{
		{ "White", Generate_White<2> },
		{ "Regular Grid", Generate2D_RegularGrid },
		{ "Hex Grid", Generate2D_HexGrid },
		{ "Stratified", Generate2D_Stratified },
		{ "Fibonacci", Generate2D_Fibonacci },
		{ "R2", Generate2D_R2 },
		{ "Halton23", Generate2D_Halton23 },
		{ "Sobol", Generate2D_Sobol },
		{ "Burley Sobol", Generate2D_BurleySobol },
		{ "Blue - Wrap", Generate2D_Blue_Wrap },
		{ "Blue - No Wrap", Generate2D_Blue_NoWrap },
		{ "Blue - No Wrap Edge", Generate2D_Blue_NoWrap_Edge },
		{ "Blue - No Wrap Half Edge", Generate2D_Blue_NoWrap_HalfEdge },
	};

	// smooth tests
	DoTests(
		"Smooth",
		noiseTypes,
		c_2DTestCount,
		c_2DTestPointCount,
		"out/2DSquareResultsSmooth.csv",
		[&](int testIndex, pcg32_random_t& rng)
		{
			// Generate a 4x4 grid of random numbers to be the control points of our bicubic bezier surface
			float controlPoints[16];
			for (int i = 0; i < 16; ++i)
				controlPoints[i] = RandomFloatRange(rng, c_2DTestControlPointMin, c_2DTestControlPointMax);

			// Calculate the definite integral of the bezier surface, in [0,1]^2
			const float c_actualValueControlPoints[4] =
			{
				Integral1DCubicBezier(controlPoints[0], controlPoints[1], controlPoints[2], controlPoints[3]),
				Integral1DCubicBezier(controlPoints[4], controlPoints[5], controlPoints[6], controlPoints[7]),
				Integral1DCubicBezier(controlPoints[8], controlPoints[9], controlPoints[10], controlPoints[11]),
				Integral1DCubicBezier(controlPoints[12], controlPoints[13], controlPoints[14], controlPoints[15])
			};
			const float c_actualValue = Integral1DCubicBezier(c_actualValueControlPoints[0], c_actualValueControlPoints[1], c_actualValueControlPoints[2], c_actualValueControlPoints[3]);

			// for each type of noise
			for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
			{
				Noise<2>& noise = noiseTypes[noiseIndex];

				// for each sample count in the test
				std::vector<float2> samples;
				for (int pointIndex = 0; pointIndex < c_2DTestPointCount; ++pointIndex)
				{
					// generate the samples
					samples = noise.Generate(rng, pointIndex + 1, samples);

					// integrate!
					float yAvg = 0.0f;
					for (int sampleIndex = 0; sampleIndex < pointIndex + 1; ++sampleIndex)
					{
						float controlPointsY[4] =
						{
							Evaluate1DCubicBezier(controlPoints[0], controlPoints[1], controlPoints[2], controlPoints[3], samples[sampleIndex][0]),
							Evaluate1DCubicBezier(controlPoints[4], controlPoints[5], controlPoints[6], controlPoints[7], samples[sampleIndex][0]),
							Evaluate1DCubicBezier(controlPoints[8], controlPoints[9], controlPoints[10], controlPoints[11], samples[sampleIndex][0]),
							Evaluate1DCubicBezier(controlPoints[12], controlPoints[13], controlPoints[14], controlPoints[15], samples[sampleIndex][0])
						};

						float y = Evaluate1DCubicBezier(controlPointsY[0], controlPointsY[1], controlPointsY[2], controlPointsY[3], samples[sampleIndex][1]);
						yAvg = Lerp(yAvg, y, 1.0f / float(sampleIndex + 1));
					}

					// store the error
					noise.error[testIndex * c_2DTestPointCount + pointIndex] = std::abs(yAvg - c_actualValue);
				}
			}
		}
	);

	// Non smooth tests
	DoTests(
		"Non Smooth",
		noiseTypes,
		c_2DTestCount,
		c_2DTestPointCount,
		"out/2DSquareResultsNonSmooth.csv",
		[&](int testIndex, pcg32_random_t& rng)
		{
			// Generate 2x2 control points for each bilinear patch in a 2x2 grid, for a total of 16 control points.
			float controlPoints[16];
			for (int i = 0; i < 16; ++i)
				controlPoints[i] = RandomFloatRange(rng, c_2DTestControlPointMin, c_2DTestControlPointMax);

			// The integral of a bilinear patch is the average of it's 4 control points multiplied by it's size.
			// Integrating muiltiple piecewise bilinear patches is just averaging their values and multiplying by total size.
			// Our square of integration is [0,1]^2 so has area of 1.
			// So, the integral of this 2x2 grid of randomly generated bilinear patches is just the average of
			// all of the control points.

			// Calculate the definite integral of the bezier surface, in [0,1]^2
			const float c_actualValue =
			(
				controlPoints[0] + controlPoints[1] + controlPoints[2] + controlPoints[3] +
				controlPoints[4] + controlPoints[5] + controlPoints[6] + controlPoints[7] +
				controlPoints[8] + controlPoints[9] + controlPoints[10] + controlPoints[11] +
				controlPoints[12] + controlPoints[13] + controlPoints[14] + controlPoints[15]
			) / 16.0f;

			// for each type of noise
			for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
			{
				Noise<2>& noise = noiseTypes[noiseIndex];

				// for each sample count in the test
				std::vector<float2> samples;
				for (int pointIndex = 0; pointIndex < c_2DTestPointCount; ++pointIndex)
				{
					// generate the samples
					samples = noise.Generate(rng, pointIndex + 1, samples);

					// integrate!
					float yAvg = 0.0f;
					for (int sampleIndex = 0; sampleIndex < pointIndex + 1; ++sampleIndex)
					{
						float cellXf = samples[sampleIndex][0] * 2.0f;
						float cellYf = samples[sampleIndex][1] * 2.0f;

						int cellXi = std::min(int(cellXf), 1);
						int cellYi = std::min(int(cellYf), 1);

						float fractX = cellXf - float(cellXi);
						float fractY = cellYf - float(cellYi);

						int cpOffset = (cellYi * 2 + cellXi) * 4;

						float cpX0 = Lerp(controlPoints[cpOffset + 0], controlPoints[cpOffset + 1], fractX);
						float cpX1 = Lerp(controlPoints[cpOffset + 2], controlPoints[cpOffset + 3], fractX);

						float y = Lerp(cpX0, cpX1, fractY);

						yAvg = Lerp(yAvg, y, 1.0f / float(sampleIndex + 1));
					}

					// store the error
					noise.error[testIndex * c_2DTestPointCount + pointIndex] = std::abs(yAvg - c_actualValue);
				}
			}
		}
	);

	// Non separable tests
	DoTests(
		"Non Separable",
		noiseTypes,
		c_2DTestCount,
		c_2DTestPointCount,
		"out/2DSquareResultsNonSeparable.csv",
		[&](int testIndex, pcg32_random_t& rng)
		{
			// Generate a random xy translation, and xy scale.
			float parameters[4];
			for (int i = 0; i < 4; ++i)
				parameters[i] = RandomFloatRange(rng, c_2DTestControlPointMin, c_2DTestControlPointMax);

			auto F = [&](const float2& p)
			{
				float dx = (p[0] - parameters[0]) * parameters[2];
				float dy = (p[1] - parameters[1]) * parameters[3];
				return std::sin(std::sqrt(dx * dx + dy * dy));
			};

			// Calculate actual value through monte carlo integration
			float actualValue = 0.0f;
			for (int i = 0; i < c_2DCircleActualValueSamples; ++i)
			{
				float2 p = float2{ RandomFloat01(rng), RandomFloat01(rng) };
				float y = F(p);
				actualValue = Lerp(actualValue, y, 1.0f / float(i + 1));
			}

			// for each type of noise
			for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
			{
				Noise<2>& noise = noiseTypes[noiseIndex];

				// for each sample count in the test
				std::vector<float2> samples;
				for (int pointIndex = 0; pointIndex < c_2DTestPointCount; ++pointIndex)
				{
					// generate the samples
					samples = noise.Generate(rng, pointIndex + 1, samples);

					// integrate!
					float yAvg = 0.0f;
					for (int sampleIndex = 0; sampleIndex < pointIndex + 1; ++sampleIndex)
					{
						float y = F(samples[sampleIndex]);
						yAvg = Lerp(yAvg, y, 1.0f / float(sampleIndex + 1));
					}

					// store the error
					noise.error[testIndex * c_2DTestPointCount + pointIndex] = std::abs(yAvg - actualValue);
				}
			}
		}
	);

	// write out example sample points
	{
		pcg32_random_t rng = GetRNG(0);

		// generate the noise types
		std::vector<std::vector<float2>> noiseSamplePoints(_countof(noiseTypes));
		for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
		{
			std::vector<float2> samples;
			noiseSamplePoints[noiseIndex] = noiseTypes[noiseIndex].Generate(rng, c_2DNumPointsReported, samples);
		}

		// create the file
		FILE* file = nullptr;
		fopen_s(&file, "out/2DPoints.csv", "wb");

		// write the header
		fprintf(file, "\"samples\"");
		for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
			fprintf(file, ",\"%s X\",\"%s Y\"", noiseTypes[noiseIndex].label, noiseTypes[noiseIndex].label);
		fprintf(file, "\n");

		// write the sample points
		for (int i = 0; i < c_2DNumPointsReported; ++i)
		{
			fprintf(file, "\"%i\"", i + 1);
			for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
				fprintf(file, ",\"%f\",\"%f\"", noiseSamplePoints[noiseIndex][i][0], noiseSamplePoints[noiseIndex][i][1]);
			fprintf(file, "\n");
		}

		fclose(file);
	}
}
