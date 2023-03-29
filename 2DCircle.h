#pragma once

#include "shared.h"
#include "1D.h"
#include "2DSquare.h"
#include <algorithm>

std::vector<float2> GenerateCircle_White(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	std::vector<float2> ret(numSamples);
	for (int i = 0; i < numSamples; ++i)
	{
		do
		{
			ret[i][0] = RandomFloat01(rng);
			ret[i][1] = RandomFloat01(rng);
		}
		while (Distance(ret[i], float2{0.5f, 0.5f}) > 0.5f);
	}
	return ret;
}

std::vector<float2> GenerateCircle_RegularGrid(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	// iteratively reject points in the circle and increase the number of points generated
	// until we have the right number of points in the circle
	std::vector<float2> ret;
	int numSamplesOrig = numSamples;
	do
	{
		ret = Generate2D_RegularGrid(rng, numSamples, ret);

		// remove any points outside of the circle
		ret.erase(std::remove_if(ret.begin(), ret.end(),
			[](const float2& p)
			{
				return Distance(p, float2{ 0.5f, 0.5f }) > 0.5f;
			}
		), ret.end());

		numSamples++;
	}
	while(ret.size() < numSamplesOrig);

	// truncate any extra samples
	ret.resize(numSamplesOrig);
	return ret;
}

std::vector<float2> GenerateCircle_RegularGridCircle(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	std::vector<float2> ret = Generate2D_RegularGrid(rng, numSamples, lastSamples);

	// convert from square to circle
	for (float2& p : ret)
	{
		float theta = p[0] * 2.0f * c_pi;
		float radius = std::sqrt(p[1]);

		p[0] = 0.5f + std::cos(theta) * 0.5f * radius;
		p[1] = 0.5f + std::sin(theta) * 0.5f * radius;
	}

	return ret;
}

std::vector<float2> GenerateCircle_HexGrid(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	// iteratively reject points in the circle and increase the number of points generated
	// until we have the right number of points in the circle
	std::vector<float2> ret;
	int numSamplesOrig = numSamples;
	do
	{
		ret = Generate2D_HexGrid(rng, numSamples, ret);

		// remove any points outside of the circle
		ret.erase(std::remove_if(ret.begin(), ret.end(),
			[](const float2& p)
			{
				return Distance(p, float2{ 0.5f, 0.5f }) > 0.5f;
			}
		), ret.end());

		numSamples++;
	}
	while(ret.size() < numSamplesOrig);

	// truncate any extra samples
	ret.resize(numSamplesOrig);
	return ret;
}

std::vector<float2> GenerateCircle_HexGridCircle(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	std::vector<float2> ret = Generate2D_HexGrid(rng, numSamples, lastSamples);

	// convert from square to circle
	for (float2& p : ret)
	{
		float theta = p[0] * 2.0f * c_pi;
		float radius = std::sqrt(p[1]);

		p[0] = 0.5f + std::cos(theta) * 0.5f * radius;
		p[1] = 0.5f + std::sin(theta) * 0.5f * radius;
	}

	return ret;
}

std::vector<float2> GenerateCircle_Stratified(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	std::vector<float2> ret = Generate2D_Stratified(rng, numSamples, lastSamples);

	// convert square to circle
	for (float2& p : ret)
	{
		float theta = p[0] * 2.0f * c_pi;
		float radius = std::sqrt(p[1]);

		p[0] = 0.5f + std::cos(theta) * 0.5f * radius;
		p[1] = 0.5f + std::sin(theta) * 0.5f * radius;
	}

	return ret;
}

std::vector<float2> GenerateCircle_StratifiedCircle(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	std::vector<float2> ret = Generate2D_Stratified(rng, numSamples, lastSamples);

	// convert from square to circle
	for (float2& p : ret)
	{
		float theta = p[0] * 2.0f * c_pi;
		float radius = std::sqrt(p[1]);

		p[0] = 0.5f + std::cos(theta) * 0.5f * radius;
		p[1] = 0.5f + std::sin(theta) * 0.5f * radius;
	}

	return ret;
}

std::vector<float2> GenerateCircle_Fibonacci(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	std::vector<float2> ret(numSamples);

	for (int i = 0; i < numSamples; ++i)
	{
		float theta = 2.0f * c_pi * std::fmodf(0.5f + float(i) * c_goldenRatioConjugate, 1.0f);
		float radius = std::sqrt((float(i)+0.5f) / float(numSamples));
		ret[i] = {0.5f + std::cos(theta) * 0.5f * radius, 0.5f + std::sin(theta) * 0.5f * radius};
	}

	return ret;
}

std::vector<float2> GenerateCircle_Blue_NoWrap(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
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
			float theta = RandomFloat01(rng) * 2.0f * c_pi;
			float radius = std::sqrt(RandomFloat01(rng));
			float2 candidate = { 0.5f + std::cos(theta) * 0.5f * radius, 0.5f + std::sin(theta) * 0.5f * radius };

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

std::vector<float2> GenerateCircle_Blue_NoWrap_Edge(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
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
			float theta = RandomFloat01(rng) * 2.0f * c_pi;
			float radius = std::sqrt(RandomFloat01(rng));
			float2 candidate = { 0.5f + std::cos(theta) * 0.5f * radius, 0.5f + std::sin(theta) * 0.5f * radius };

			// initialize the score to be the distance to the edge
			float candidateScore = 0.5f - radius * 0.5f;
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

std::vector<float2> GenerateCircle_Blue_NoWrap_HalfEdge(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
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
			float theta = RandomFloat01(rng) * 2.0f * c_pi;
			float radius = std::sqrt(RandomFloat01(rng));
			float2 candidate = { 0.5f + std::cos(theta) * 0.5f * radius, 0.5f + std::sin(theta) * 0.5f * radius };

			// initialize the score to be twice the distance to the edge (make distance to edge count half as much)
			float candidateScore = 2.0f * (0.5f - radius * 0.5f);
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

std::vector<float2> GenerateCircle_R2(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	static const float g = 1.32471795724474602596f;
	static const float a1 = 1.0f / g;
	static const float a2 = 1.0f / (g * g);

	float lastX = 0.5f;
	float lastY = 0.5f;

	std::vector<float2> ret;

	while(ret.size() < numSamples)
	{
		if (Distance(float2{ lastX, lastY }, float2{0.5f, 0.5f}) < 0.5f)
			ret.push_back(float2{ lastX, lastY });

		lastX = std::fmodf(lastX + a1, 1.0f);
		lastY = std::fmodf(lastY + a2, 1.0f);
	}

	return ret;
}

std::vector<float2> GenerateCircle_R2Circle(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	std::vector<float2> ret = Generate2D_R2(rng, numSamples, lastSamples);

	// convert from square to circle
	for (float2& p : ret)
	{
		float theta = p[0] * 2.0f * c_pi;
		float radius = std::sqrt(p[1]);

		p[0] = 0.5f + std::cos(theta) * 0.5f * radius;
		p[1] = 0.5f + std::sin(theta) * 0.5f * radius;
	}

	return ret;
}

std::vector<float2> GenerateCircle_Halton23(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	std::vector<float2> ret;

	int index = 1;
	while (ret.size() < numSamples)
	{
		float2 p = {
			Halton(index, 2),
			Halton(index, 3)
		};
		index++;

		if (Distance(p, float2{ 0.5f, 0.5f }) < 0.5f)
			ret.push_back(p);
	}

	return ret;
}

std::vector<float2> GenerateCircle_Halton23Circle(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	std::vector<float2> ret = Generate2D_Halton23(rng, numSamples, lastSamples);

	// convert from square to circle
	for (float2& p : ret)
	{
		float theta = p[0] * 2.0f * c_pi;
		float radius = std::sqrt(p[1]);

		p[0] = 0.5f + std::cos(theta) * 0.5f * radius;
		p[1] = 0.5f + std::sin(theta) * 0.5f * radius;
	}

	return ret;
}

std::vector<float2> GenerateCircle_Sobol(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	std::vector<float2> ret;

	int index = 0;
	while (ret.size() < numSamples)
	{
		uint4 s = Sobol(index);

		float2 p = {
			float(s[0]) / pow(2.0f, 32.0f),
			float(s[1]) / pow(2.0f, 32.0f)
		};
		index++;

		if (Distance(p, float2{ 0.5f, 0.5f }) < 0.5f)
			ret.push_back(p);
	}

	return ret;
}

std::vector<float2> GenerateCircle_SobolCircle(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	std::vector<float2> ret = Generate2D_Sobol(rng, numSamples, lastSamples);

	// convert from square to circle
	for (float2& p : ret)
	{
		float theta = p[0] * 2.0f * c_pi;
		float radius = std::sqrt(p[1]);

		p[0] = 0.5f + std::cos(theta) * 0.5f * radius;
		p[1] = 0.5f + std::sin(theta) * 0.5f * radius;
	}

	return ret;
}

std::vector<float2> GenerateCircle_BurleySobol(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	std::vector<float2> ret;

	uint32_t seed = pcg32_random_r(&rng);

	int index = 0;
	while (ret.size() < numSamples)
	{
		uint4 s = ShuffledScrambledSobol(index, seed);

		float2 p = {
			float(s[0]) / pow(2.0f, 32.0f),
			float(s[1]) / pow(2.0f, 32.0f)
		};
		index++;

		if (Distance(p, float2{ 0.5f, 0.5f }) < 0.5f)
			ret.push_back(p);
	}

	return ret;
}

std::vector<float2> GenerateCircle_BurleySobolCircle(pcg32_random_t& rng, int numSamples, std::vector<float2>& lastSamples)
{
	std::vector<float2> ret = Generate2D_BurleySobol(rng, numSamples, lastSamples);

	// convert from square to circle
	for (float2& p : ret)
	{
		float theta = p[0] * 2.0f * c_pi;
		float radius = std::sqrt(p[1]);

		p[0] = 0.5f + std::cos(theta) * 0.5f * radius;
		p[1] = 0.5f + std::sin(theta) * 0.5f * radius;
	}

	return ret;
}

void Do2DCircleTests()
{
	printf("==================== 2D Circle ====================\n");

	Noise<2> noiseTypes[] =
	{
		{ "White", GenerateCircle_White },
		{ "Regular Grid", GenerateCircle_RegularGrid },
		{ "Regular Grid Circle", GenerateCircle_RegularGridCircle },
		{ "Hex Grid", GenerateCircle_HexGrid },
		{ "Hex Grid Circle", GenerateCircle_HexGridCircle },
		{ "Stratified", GenerateCircle_Stratified },
		{ "Stratified Circle", GenerateCircle_StratifiedCircle },
		{ "R2", GenerateCircle_R2 },
		{ "R2 Circle", GenerateCircle_R2Circle },
		{ "Halton23", GenerateCircle_Halton23 },
		{ "Halton23 Circle", GenerateCircle_Halton23Circle },
		{ "Sobol", GenerateCircle_Sobol},
		{ "Sobol Circle", GenerateCircle_SobolCircle },
		{ "Burley Sobol", GenerateCircle_BurleySobol},
		{ "Burley Sobol Circle", GenerateCircle_BurleySobolCircle },
		{ "Fibonacci", GenerateCircle_Fibonacci },
		{ "Blue - No Wrap", GenerateCircle_Blue_NoWrap },
		{ "Blue - No Wrap Edge", GenerateCircle_Blue_NoWrap_Edge },
		{ "Blue - No Wrap Half Edge", GenerateCircle_Blue_NoWrap_HalfEdge },
	};

	// smooth tests
	DoTests(
		"Smooth",
		noiseTypes,
		c_2DTestCount,
		c_2DTestPointCount,
		"out/2DCircleResultsSmooth.csv",
		[&](int testIndex, pcg32_random_t& rng)
		{
			// Generate a 4x4 grid of random numbers to be the control points of our biquadratic bezier surface
			float controlPoints[16];
			for (int i = 0; i < 16; ++i)
				controlPoints[i] = RandomFloatRange(rng, c_2DTestControlPointMin, c_2DTestControlPointMax);

			// Calculate actual value through monte carlo integration
			float actualValue = 0.0f;
			for (int i = 0; i < c_2DCircleActualValueSamples; ++i)
			{
				float theta = 2.0f * c_pi * RandomFloat01(rng);
				float radius = std::sqrt(RandomFloat01(rng));

				float2 p = float2{
					0.5f + std::cos(theta) * 0.5f * radius,
					0.5f + std::sin(theta) * 0.5f * radius
				};

				float controlPointsY[4] =
				{
					Evaluate1DCubicBezier(controlPoints[0], controlPoints[1], controlPoints[2], controlPoints[3], p[0]),
					Evaluate1DCubicBezier(controlPoints[4], controlPoints[5], controlPoints[6], controlPoints[7], p[0]),
					Evaluate1DCubicBezier(controlPoints[8], controlPoints[9], controlPoints[10], controlPoints[11], p[0]),
					Evaluate1DCubicBezier(controlPoints[12], controlPoints[13], controlPoints[14], controlPoints[15], p[0])
				};

				float y = Evaluate1DCubicBezier(controlPointsY[0], controlPointsY[1], controlPointsY[2], controlPointsY[3], p[1]);

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
					noise.error[testIndex * c_2DTestPointCount + pointIndex] = std::abs(yAvg - actualValue);
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
		"out/2DCircleResultsNonSmooth.csv",
		[&](int testIndex, pcg32_random_t& rng)
		{
			// Generate 2x2 control points for each bilinear patch in a 2x2 grid, for a total of 16 control points.
			float controlPoints[16];
			for (int i = 0; i < 16; ++i)
				controlPoints[i] = RandomFloatRange(rng, c_2DTestControlPointMin, c_2DTestControlPointMax);

			// Calculate actual value through monte carlo integration
			float actualValue = 0.0f;
			for (int i = 0; i < c_2DCircleActualValueSamples; ++i)
			{
				float theta = 2.0f * c_pi * RandomFloat01(rng);
				float radius = std::sqrt(RandomFloat01(rng));

				float2 p = float2{
					0.5f + std::cos(theta) * 0.5f * radius,
					0.5f + std::sin(theta) * 0.5f * radius
				};

				float cellXf = p[0] * 2.0f;
				float cellYf = p[1] * 2.0f;

				int cellXi = std::min(int(cellXf), 1);
				int cellYi = std::min(int(cellYf), 1);

				float fractX = cellXf - float(cellXi);
				float fractY = cellYf - float(cellYi);

				int cpOffset = (cellYi * 2 + cellXi) * 4;

				float cpX0 = Lerp(controlPoints[cpOffset + 0], controlPoints[cpOffset + 1], fractX);
				float cpX1 = Lerp(controlPoints[cpOffset + 2], controlPoints[cpOffset + 3], fractX);

				float y = Lerp(cpX0, cpX1, fractY);

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
					noise.error[testIndex * c_2DTestPointCount + pointIndex] = std::abs(yAvg - actualValue);
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
		"out/2DCircleResultsNonSeparable.csv",
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
				float theta = 2.0f * c_pi * RandomFloat01(rng);
				float radius = std::sqrt(RandomFloat01(rng));

				float2 p = float2{
					0.5f + std::cos(theta) * 0.5f * radius,
					0.5f + std::sin(theta) * 0.5f * radius
				};

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
		fopen_s(&file, "out/2DCirclePoints.csv", "wb");

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
