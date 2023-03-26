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

	// TODO: take the above sqrt out and verify it is working correctly as is.

	return ret;
}

void Do2DCircleTests()
{
	printf("==================== 2D Circle ====================\n");

	Noise<2> noiseTypes[] =
	{
		{ "White", GenerateCircle_White },
		{ "Regular Grid", GenerateCircle_RegularGrid },
		{ "Hex Grid", GenerateCircle_HexGrid },
		{ "Stratified", GenerateCircle_Stratified },

/*
		{ "Blue - Wrap", Generate2D_Blue_Wrap },
		{ "Blue - No Wrap", Generate2D_Blue_NoWrap },
		{ "Blue - No Wrap Edge", Generate2D_Blue_NoWrap_Edge },
		{ "Blue - No Wrap Half Edge", Generate2D_Blue_NoWrap_HalfEdge },
		*/
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
		"out/2DCircleResultsNonSmooth.csv",
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
