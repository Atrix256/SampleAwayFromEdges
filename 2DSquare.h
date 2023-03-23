#pragma once

#include "shared.h"

void Do2DSquareTests()
{
	printf("==================== 2D Square ====================\n");

	Noise<2> noiseTypes[] =
	{
		{ "White", Generate_White<2> },
	};

	// allocate space for the results of our test.
	// store them all out so we can work multithreadedly, then deterministically combine the results together.
	for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
	{
		Noise<2>& noise = noiseTypes[noiseIndex];
		noise.error.resize(c_1DTestCount * c_1DTestPointCount);
		noise.avgErrorAtSamples.resize(c_1DTestPointCount);
		noise.avgErrorSqAtSamples.resize(c_1DTestPointCount);
	}



}
