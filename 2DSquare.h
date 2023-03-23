#pragma once

#include "shared.h"

void Do2DSquareTests()
{
	printf("==================== 2D Square ====================\n");

	Noise<2> noiseTypes[] =
	{
		{ "White", Generate_White<2> },
	};

	// smooth tests
	DoTests(
		"Smooth",
		noiseTypes,
		c_1DTestCount,
		c_1DTestPointCount,
		"out/2DSquareResultsSmooth.csv",
		[&](int testIndex, pcg32_random_t& rng)
		{
		}
	);
}
