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

static const int	c_2DTestCount = 1000;
static const int	c_2DTestPointCount = 200;
static const float	c_2DTestControlPointMin = -10.0f;
static const float	c_2DTestControlPointMax = 10.0f;
static const int	c_2DNumPointsReported = 256;

static const int	c_2DCircleActualValueSamples = 1000000;

#include "1D.h"
#include "2DSquare.h"
#include "2DCircle.h"

int main(int argc, char** argv)
{
	_mkdir("out");

	Do1DTests();
	Do2DSquareTests();
	Do2DCircleTests();

	return 0;
}
