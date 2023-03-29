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

/*
TODO:
* break the LDS into their own graph. keep stratified and white to have them common to all graphs
* also probably should have some function for 2d square and circle that isn't separable.
 * like sine of distance maybe.

Notes:
* this started it: https://mastodon.gamedev.place/@jkaniarz/110032776950329500
* "ok so if you do 2 samples on the unit numberline, you can do it at 0, 0.5.  or 0.5 and 1. or you can center them, which gives you the 1/4, 3/4 setup."
* show the math for how easy it is to integrate tensor product bezier surfaces.
* may be better to find factors of numPoints instead of doing the square root dance. not doing it here though.
* 2d circle stratified stratifies the square, then maps square to circle
* 2d regular and hex grid clip to the circle, and then add random points in circle to fill up the rest.
* Note that for more complex shapes, an SDF would work for blue noise candidate scoring

*/