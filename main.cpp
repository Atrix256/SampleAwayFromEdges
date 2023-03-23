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

#include "1D.h"
#include "2DSquare.h"

int main(int argc, char** argv)
{
	_mkdir("out");

	Do1DTests();

	Do2DSquareTests();

	return 0;
}

/*
TODO:

----- 2D sampling -----

* Sample Types:
 * r2, sobol, halton(2,3)?
 * blue wrap
 * blue no wrap
 * blue edge
 * blue half edge
 * regular grid?
 * hexagon type sampling by half offsetting each row?
 * white
 * stratified

 Three sampling tests total:
 0) 1d tests
 1) bezier in a square
 2) bezier in a shape?
   * rejection sampling for blue and white.
   ? what about stratified? and everything else! r2, sobol, halton?
 ! smooth and non smooth.

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