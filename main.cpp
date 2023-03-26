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
static const int	c_2DNumPointsReported = 64;

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

* GenerateCircle_RegularGrid is how all the 2d things should work
* maybe do fibonacci spiral too for point in circle?
* could also do uniform point in circle transformation from square to circle

* calculate the correct actual value for 2d circle
* hide the square in the 2d circle plots and draw a circle to show where the points are constrained to!

? is it meaningful that the hex grid isn't centered in the circle? i think it might be
* 2d circle stratified isn't quite right.

----- 2D sampling -----

* Sample Types:
 * r2, sobol, halton(2,3)?
 
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
* show the math for how easy it is to integrate tensor product bezier surfaces.
* may be better to find factors of numPoints instead of doing the square root dance. not doing it here though.

*/