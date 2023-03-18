#define _CRT_SECURE_NO_WARNINGS

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

#include "pcg/pcg_basic.h"

#include <array>
#include <vector>
#include <string>
#include <random>

#define DETERMINISTIC() false

static const int c_1DTestCount = 100;
static const int c_1DTestPointCount = 1000;
static const float c_1DTestControlPointMin = -10.0f;
static const float c_1DTestControlPointMax = 10.0f;

typedef std::array<float, 2> float2;

float Lerp(float A, float B, float t)
{
	return A * (1.0f - t) + B * t;
}

pcg32_random_t GetRNG()
{
	pcg32_random_t rng;

#if DETERMINISTIC()
	pcg32_srandom_r(&rng, 0x1337b337, 0xbeefcafe);
#else
	std::random_device rd;
	pcg32_srandom_r(&rng, rd(), rd());
#endif

	return rng;
}

float RandomFloat01(pcg32_random_t& rng)
{
	return ldexpf((float)pcg32_random_r(&rng), -32);
}

float RandomFloatRange(pcg32_random_t& rng, float min, float max)
{
	return Lerp(min, max, RandomFloat01(rng));
}

// https://blog.demofox.org/2017/10/01/calculating-the-distance-between-points-in-wrap-around-toroidal-space/
float DistanceWrap(const float2& A, const float2& B)
{
	float dx = std::abs(A[0] - B[0]);
	dx = std::min(dx, 1.0f - dx);
	float dy = std::abs(A[1] - B[1]);
	dy = std::min(dy, 1.0f - dy);
	return std::sqrt(dx * dx + dy * dy);
}

// A,B,C,D are control points.
// t is a value between 0 and 1.
float Evaluate1DCubicBezier(float A, float B, float C, float D, float t)
{
	float s = 1.0f - t;
	return
		A * s * s * s * 1.0f +
		B * s * s * t * 3.0f +
		C * s * t * t * 3.0f +
		D * t * t * t * 1.0f;
}

void Do1DTests()
{
	pcg32_random_t rng = GetRNG();

	struct Noise
	{
		const char* label;
	};

	Noise noiseTypes[] =
	{
		{ "Regular - Ends" },
		{ "Regular - Left" },
		{ "Regular - Right" },
		{ "Regular - Center" },
		{ "Regular - Equal" },
		{ "White" },
		{ "Stratified" },
		{ "Blue - Wrap" },
		{ "Blue - No Wrap" },
		{ "Blue - Edge" },
		{ "Blue - Half Edge" },
	};

	for (int testIndex = 0; testIndex < c_1DTestCount; ++testIndex)
	{
		// Generate a random 1d Bezier curve
		float A = RandomFloatRange(rng, c_1DTestControlPointMin, c_1DTestControlPointMax);
		float B = RandomFloatRange(rng, c_1DTestControlPointMin, c_1DTestControlPointMax);
		float C = RandomFloatRange(rng, c_1DTestControlPointMin, c_1DTestControlPointMax);
		float D = RandomFloatRange(rng, c_1DTestControlPointMin, c_1DTestControlPointMax);

		for (int noiseIndex = 0; noiseIndex < _countof(noiseTypes); ++noiseIndex)
		{
			for (int pointIndex = 0; pointIndex < c_1DTestPointCount; ++pointIndex)
			{

			}
		}
	}
}

int main(int argc, char** argv)
{
	return 0;
}

/*
TODO:

----- 1D sampling -----

* Sample Types:
 * 0.0, 1.0 style
 * 0.0, 0.5 style
 * 0.5, 1.0 style
 * 0.25, 0.75 style
 * 0.3, 0.6 style
 * white noise?
 * stratified?
 * golden ratio?
 * MBC? (EBC?) with and without wrap around distance

* Integrate random Bezier functions
* show RMSE at each sample

* make CSV
* draw some low point count numberline for each type of noise, to show it on the blog post.

----- 2D sampling -----

* Sample Types:
 * TODO: fill out

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