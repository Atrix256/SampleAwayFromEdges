#define _CRT_SECURE_NO_WARNINGS

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb/stb_image_write.h"

#include "pcg/pcg_basic.h"

#include <array>
#include <vector>

#define DETERMINISTIC() false

typedef std::array<float, 2> float2;

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

// https://blog.demofox.org/2017/10/01/calculating-the-distance-between-points-in-wrap-around-toroidal-space/
float DistanceWrap(const float2& A, const float2& B)
{
	float dx = std::abs(A[0] - B[0]);
	dx = std::min(dx, 1.0f - dx);
	float dy = std::abs(A[1] - B[1]);
	dy = std::min(dy, 1.0f - dy);
	return std::sqrt(dx * dx + dy * dy);
}

double Lerp(double A, double B, double t)
{
	return A * (1.0 - t) + B * t;
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