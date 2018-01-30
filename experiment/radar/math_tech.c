#include <math.h>
#include "math_tech.h"
#include "radar_basic.h"

float angle2radian(float ang)
{
        return ang * M_PI / 180.0;
}

float get_wavelength(float freq)
{
        return LIGHT_SPEED / freq;
}
