#include "radar_basic.h"

void pulse_train(float tau, float pri, float p_peak,
		float *dt, float *prf, float *pav, float *ep, float *ru)
{
	*dt = tau / pri;
	*prf = 1.0 / pri;
	*pav = p_peak * (*dt);
	*ep = p_peak * tau;
	*ru = 0.001 * LIGHT_SPEED * pri / 2.0;
}
