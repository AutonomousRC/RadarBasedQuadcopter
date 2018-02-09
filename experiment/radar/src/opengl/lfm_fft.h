#ifndef __LFM_FFT_H__
#define __LFM_FFT_H__

#include <time.h>

#include <GL/glut.h>
#include <GL/gl.h>

void lfm_fft(void);

typedef struct complex
{
        double re;
        double im;
} c;

#endif
