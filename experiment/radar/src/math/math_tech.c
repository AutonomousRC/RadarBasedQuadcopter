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

void slice_section(double s, double e, double itv, double *arr)
{
        int i;
        double tmp = e / itv;

        for(i = 0; i < (int)(tmp + 2); i++)
        {
                arr[i] = itv * i;
#if 0
                printf("arr[%d] = %.9lf\n", i, arr[i]);
#endif
        }
}

void linear_slice(double start, double end, int num, double *arr)
{
        int i;
        double interval = (end - start) / (num - 1);
	//double interval = (end - start) / num;

        for(i = 0; i < num; i++)
                arr[i] = start + interval * i;
}

void pow_vec(double *x, double *y, int data_num, double series)
{
	int i;

	for(i = 0; i < data_num; i++)
		y[i] = pow(x[i], series);
}

void p_arr(double *arr, int num)
{
        int i;

        for(i = 0; i < num; i++)
                printf("arr[%d] = %.12lf\n", i, arr[i]);
}

void swap(double *arr, int left, int right)
{
	int tmp;
	tmp = arr[left];
	arr[left] = arr[right];
	arr[right] = tmp;
}

int partition(double *arr, int left, int right, int pivot_idx)
{
	int pivot_val = arr[pivot_idx];
	int store_idx = left;
	int i;

	swap(arr, pivot_idx, right);

	for(i = left; i < right; i++)
		if(arr[i] <= pivot_val)
		{
			swap(arr, i, store_idx);
			++store_idx;
		}

	swap(arr, store_idx, right);
	return store_idx;
}

void quicksort(double *arr, int left, int right)
{
	int pivot_idx = left;
	int pivot_new_idx;

loop:
	if(right > left)
	{
		pivot_new_idx = partition(arr, left, right, pivot_idx);
		quicksort(arr, left, pivot_new_idx - 1);
		pivot_idx = left = pivot_new_idx + 1;
		goto loop;
	}
}

double find_max(double *data, int num)
{
	quicksort(data, 0, num - 1);

	return data[num - 1];
}

double find_min(double *data, int num)
{
	quicksort(data, 0, num - 1);

	return data[0];
}

void divide_vec(double *data, double div, int num)
{
	int i;

	for(i = 0; i < num; i++)
		data[i] /= div;
}

void volt_base_log(double *data, double eps, int num)
{
	int i;

	for(i = 0; i < num; i++)
		data[i] = 20 * log10(data[i] + eps);
}
