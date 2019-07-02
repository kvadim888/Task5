#include "libfix.h"

int32_t	float_to_fix(double num)
{
	int32_t	fixed;
	if (num < -2)
		return INT32_MIN;
	if (num >= 2)
		return INT32_MAX;
	return num * SCALE31;
}

double	fix_to_float(int32_t num)
{
	return num / SCALE31;
}

int32_t	fix_saturate(int64_t num)
{
	t_sample tmp;
	if (num > (int64_t)INT32_MAX)
		return (INT32_MAX);
	if (num < (int64_t)INT32_MIN)
		return (INT32_MIN);
	tmp.int64 = num;
	return tmp.int32[0];
}

int32_t	fix_round(int64_t num)
{
	t_sample tmp;
	tmp.int64 = num + ((uint64_t)1 << PRECISION31);
	return tmp.int32[1];
}

int32_t fix_add(int32_t a, int32_t b)
{
	int64_t sum = (int64_t)a + b;
	return fix_saturate(sum);
}

int32_t fix_sub(int32_t a, int32_t b)
{
	int64_t sub = (int64_t)a - b;
	return fix_saturate(sub);
}

int32_t fix_mul(int32_t a, int32_t b)
{
	int64_t	mul = ((int64_t)a * b) << 1;
	return fix_round(mul);
}

int32_t fix_mac(int64_t *acc, int32_t a, int32_t b)
 {
	*acc += ((int64_t)a * b) << 1;
 	return fix_round(*acc);
 }
 
int32_t fix_msub(int64_t *acc, int32_t a, int32_t b)
 {
	*acc -= ((int64_t)a * b) << 1;
 	return fix_round(*acc);
 }

int32_t	fix_leftshift(int32_t num, int8_t shift)
{
	int64_t res = num;
	return fix_saturate(res << shift);
}

int32_t	fix_rightshift(int32_t num, int8_t shift)
{
	t_sample res;
	res.int32[1] = num;
	return fix_round(res.int64 >> shift);
}
