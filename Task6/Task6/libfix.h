#ifndef LIBFIX_H
#define LIBFIX_H

#include <stdint.h>
#include <stdio.h>

#include <assert.h>

#define PRECISION31		31 
#define PRECISION30		30 

#define SCALE31			(double)(1LL << PRECISION31)
#define SCALE30			(double)(1LL << PRECISION30)

typedef union
{
	uint64_t	int64;
	uint32_t	int32[2];
	uint16_t	int16[4];
	uint8_t		int8[8];
}				t_sample;

int32_t	float_to_fix(double num);
double	fix_to_float(int32_t num);

int32_t	fix_saturate(int64_t num);
int32_t	fix_round(int64_t num);

int32_t fix_add(int32_t a, int32_t b);
int32_t fix_sub(int32_t a, int32_t b);
int32_t	fix_mul(int32_t a, int32_t b);

int32_t	fix_mac(int64_t *acc, int32_t a, int32_t b);
int32_t	fix_msub(int64_t *acc, int32_t a, int32_t b);

int32_t	fix_leftshift(int32_t num, int8_t shift);
int32_t	fix_rightshift(int32_t num, int8_t shift);

#endif	//LIBFIX_H
