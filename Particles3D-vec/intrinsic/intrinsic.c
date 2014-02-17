#include <immintrin.h>

void foo()
{
	__m512d v2;

	_mm512_reduce_add_pd(v2);
}
