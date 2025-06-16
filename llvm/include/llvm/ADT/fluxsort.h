#ifndef FLUXSORT_H
#define FLUXSORT_H

#include <assert.h>
#include <errno.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "quadsort.h"
#include "fluxsort.inc"

template <class T>
void fluxsort(T *array, size_t nmemb, size_t size, CMPFUNC<T> cmp) noexcept {
  if (nmemb < 2) {
    return;
  }
  if (nmemb <= 132) {
    quadsort(array, nmemb, size, cmp);
  } else {
    T *pta = (T *)array;
    T swap[nmemb];
    flux_analyze(pta, swap, nmemb, nmemb, cmp);
  }
}

#endif
