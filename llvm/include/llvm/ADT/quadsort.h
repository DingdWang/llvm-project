#ifndef QUADSORT_H
#define QUADSORT_H

#include <assert.h>
#include <errno.h>
#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define head_branchless_merge(ptd, x, ptl, ptr, cmp)                           \
  *ptd++ = cmp(ptl, ptr) <= 0 ? *ptl++ : *ptr++;

#define tail_branchless_merge(tpd, x, tpl, tpr, cmp)                           \
  *tpd-- = cmp(tpl, tpr) > 0 ? *tpl-- : *tpr--;

#define parity_merge_two(array, swap, x, ptl, ptr, pts, cmp)                   \
  ptl = array;                                                                 \
  ptr = array + 2;                                                             \
  pts = swap;                                                                  \
  head_branchless_merge(pts, x, ptl, ptr, cmp);                                \
  *pts = cmp(ptl, ptr) <= 0 ? *ptl : *ptr;                                     \
                                                                               \
  ptl = array + 1;                                                             \
  ptr = array + 3;                                                             \
  pts = swap + 3;                                                              \
  tail_branchless_merge(pts, x, ptl, ptr, cmp);                                \
  *pts = cmp(ptl, ptr) > 0 ? *ptl : *ptr;

#define parity_merge_four(array, swap, x, ptl, ptr, pts, cmp)                  \
  ptl = array + 0;                                                             \
  ptr = array + 4;                                                             \
  pts = swap;                                                                  \
  head_branchless_merge(pts, x, ptl, ptr, cmp);                                \
  head_branchless_merge(pts, x, ptl, ptr, cmp);                                \
  head_branchless_merge(pts, x, ptl, ptr, cmp);                                \
  *pts = cmp(ptl, ptr) <= 0 ? *ptl : *ptr;                                     \
                                                                               \
  ptl = array + 3;                                                             \
  ptr = array + 7;                                                             \
  pts = swap + 7;                                                              \
  tail_branchless_merge(pts, x, ptl, ptr, cmp);                                \
  tail_branchless_merge(pts, x, ptl, ptr, cmp);                                \
  tail_branchless_merge(pts, x, ptl, ptr, cmp);                                \
  *pts = cmp(ptl, ptr) > 0 ? *ptl : *ptr;

#define branchless_swap(pta, swap, x, cmp)                                     \
  x = 0;                                                                       \
  swap = cmp(pta, pta + 1) > 0 ? pta[x++] : pta[1];                            \
  pta[0] = pta[x];                                                             \
  pta[1] = swap;

#define swap_branchless(pta, swap, x, y, cmp)                                  \
  x = cmp(pta, pta + 1) > 0;                                                   \
  y = !x;                                                                      \
  swap = pta[y];                                                               \
  pta[0] = pta[x];                                                             \
  pta[1] = swap;

template <typename T>
using CMPFUNC = int(*)(const T *a, const T *b);
#include "quadsort.inc"


template <class T>
void quadsort(T *array, size_t nmemb, size_t size, CMPFUNC<T> cmp) noexcept {
  // can not use quadsort direct
  // if (nmemb < 2) {
  //   return;
  // }

  T *pta = array;

  if (nmemb < 32) {
    T swap[nmemb];

    (tail_swap)(pta, swap, nmemb, cmp);
  } else if (quad_swap(pta, nmemb, cmp) == 0) {
    size_t block, swap_size = nmemb;

    if (nmemb > 4194304)
      for (swap_size = 4194304; swap_size * 8 <= nmemb; swap_size *= 4) {
      }

    T swap[nmemb];

    block = quad_merge(pta, swap, swap_size, nmemb, 32, cmp);

    rotate_merge(pta, swap, swap_size, nmemb, block, cmp);
  }
}

#endif
