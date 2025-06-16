// fluxsort 1.2.1.3 - Igor van den Hoven ivdhoven@gmail.com

#define FLUX_OUT 96
template <class T>
void flux_partition(T *array, T *swap, T *ptx, T *ptp, size_t nmemb, CMPFUNC<T> cmp) noexcept;

// Determine whether to use mergesort or quicksort
template <class T>
void flux_analyze(T *array, T *swap, size_t swap_size, size_t nmemb, CMPFUNC<T> cmp) noexcept {
  unsigned char loop, asum, bsum, csum, dsum;
  unsigned int astreaks, bstreaks, cstreaks, dstreaks;
  size_t quad1, quad2, quad3, quad4, half1, half2;
  size_t cnt, abalance, bbalance, cbalance, dbalance;
  T *pta, *ptb, *ptc, *ptd;

  half1 = nmemb / 2;
  quad1 = half1 / 2;
  quad2 = half1 - quad1;
  half2 = nmemb - half1;
  quad3 = half2 / 2;
  quad4 = half2 - quad3;

  pta = array;
  ptb = array + quad1;
  ptc = array + half1;
  ptd = array + half1 + quad3;

  astreaks = bstreaks = cstreaks = dstreaks = 0;
  abalance = bbalance = cbalance = dbalance = 0;

  if (quad1 < quad2) {
    bbalance += cmp(ptb, ptb + 1) > 0;
    ptb++;
  }
  if (quad1 < quad3) {
    cbalance += cmp(ptc, ptc + 1) > 0;
    ptc++;
  }
  if (quad1 < quad4) {
    dbalance += cmp(ptd, ptd + 1) > 0;
    ptd++;
  }

  for (cnt = nmemb; cnt > 132; cnt -= 128) {
    for (asum = bsum = csum = dsum = 0, loop = 32; loop; loop--) {
      asum += cmp(pta, pta + 1) > 0;
      pta++;
      bsum += cmp(ptb, ptb + 1) > 0;
      ptb++;
      csum += cmp(ptc, ptc + 1) > 0;
      ptc++;
      dsum += cmp(ptd, ptd + 1) > 0;
      ptd++;
    }
    abalance += asum;
    astreaks += asum = (asum == 0) | (asum == 32);
    bbalance += bsum;
    bstreaks += bsum = (bsum == 0) | (bsum == 32);
    cbalance += csum;
    cstreaks += csum = (csum == 0) | (csum == 32);
    dbalance += dsum;
    dstreaks += dsum = (dsum == 0) | (dsum == 32);

    if (cnt > 516 && asum + bsum + csum + dsum == 0) {
      abalance += 48;
      pta += 96;
      bbalance += 48;
      ptb += 96;
      cbalance += 48;
      ptc += 96;
      dbalance += 48;
      ptd += 96;
      cnt -= 384;
    }
  }

  for (; cnt > 7; cnt -= 4) {
    abalance += cmp(pta, pta + 1) > 0;
    pta++;
    bbalance += cmp(ptb, ptb + 1) > 0;
    ptb++;
    cbalance += cmp(ptc, ptc + 1) > 0;
    ptc++;
    dbalance += cmp(ptd, ptd + 1) > 0;
    ptd++;
  }

  cnt = abalance + bbalance + cbalance + dbalance;

  if (cnt == 0) {
    if (cmp(pta, pta + 1) <= 0 && cmp(ptb, ptb + 1) <= 0 &&
        cmp(ptc, ptc + 1) <= 0) {
      return;
    }
  }

  asum = quad1 - abalance == 1;
  bsum = quad2 - bbalance == 1;
  csum = quad3 - cbalance == 1;
  dsum = quad4 - dbalance == 1;

  if (asum | bsum | csum | dsum) {
    unsigned char span1 = (asum && bsum) * (cmp(pta, pta + 1) > 0);
    unsigned char span2 = (bsum && csum) * (cmp(ptb, ptb + 1) > 0);
    unsigned char span3 = (csum && dsum) * (cmp(ptc, ptc + 1) > 0);

    switch (span1 | span2 * 2 | span3 * 4) {
    case 0:
      break;
    case 1:
      quad_reversal(array, ptb);
      abalance = bbalance = 0;
      break;
    case 2:
      quad_reversal(pta + 1, ptc);
      bbalance = cbalance = 0;
      break;
    case 3:
      quad_reversal(array, ptc);
      abalance = bbalance = cbalance = 0;
      break;
    case 4:
      quad_reversal(ptb + 1, ptd);
      cbalance = dbalance = 0;
      break;
    case 5:
      quad_reversal(array, ptb);
      quad_reversal(ptb + 1, ptd);
      abalance = bbalance = cbalance = dbalance = 0;
      break;
    case 6:
      quad_reversal(pta + 1, ptd);
      bbalance = cbalance = dbalance = 0;
      break;
    case 7:
      quad_reversal(array, ptd);
      return;
    }
    if (asum && abalance) {
      quad_reversal(array, pta);
      abalance = 0;
    }
    if (bsum && bbalance) {
      quad_reversal(pta + 1, ptb);
      bbalance = 0;
    }
    if (csum && cbalance) {
      quad_reversal(ptb + 1, ptc);
      cbalance = 0;
    }
    if (dsum && dbalance) {
      quad_reversal(ptc + 1, ptd);
      dbalance = 0;
    }
  }

  cnt = nmemb / 512; // switch to quadsort if at least 25% ordered

  asum = astreaks > cnt;
  bsum = bstreaks > cnt;
  csum = cstreaks > cnt;
  dsum = dstreaks > cnt;

  if (quad1 > 262144) {
    asum = bsum = csum = dsum = 1;
  }

  switch (asum + bsum * 2 + csum * 4 + dsum * 8) {
  case 0:
    flux_partition(array, swap, array, swap + nmemb, nmemb, cmp);
    return;
  case 1:
    if (abalance)
      quadsort_swap(array, swap, swap_size, quad1, cmp);
    flux_partition(pta + 1, swap, pta + 1, swap + quad2 + half2, quad2 + half2,
                   cmp);
    break;
  case 2:
    flux_partition(array, swap, array, swap + quad1, quad1, cmp);
    if (bbalance)
      quadsort_swap(pta + 1, swap, swap_size, quad2, cmp);
    flux_partition(ptb + 1, swap, ptb + 1, swap + half2, half2, cmp);
    break;
  case 3:
    if (abalance)
      quadsort_swap(array, swap, swap_size, quad1, cmp);
    if (bbalance)
      quadsort_swap(pta + 1, swap, swap_size, quad2, cmp);
    flux_partition(ptb + 1, swap, ptb + 1, swap + half2, half2, cmp);
    break;
  case 4:
    flux_partition(array, swap, array, swap + half1, half1, cmp);
    if (cbalance)
      quadsort_swap(ptb + 1, swap, swap_size, quad3, cmp);
    flux_partition(ptc + 1, swap, ptc + 1, swap + quad4, quad4, cmp);
    break;
  case 8:
    flux_partition(array, swap, array, swap + half1 + quad3, half1 + quad3,
                   cmp);
    if (dbalance)
      quadsort_swap(ptc + 1, swap, swap_size, quad4, cmp);
    break;
  case 9:
    if (abalance)
      quadsort_swap(array, swap, swap_size, quad1, cmp);
    flux_partition(pta + 1, swap, pta + 1, swap + quad2 + quad3, quad2 + quad3,
                   cmp);
    if (dbalance)
      quadsort_swap(ptc + 1, swap, swap_size, quad4, cmp);
    break;
  case 12:
    flux_partition(array, swap, array, swap + half1, half1, cmp);
    if (cbalance)
      quadsort_swap(ptb + 1, swap, swap_size, quad3, cmp);
    if (dbalance)
      quadsort_swap(ptc + 1, swap, swap_size, quad4, cmp);
    break;
  case 5:
  case 6:
  case 7:
  case 10:
  case 11:
  case 13:
  case 14:
  case 15:
    if (asum) {
      if (abalance)
        quadsort_swap(array, swap, swap_size, quad1, cmp);
    } else
      flux_partition(array, swap, array, swap + quad1, quad1, cmp);
    if (bsum) {
      if (bbalance)
        quadsort_swap(pta + 1, swap, swap_size, quad2, cmp);
    } else
      flux_partition(pta + 1, swap, pta + 1, swap + quad2, quad2, cmp);
    if (csum) {
      if (cbalance)
        quadsort_swap(ptb + 1, swap, swap_size, quad3, cmp);
    } else
      flux_partition(ptb + 1, swap, ptb + 1, swap + quad3, quad3, cmp);
    if (dsum) {
      if (dbalance)
        quadsort_swap(ptc + 1, swap, swap_size, quad4, cmp);
    } else
      flux_partition(ptc + 1, swap, ptc + 1, swap + quad4, quad4, cmp);
    break;
  }

  if (cmp(pta, pta + 1) <= 0) {
    if (cmp(ptc, ptc + 1) <= 0) {
      if (cmp(ptb, ptb + 1) <= 0) {
        return;
      }
      memcpy(swap, array, nmemb * sizeof(T));
    } else {
      cross_merge(swap + half1, array + half1, quad3, quad4, cmp);
      memcpy(swap, array, half1 * sizeof(T));
    }
  } else {
    if (cmp(ptc, ptc + 1) <= 0) {
      memcpy(swap + half1, array + half1, half2 * sizeof(T));
      cross_merge(swap, array, quad1, quad2, cmp);
    } else {
      cross_merge(swap + half1, ptb + 1, quad3, quad4, cmp);
      cross_merge(swap, array, quad1, quad2, cmp);
    }
  }
  cross_merge(array, swap, half1, half2, cmp);
}

// The next 4 functions are used for pivot selection
template <class T> T binary_median(T *pta, T *ptb, size_t len, CMPFUNC<T> cmp) noexcept {
  while (len /= 2) {
    if (cmp(pta + len, ptb + len) <= 0)
      pta += len;
    else
      ptb += len;
  }
  return cmp(pta, ptb) > 0 ? *pta : *ptb;
}

template <class T> void trim_four(T *pta, CMPFUNC<T> cmp) noexcept {
  T swap;
  size_t x;

  x = cmp(pta, pta + 1) > 0;
  swap = pta[!x];
  pta[0] = pta[x];
  pta[1] = swap;
  pta += 2;
  x = cmp(pta, pta + 1) > 0;
  swap = pta[!x];
  pta[0] = pta[x];
  pta[1] = swap;
  pta -= 2;

  x = (cmp(pta, pta + 2) <= 0) * 2;
  pta[2] = pta[x];
  pta++;
  x = (cmp(pta, pta + 2) > 0) * 2;
  pta[0] = pta[x];
}

template <class T> T median_of_nine(T *array, size_t nmemb, CMPFUNC<T> cmp) noexcept {
  T *pta, swap[9];
  size_t x, y, z;

  z = nmemb / 9;

  pta = array;

  for (x = 0; x < 9; x++) {
    swap[x] = *pta;

    pta += z;
  }

  trim_four(swap, cmp);
  trim_four(swap + 4, cmp);

  swap[0] = swap[5];
  swap[3] = swap[8];

  trim_four(swap, cmp);

  swap[0] = swap[6];

  x = cmp(swap + 0, swap + 1) > 0;
  y = cmp(swap + 0, swap + 2) > 0;
  z = cmp(swap + 1, swap + 2) > 0;

  return swap[(x == y) + (y ^ z)];
}

template <class T>
T median_of_cbrt(T *array, T *swap, T *ptx, size_t nmemb, int *generic, CMPFUNC<T> cmp) noexcept {
  T *pta, *pts;
  size_t cnt, div, cbrt;

  for (cbrt = 32; nmemb > cbrt * cbrt * cbrt; cbrt *= 2) {
  }

  div = nmemb / cbrt;

  pta = ptx + (size_t)&div / 16 % div;
  pts = ptx == array ? swap : array;

  for (cnt = 0; cnt < cbrt; cnt++) {
    pts[cnt] = *pta;

    pta += div;
  }
  cbrt /= 2;

  quadsort_swap(pts, pts + cbrt * 2, cbrt, cbrt, cmp);
  quadsort_swap(pts + cbrt, pts + cbrt * 2, cbrt, cbrt, cmp);

  *generic =
      (cmp(pts + cbrt * 2 - 1, pts) <= 0) & (cmp(pts + cbrt - 1, pts) <= 0);

  return binary_median(pts, pts + cbrt, cbrt, cmp);
}

// As per suggestion by Marshall Lochbaum to improve generic data handling by
// mimicking dual-pivot quicksort
template <class T>
void flux_reverse_partition(T *array, T *swap, T *ptx, T *piv, size_t nmemb,
                            CMPFUNC<T> cmp) noexcept{
  size_t a_size, s_size;

  {
    size_t cnt;
    T *tmp, *pta = array, *pts = swap;

    for (cnt = nmemb / 8; cnt; cnt--) {
      tmp = cmp(piv, ptx) > 0 ? pta++ : pts++;
      *tmp = *ptx++;
      tmp = cmp(piv, ptx) > 0 ? pta++ : pts++;
      *tmp = *ptx++;
      tmp = cmp(piv, ptx) > 0 ? pta++ : pts++;
      *tmp = *ptx++;
      tmp = cmp(piv, ptx) > 0 ? pta++ : pts++;
      *tmp = *ptx++;
      tmp = cmp(piv, ptx) > 0 ? pta++ : pts++;
      *tmp = *ptx++;
      tmp = cmp(piv, ptx) > 0 ? pta++ : pts++;
      *tmp = *ptx++;
      tmp = cmp(piv, ptx) > 0 ? pta++ : pts++;
      *tmp = *ptx++;
      tmp = cmp(piv, ptx) > 0 ? pta++ : pts++;
      *tmp = *ptx++;
    }

    for (cnt = nmemb % 8; cnt; cnt--) {
      tmp = cmp(piv, ptx) > 0 ? pta++ : pts++;
      *tmp = *ptx++;
    }
    a_size = pta - array;
    s_size = pts - swap;
  }

  memcpy(array + a_size, swap, s_size * sizeof(T));

  if (s_size <= a_size / 16 || a_size <= FLUX_OUT) {
    quadsort_swap(array, swap, a_size, a_size, cmp);
    return;
  }
  flux_partition(array, swap, array, piv, a_size, cmp);
}

template <class T>
size_t flux_default_partition(T *array, T *swap, T *ptx, T *piv, size_t nmemb,
                              CMPFUNC<T> cmp) noexcept{
  size_t run = 0, a = 0, m = 0;

  T *tmp, *pta = array, *pts = swap;

  for (a = 8; a <= nmemb; a += 8) {
    tmp = cmp(ptx, piv) <= 0 ? pta++ : pts++;
    *tmp = *ptx++;
    tmp = cmp(ptx, piv) <= 0 ? pta++ : pts++;
    *tmp = *ptx++;
    tmp = cmp(ptx, piv) <= 0 ? pta++ : pts++;
    *tmp = *ptx++;
    tmp = cmp(ptx, piv) <= 0 ? pta++ : pts++;
    *tmp = *ptx++;

    tmp = cmp(ptx, piv) <= 0 ? pta++ : pts++;
    *tmp = *ptx++;
    tmp = cmp(ptx, piv) <= 0 ? pta++ : pts++;
    *tmp = *ptx++;
    tmp = cmp(ptx, piv) <= 0 ? pta++ : pts++;
    *tmp = *ptx++;
    tmp = cmp(ptx, piv) <= 0 ? pta++ : pts++;
    *tmp = *ptx++;

    if (pta == array || pts == swap)
      run = a;
  }

  for (a = nmemb % 8; a; a--) {
    tmp = cmp(ptx, piv) <= 0 ? pta++ : pts++;
    *tmp = *ptx++;
  }
  m = pta - array;

  if (run <= nmemb / 4) {
    return m;
  }

  if (m == nmemb) {
    return m;
  }

  a = nmemb - m;

  memcpy(array + m, swap, a * sizeof(T));

  quadsort_swap(array + m, swap, a, a, cmp);
  quadsort_swap(array, swap, m, m, cmp);

  return 0;
}

template <class T>
void flux_partition(T *array, T *swap, T *ptx, T *piv, size_t nmemb, CMPFUNC<T> cmp) noexcept{
  size_t a_size = 0, s_size;
  int generic = 0;

  while (1) {
    --piv;

    if (nmemb <= 2048) {
      *piv = median_of_nine(ptx, nmemb, cmp);
    } else {
      *piv = median_of_cbrt(array, swap, ptx, nmemb, &generic, cmp);

      if (generic) {
        if (ptx == swap) {
          memcpy(array, swap, nmemb * sizeof(T));
        }
        quadsort_swap(array, swap, nmemb, nmemb, cmp);
        return;
      }
    }

    if (a_size && cmp(piv + 1, piv) <= 0) {
      flux_reverse_partition(array, swap, array, piv, nmemb, cmp);
      return;
    }
    a_size = flux_default_partition(array, swap, ptx, piv, nmemb, cmp);
    s_size = nmemb - a_size;

    if (a_size <= s_size / 32 || s_size <= FLUX_OUT) {
      if (a_size == 0) {
        return;
      }
      if (s_size == 0) {
        flux_reverse_partition(array, swap, array, piv, a_size, cmp);
        return;
      }
      memcpy(array + a_size, swap, s_size * sizeof(T));
      quadsort_swap(array + a_size, swap, s_size, s_size, cmp);
    } else {
      flux_partition(array + a_size, swap, swap, piv, s_size, cmp);
    }

    if (s_size <= a_size / 32 || a_size <= FLUX_OUT) {
      if (a_size <= FLUX_OUT) {
        quadsort_swap(array, swap, a_size, a_size, cmp);
      } else {
        flux_reverse_partition(array, swap, array, piv, a_size, cmp);
      }
      return;
    }
    nmemb = a_size;
    ptx = array;
  }
}
