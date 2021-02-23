/*
**                    TComb v2.x for Avisynth 2.6 and Avisynth+
**
**   TComb is a temporal comb filter (it reduces cross-luminance (rainbowing)
**   and cross-chrominance (dot crawl) artifacts in static areas of the picture).
**   It will ONLY work with NTSC material, and WILL NOT work with telecined material
**   where the rainbowing/dotcrawl was introduced prior to the telecine process!
**   It must be used before ivtc or deinterlace.
**
**   Copyright (C) 2021 Ferenc Pintér
**
**   Copyright (C) 2015 Shane Panke
**
**   Copyright (C) 2005-2006 Kevin Stone
**
**   This program is free software; you can redistribute it and/or modify
**   it under the terms of the GNU General Public License as published by
**   the Free Software Foundation; either version 2 of the License, or
**   (at your option) any later version.
**
**   This program is distributed in the hope that it will be useful,
**   but WITHOUT ANY WARRANTY; without even the implied warranty of
**   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**   GNU General Public License for more details.
**
**   You should have received a copy of the GNU General Public License
**   along with this program; if not, write to the Free Software
**   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include "TComb.h"
#include <stdint.h>

#ifdef INTEL_INTRINSICS
#include <xmmintrin.h>
#include <emmintrin.h>
#endif
#include <algorithm>

template<typename pixel_t>
void checkSceneChangePlanar_1_c(const pixel_t* srcp, const pixel_t* nxtp,
  int height, int width, int src_pitch, int nxt_pitch, uint64_t& diff)
{
  for (int y = 0; y < height; ++y)
  {
    uint32_t rowdiff = 0;
    for (int x = 0; x < width; x += 4)
    {
      rowdiff += abs(srcp[x + 0] - nxtp[x + 0]);
      rowdiff += abs(srcp[x + 1] - nxtp[x + 1]);
      rowdiff += abs(srcp[x + 2] - nxtp[x + 2]);
      rowdiff += abs(srcp[x + 3] - nxtp[x + 3]);
    }
    diff += rowdiff;
    srcp += src_pitch;
    nxtp += nxt_pitch;
  }
}

// instantiate
template void checkSceneChangePlanar_1_c<uint8_t>(const uint8_t* srcp, const uint8_t* nxtp,
  int height, int width, int src_pitch, int nxt_pitch, uint64_t& diff);

#ifdef INTEL_INTRINSICS
void checkSceneChangePlanar_1_SSE2_simd(const uint8_t* prvp, const uint8_t* srcp,
  int height, int width, int prv_pitch, int src_pitch, uint64_t& diffp)
{
  __m128i sum = _mm_setzero_si128();
  while (height--) {
    for (int x = 0; x < width; x += 16)
    {
      __m128i src1 = _mm_load_si128(reinterpret_cast<const __m128i*>(prvp + x));
      __m128i src2 = _mm_load_si128(reinterpret_cast<const __m128i*>(srcp + x));
      __m128i sad = _mm_sad_epu8(src1, src2);
      sum = _mm_add_epi32(sum, sad);
    }
    prvp += prv_pitch;
    srcp += src_pitch;
  }
  __m128i res = _mm_add_epi32(sum, _mm_srli_si128(sum, 8));
  diffp = _mm_cvtsi128_si32(res);
}
#endif

#ifdef INTEL_INTRINSICS
void andMasks_SSE2_simd(const uint8_t* s1p, const uint8_t* s2p, uint8_t* dstp, int stride, int width, int height)
{
  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; x += 16)
    {
      __m128i src1 = _mm_load_si128(reinterpret_cast<const __m128i*>(s1p + x));
      __m128i src2 = _mm_load_si128(reinterpret_cast<const __m128i*>(s2p + x));
      __m128i result = _mm_and_si128(src1, src2);
      _mm_store_si128(reinterpret_cast<__m128i*>(dstp + x), result);
    }

    s1p += stride;
    s2p += stride;
    dstp += stride;
  }
}
#endif

void andMasks_c(const uint8_t* s1p, const uint8_t* s2p, uint8_t* dstp, int stride, int width, int height)
{
  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; ++x)
      dstp[x] = (s1p[x] & s2p[x]);

    s1p += stride;
    s2p += stride;
    dstp += stride;
  }
}

#ifdef INTEL_INTRINSICS
void orAndMasks_SSE2_simd(const uint8_t* s1p, const uint8_t* s2p, uint8_t* dstp, int stride, int width, int height)
{
  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; x += 16)
    {
      __m128i src1 = _mm_load_si128(reinterpret_cast<const __m128i*>(s1p + x));
      __m128i src2 = _mm_load_si128(reinterpret_cast<const __m128i*>(s2p + x));
      __m128i dst = _mm_load_si128(reinterpret_cast<const __m128i*>(dstp + x));
      __m128i result = _mm_or_si128(dst, _mm_and_si128(src1, src2));
      _mm_store_si128(reinterpret_cast<__m128i*>(dstp + x), result);
    }

    s1p += stride;
    s2p += stride;
    dstp += stride;
  }
}
#endif

void orAndMasks_c(const uint8_t* s1p, const uint8_t* s2p, uint8_t* dstp, int stride, int width, int height)
{
  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; ++x)
      dstp[x] |= (s1p[x] & s2p[x]);

    s1p += stride;
    s2p += stride;
    dstp += stride;
  }
}

#ifdef INTEL_INTRINSICS
void or3Masks_SSE2_simd(const uint8_t* s1p, const uint8_t* s2p, const uint8_t* s3p, uint8_t* dstp, int stride, int width, int height)
{
  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; x += 16)
    {
      __m128i src1 = _mm_load_si128(reinterpret_cast<const __m128i*>(s1p + x));
      __m128i src2 = _mm_load_si128(reinterpret_cast<const __m128i*>(s2p + x));
      __m128i src3 = _mm_load_si128(reinterpret_cast<const __m128i*>(s3p + x));
      __m128i result = _mm_or_si128(src1, _mm_or_si128(src2, src3));
      _mm_store_si128(reinterpret_cast<__m128i*>(dstp + x), result);
    }

    s1p += stride;
    s2p += stride;
    s3p += stride;
    dstp += stride;
  }
}
#endif

void or3Masks_c(const uint8_t* s1p, const uint8_t* s2p, const uint8_t* s3p, uint8_t* dstp, int stride, int width, int height)
{
  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; ++x)
      dstp[x] = (s1p[x] | s2p[x] | s3p[x]);

    s1p += stride;
    s2p += stride;
    s3p += stride;
    dstp += stride;
  }
}

#ifdef INTEL_INTRINSICS
void calcAverages_SSE2_simd(const uint8_t* s1p, const uint8_t* s2p, uint8_t* dstp, int stride, int width, int height)
{
  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; x += 16)
    {
      __m128i src1 = _mm_load_si128(reinterpret_cast<const __m128i*>(s1p + x));
      __m128i src2 = _mm_load_si128(reinterpret_cast<const __m128i*>(s2p + x));
      __m128i result = _mm_avg_epu8(src1, src2);
      _mm_store_si128(reinterpret_cast<__m128i*>(dstp + x), result);
    }

    s1p += stride;
    s2p += stride;
    dstp += stride;
  }
}
#endif

void calcAverages_c(const uint8_t* s1p, const uint8_t* s2p, uint8_t* dstp, int stride, int width, int height)
{
  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; ++x)
      dstp[x] = (s1p[x] + s2p[x] + 1) >> 1;

    s1p += stride;
    s2p += stride;
    dstp += stride;
  }
}

#ifdef INTEL_INTRINSICS
void MinMax_SSE2_simd(const uint8_t* srcp, uint8_t* dstpMin, uint8_t* dstpMax, int src_stride, int dmin_stride, int width, int height, int thresh)
{
  const uint8_t* srcpp = srcp - src_stride;
  const uint8_t* srcpn = srcp + src_stride;

  const auto threshp = _mm_set1_epi8(thresh);

  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; x += 16)
    {
      __m128i srcpp_m_1 = _mm_load_si128(reinterpret_cast<const __m128i*>(srcpp + x - 1));
      __m128i srcpp_0 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(srcpp + x));
      __m128i srcpp_p_1 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(srcpp + x + 1));

      __m128i srcp_m_1 = _mm_load_si128(reinterpret_cast<const __m128i*>(srcp + x - 1));
      __m128i srcp_0 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(srcp + x));
      __m128i srcp_p_1 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(srcp + x + 1));

      __m128i srcpn_m_1 = _mm_load_si128(reinterpret_cast<const __m128i*>(srcpn + x - 1));
      __m128i srcpn_0 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(srcpn + x));
      __m128i srcpn_p_1 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(srcpn + x + 1));

      auto tmpmin = _mm_min_epu8(_mm_min_epu8(_mm_min_epu8(_mm_min_epu8(srcpp_m_1, srcpp_0),
        _mm_min_epu8(srcpp_p_1, srcp_m_1)),
        _mm_min_epu8(_mm_min_epu8(srcp_0, srcp_p_1),
          _mm_min_epu8(srcpn_m_1, srcpn_0))), srcpn_p_1);

      auto min = _mm_subs_epu8(tmpmin, threshp);

      _mm_store_si128(reinterpret_cast<__m128i*>(dstpMin + x), min);

      auto tmpmax = _mm_max_epu8(_mm_max_epu8(_mm_max_epu8(_mm_max_epu8(srcpp_m_1, srcpp_0),
        _mm_max_epu8(srcpp_p_1, srcp_m_1)),
        _mm_max_epu8(_mm_max_epu8(srcp_0, srcp_p_1),
          _mm_max_epu8(srcpn_m_1, srcpn_0))), srcpn_p_1);

      auto max = _mm_adds_epu8(tmpmax, threshp); // future warning: 10-14 bitss

      _mm_store_si128(reinterpret_cast<__m128i*>(dstpMax + x), max);
    }

    srcpp += src_stride;
    srcp += src_stride;
    srcpn += src_stride;
    dstpMin += dmin_stride;
    dstpMax += dmin_stride;
  }
}
#endif

void MinMax_c(const uint8_t* srcp, uint8_t* dstpMin, uint8_t* dstpMax, int src_stride, int dmin_stride, int width, int height, int thresh)
{
  const uint8_t* srcpp = srcp - src_stride;
  const uint8_t* srcpn = srcp + src_stride;

  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; ++x)
    {
      dstpMin[x] = std::max(std::min(std::min(std::min(std::min(srcpp[x - 1], srcpp[x]),
        std::min(srcpp[x + 1], srcp[x - 1])),
        std::min(std::min(srcp[x], srcp[x + 1]),
          std::min(srcpn[x - 1], srcpn[x]))), srcpn[x + 1]) - thresh, 0);
      dstpMax[x] = std::min(std::max(std::max(std::max(std::max(srcpp[x - 1], srcpp[x]),
        std::max(srcpp[x + 1], srcp[x - 1])),
        std::max(std::max(srcp[x], srcp[x + 1]),
          std::max(srcpn[x - 1], srcpn[x]))), srcpn[x + 1]) + thresh, 255);
    }

    srcpp += src_stride;
    srcp += src_stride;
    srcpn += src_stride;
    dstpMin += dmin_stride;
    dstpMax += dmin_stride;
  }
}

#ifdef INTEL_INTRINSICS
void absDiff_SSE2_simd(const uint8_t* srcp1, const uint8_t* srcp2, uint8_t* dstp, int stride, int width, int height)
{
  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; x += 16) {
      auto src1 = _mm_load_si128(reinterpret_cast<const __m128i*>(srcp1 + x));
      auto src2 = _mm_load_si128(reinterpret_cast<const __m128i*>(srcp2 + x));
      auto diff12 = _mm_subs_epu8(src1, src2);
      auto diff21 = _mm_subs_epu8(src2, src1);
      auto diff = _mm_or_si128(diff12, diff21);
      _mm_store_si128(reinterpret_cast<__m128i*>(dstp + x), diff);
    }

    srcp1 += stride;
    srcp2 += stride;
    dstp += stride;
  }
}
#endif

void absDiff_c(const uint8_t* srcp1, const uint8_t* srcp2, uint8_t* dstp, int stride, int width, int height)
{
  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; ++x)
      dstp[x] = abs(srcp1[x] - srcp2[x]);

    srcp1 += stride;
    srcp2 += stride;
    dstp += stride;
  }
}

#ifdef INTEL_INTRINSICS
void buildFinalMask_SSE2_simd(const uint8_t* s1p, const uint8_t* s2p, const uint8_t* m1p, uint8_t* dstp, int stride, int width, int height, int thresh)
{
  auto thresh_minus1 = _mm_set1_epi8(thresh-1);
  auto zero = _mm_setzero_si128();

  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; x += 16)
    {
      auto src1 = _mm_load_si128(reinterpret_cast<const __m128i*>(s1p + x));
      auto src2 = _mm_load_si128(reinterpret_cast<const __m128i*>(s2p + x));
      auto diff12 = _mm_subs_epu8(src1, src2);
      auto diff21 = _mm_subs_epu8(src2, src1);
      auto diff = _mm_or_si128(diff12, diff21);
      auto addedsthresh = _mm_subs_epu8(diff, thresh_minus1);
      auto cmpresult = _mm_cmpeq_epi8(addedsthresh, zero);
      auto m1 = _mm_load_si128(reinterpret_cast<const __m128i*>(m1p + x));
      auto tmp = _mm_and_si128(cmpresult, m1);
      _mm_store_si128(reinterpret_cast<__m128i*>(dstp + x), tmp);

      /*
      if (m1p[x] && abs(s1p[x] - s2p[x]) < thresh)
        dstp[x] = 0xFF;
      else
        dstp[x] = 0;
      */
    }

    m1p += stride;
    s1p += stride;
    s2p += stride;
    dstp += stride;
  }
}
#endif

void buildFinalMask_c(const uint8_t* s1p, const uint8_t* s2p, const uint8_t* m1p, uint8_t* dstp, int stride, int width, int height, int thresh)
{
  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; ++x)
    {
      if (m1p[x] && abs(s1p[x] - s2p[x]) < thresh)
        dstp[x] = 0xFF;
      else
        dstp[x] = 0;
    }

    m1p += stride;
    s1p += stride;
    s2p += stride;
    dstp += stride;
  }
}

#ifdef INTEL_INTRINSICS
void checkOscillation5_SSE2_simd(const uint8_t* p2p, const uint8_t* p1p, const uint8_t* s1p, const uint8_t* n1p, const uint8_t* n2p, uint8_t* dstp, int stride, int width, int height, int thresh)
{
  int threshm1 = std::min(std::max(thresh - 1, 0), 255);
  auto thresh_minus1 = _mm_set1_epi8(threshm1);
  auto one = _mm_set1_epi8(1);
  auto zero = _mm_setzero_si128();

  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; x += 16)
    {
      // trick: x < thresh ==> x <= (thresh - 1) ==> x - (thresh - 1) <= 0 ==> sub_sat(x, thresh - 1) == 0
      // pcmpeqb(psubusb(x, thresh - 1), zero) : 0xFF where x < thresh

      __m128i src_p2p = _mm_load_si128(reinterpret_cast<const __m128i*>(p2p + x));
      __m128i src_s1p = _mm_load_si128(reinterpret_cast<const __m128i*>(s1p + x));
      __m128i src_n2p = _mm_load_si128(reinterpret_cast<const __m128i*>(n2p + x));
      __m128i src_p1p = _mm_load_si128(reinterpret_cast<const __m128i*>(p1p + x));
      __m128i src_n1p = _mm_load_si128(reinterpret_cast<const __m128i*>(n1p + x));

      auto min31 = _mm_min_epu8(_mm_min_epu8(src_p2p, src_s1p), src_n2p);
      auto max31 = _mm_max_epu8(_mm_max_epu8(src_p2p, src_s1p), src_n2p);
      auto min22 = _mm_min_epu8(src_p1p, src_n1p);
      auto max22 = _mm_max_epu8(src_p1p, src_n1p);

      auto cmp1 = _mm_cmpeq_epi8(_mm_subs_epu8(max22, _mm_subs_epu8(min31, one)), zero);
      auto cmp2 = _mm_cmpeq_epi8(_mm_subs_epu8(max31, _mm_subs_epu8(min22, one)), zero);
      // No check for (max22 == 0) or (max31 == 0), like in C, sub_sat handles automatically
      auto maxmindiff31 = _mm_subs_epu8(max31, min31);
      auto cmp3 = _mm_cmpeq_epi8(_mm_subs_epu8(maxmindiff31, thresh_minus1), zero);
      auto maxmindiff22 = _mm_subs_epu8(max22, min22);
      auto cmp4 = _mm_cmpeq_epi8(_mm_subs_epu8(maxmindiff22, thresh_minus1), zero);

      auto result = _mm_and_si128(_mm_or_si128(cmp1, cmp2), _mm_and_si128(cmp3, cmp4));
      /*
      if (((max22 < min31) || max22 == 0 || (max31 < min22) || max31 == 0) &&
        max31 - min31 < thresh && max22 - min22 < thresh)
        dstp[x] = 0xFF;
      else dstp[x] = 0;
      */
      _mm_store_si128(reinterpret_cast<__m128i*>(dstp + x), result);
    }

    p2p += stride;
    p1p += stride;
    s1p += stride;
    n1p += stride;
    n2p += stride;
    dstp += stride;
  }
}
#endif

void checkOscillation5_c(const uint8_t* p2p, const uint8_t* p1p, const uint8_t* s1p, const uint8_t* n1p, const uint8_t* n2p, uint8_t* dstp, int stride, int width, int height, int thresh)
{
  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; ++x)
    {
      const int min31 = min3(p2p[x], s1p[x], n2p[x]);
      const int max31 = max3(p2p[x], s1p[x], n2p[x]);
      const int min22 = std::min(p1p[x], n1p[x]);
      const int max22 = std::max(p1p[x], n1p[x]);
      if (((max22 < min31) || max22 == 0 || (max31 < min22) || max31 == 0) &&
        max31 - min31 < thresh && max22 - min22 < thresh)
        dstp[x] = 0xFF;
      else dstp[x] = 0;
    }

    p2p += stride;
    p1p += stride;
    s1p += stride;
    n1p += stride;
    n2p += stride;
    dstp += stride;
  }
}

#ifdef INTEL_INTRINSICS
void absDiffAndMinMaskThresh_SSE2_simd(const uint8_t* srcp1, const uint8_t* srcp2, uint8_t* dstp, int stride, int width, int height, int thresh)
{
  int threshm1 = std::min(std::max(thresh - 1, 0), 255);
  auto thresh_minus1 = _mm_set1_epi8(threshm1);
  auto zero = _mm_setzero_si128();

  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; x += 16)
    {
      __m128i src1 = _mm_load_si128(reinterpret_cast<const __m128i*>(srcp1 + x));
      __m128i src2 = _mm_load_si128(reinterpret_cast<const __m128i*>(srcp2 + x));
      __m128i dst = _mm_load_si128(reinterpret_cast<const __m128i*>(dstp + x));
      auto diff12 = _mm_subs_epu8(src1, src2);
      auto diff21 = _mm_subs_epu8(src2, src1);
      auto diff = _mm_or_si128(diff12, diff21);

      auto tmp_min = _mm_min_epu8(diff, dst);
      auto result = _mm_cmpeq_epi8(_mm_subs_epu8(tmp_min, thresh_minus1), zero);
      /*
      if (diff < dstp[x]) dstp[x] = diff; // min
      if (dstp[x] < thresh)
        dstp[x] = 0xFF;
      else
        dstp[x] = 0;
      */

      _mm_store_si128(reinterpret_cast<__m128i*>(dstp + x), result);

    }

    srcp1 += stride;
    srcp2 += stride;
    dstp += stride;
  }
}
#endif

void absDiffAndMinMaskThresh_c(const uint8_t* srcp1, const uint8_t* srcp2, uint8_t* dstp, int stride, int width, int height, int thresh)
{
  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; ++x)
    {
      const int diff = abs(srcp1[x] - srcp2[x]);
      if (diff < dstp[x])
        dstp[x] = diff;
      if (dstp[x] < thresh)
        dstp[x] = 0xFF;
      else
        dstp[x] = 0;
    }

    srcp1 += stride;
    srcp2 += stride;
    dstp += stride;
  }
}

#ifdef INTEL_INTRINSICS
void absDiffAndMinMask_SSE2_simd(const uint8_t* srcp1, const uint8_t* srcp2, uint8_t* dstp, int stride, int width, int height)
{
  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; x += 16)
    {
      __m128i src1 = _mm_load_si128(reinterpret_cast<const __m128i*>(srcp1 + x));
      __m128i src2 = _mm_load_si128(reinterpret_cast<const __m128i*>(srcp2 + x));
      __m128i dst = _mm_load_si128(reinterpret_cast<const __m128i*>(dstp + x));
      auto diff12 = _mm_subs_epu8(src1, src2);
      auto diff21 = _mm_subs_epu8(src2, src1);
      auto diff = _mm_or_si128(diff12, diff21);

      auto tmp_min = _mm_min_epu8(diff, dst);
      _mm_store_si128(reinterpret_cast<__m128i*>(dstp + x), tmp_min);

      /*
      const int diff = abs(srcp1[x] - srcp2[x]);
      if (diff < dstp[x])
        dstp[x] = diff;
      */
    }

    srcp1 += stride;
    srcp2 += stride;
    dstp += stride;
  }
}
#endif

void absDiffAndMinMask_c(const uint8_t* srcp1, const uint8_t* srcp2, uint8_t* dstp, int stride, int width, int height)
{
  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; ++x)
    {
      const int diff = abs(srcp1[x] - srcp2[x]);
      if (diff < dstp[x])
        dstp[x] = diff;
    }

    srcp1 += stride;
    srcp2 += stride;
    dstp += stride;
  }
}

#ifdef INTEL_INTRINSICS
void checkAvgOscCorrelation_SSE2_simd(const uint8_t* s1p, const uint8_t* s2p, const uint8_t* s3p, const uint8_t* s4p, uint8_t* dstp, int stride, int width, int height, int thresh)
{
  int threshm1 = std::min(std::max(thresh - 1, 0), 255);
  auto thresh_minus1 = _mm_set1_epi8(threshm1);
  auto zero = _mm_setzero_si128();

  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; x += 16)
    {
      __m128i s1 = _mm_load_si128(reinterpret_cast<const __m128i*>(s1p + x));
      __m128i s2 = _mm_load_si128(reinterpret_cast<const __m128i*>(s2p + x));
      __m128i s3 = _mm_load_si128(reinterpret_cast<const __m128i*>(s3p + x));
      __m128i s4 = _mm_load_si128(reinterpret_cast<const __m128i*>(s4p + x));

      auto min = _mm_min_epu8(_mm_min_epu8(_mm_min_epu8(s1, s2), s3), s4);
      auto max = _mm_max_epu8(_mm_max_epu8(_mm_max_epu8(s1, s2), s3), s4);

      auto diffmaxmin = _mm_subs_epu8(max, min);
      auto cmp = _mm_cmpeq_epi8(_mm_subs_epu8(diffmaxmin, thresh_minus1), zero);

      __m128i dst = _mm_load_si128(reinterpret_cast<const __m128i*>(dstp + x));
      auto result = _mm_and_si128(cmp, dst);
      _mm_store_si128(reinterpret_cast<__m128i*>(dstp + x), result);

      /*
      if (max4(s1p[x], s2p[x], s3p[x], s4p[x]) - min4(s1p[x], s2p[x], s3p[x], s4p[x]) >= thresh)
        dstp[x] = 0;
      that is: 
      if(max-min < thresh) dstp[x] = dstp[x] else 0    (dst=dst&FF   0=dst&00)
      */
    }

    s1p += stride;
    s2p += stride;
    s3p += stride;
    s4p += stride;
    dstp += stride;
  }
}
#endif

void checkAvgOscCorrelation_c(const uint8_t* s1p, const uint8_t* s2p, const uint8_t* s3p, const uint8_t* s4p, uint8_t* dstp, int stride, int width, int height, int thresh)
{
  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; ++x)
    {
      if (max4(s1p[x], s2p[x], s3p[x], s4p[x]) -
        min4(s1p[x], s2p[x], s3p[x], s4p[x]) >= thresh)
        dstp[x] = 0;
    }

    s1p += stride;
    s2p += stride;
    s3p += stride;
    s4p += stride;
    dstp += stride;
  }
}

#ifdef INTEL_INTRINSICS
void VerticalBlur3_SSE2_simd(const uint8_t* srcp, uint8_t* dstp, int stride, int width, int height)
{
  const uint8_t* srcpp = srcp - stride;
  const uint8_t* srcpn = srcp + stride;

  auto zero = _mm_setzero_si128();
  auto two = _mm_set1_epi16(2);

  // top line
  for (int x = 0; x < width; x += 16) {
    __m128i s1 = _mm_load_si128(reinterpret_cast<const __m128i*>(srcp + x));
    __m128i s2 = _mm_load_si128(reinterpret_cast<const __m128i*>(srcpn + x));
    auto avg = _mm_avg_epu8(s1, s2);
    _mm_store_si128(reinterpret_cast<__m128i*>(dstp + x), avg);
    // dstp[x] = (srcp[x] + srcpn[x] + 1) >> 1;
  }

  srcpp += stride;
  srcp += stride;
  srcpn += stride;
  dstp += stride;

  for (int y = 1; y < height - 1; ++y)
  {
    for (int x = 0; x < width; x += 16) {
      __m128i p = _mm_load_si128(reinterpret_cast<const __m128i*>(srcpp + x));
      __m128i s = _mm_load_si128(reinterpret_cast<const __m128i*>(srcp + x));
      __m128i n = _mm_load_si128(reinterpret_cast<const __m128i*>(srcpn + x));

      auto p_lo = _mm_unpacklo_epi8(p, zero);
      auto p_hi = _mm_unpackhi_epi8(p, zero);
      auto s_lo = _mm_unpacklo_epi8(s, zero);
      auto s_hi = _mm_unpackhi_epi8(s, zero);
      auto n_lo = _mm_unpacklo_epi8(n, zero);
      auto n_hi = _mm_unpackhi_epi8(n, zero);
      auto res_lo = _mm_add_epi16(_mm_add_epi16(p_lo, _mm_slli_epi16(s_lo, 1)), n_lo);
      auto res_hi = _mm_add_epi16(_mm_add_epi16(p_hi, _mm_slli_epi16(s_hi, 1)), n_hi);
      res_lo = _mm_srli_epi16(_mm_add_epi16(res_lo, two), 2);
      res_hi = _mm_srli_epi16(_mm_add_epi16(res_hi, two), 2);
      auto result = _mm_packus_epi16(res_lo, res_hi);
      _mm_store_si128(reinterpret_cast<__m128i*>(dstp + x), result);
      // dstp[x] = (srcpp[x] + (srcp[x] << 1) + srcpn[x] + 2) >> 2;
    }

    srcpp += stride;
    srcp += stride;
    srcpn += stride;
    dstp += stride;
  }

  // bottom
  for (int x = 0; x < width; x += 16) {
    __m128i s1 = _mm_load_si128(reinterpret_cast<const __m128i*>(srcpp + x));
    __m128i s2 = _mm_load_si128(reinterpret_cast<const __m128i*>(srcp + x));
    auto avg = _mm_avg_epu8(s1, s2);
    _mm_store_si128(reinterpret_cast<__m128i*>(dstp + x), avg);
    //dstp[x] = (srcpp[x] + srcp[x] + 1) >> 1;
  }

}
#endif

void VerticalBlur3_c(const uint8_t* srcp, uint8_t* dstp, int stride, int width, int height)
{
  const uint8_t* srcpp = srcp - stride;
  const uint8_t* srcpn = srcp + stride;

  for (int x = 0; x < width; ++x)
    dstp[x] = (srcp[x] + srcpn[x] + 1) >> 1;

  srcpp += stride;
  srcp += stride;
  srcpn += stride;
  dstp += stride;

  for (int y = 1; y < height - 1; ++y)
  {
    for (int x = 0; x < width; ++x)
      dstp[x] = (srcpp[x] + (srcp[x] << 1) + srcpn[x] + 2) >> 2;

    srcpp += stride;
    srcp += stride;
    srcpn += stride;
    dstp += stride;
  }

  for (int x = 0; x < width; ++x)
    dstp[x] = (srcpp[x] + srcp[x] + 1) >> 1;
}

#ifdef INTEL_INTRINSICS
// width mod 16 and srcp alignment guaranteed
void HorizontalBlur3_SSE2_simd(const uint8_t* srcp, uint8_t* dstp, int stride, int width, int height)
{
  auto zero = _mm_setzero_si128();
  auto two = _mm_set1_epi16(2);

  for (int y = 0; y < height; ++y)
  {
    for (int x = 0; x < width; x += 16)
    {
      __m128i p = _mm_loadu_si128(reinterpret_cast<const __m128i*>(srcp + x - 1));
      __m128i s = _mm_load_si128(reinterpret_cast<const __m128i*>(srcp + x));
      __m128i n = _mm_loadu_si128(reinterpret_cast<const __m128i*>(srcp + x + 1));

      auto p_lo = _mm_unpacklo_epi8(p, zero);
      auto p_hi = _mm_unpackhi_epi8(p, zero);
      auto s_lo = _mm_unpacklo_epi8(s, zero);
      auto s_hi = _mm_unpackhi_epi8(s, zero);
      auto n_lo = _mm_unpacklo_epi8(n, zero);
      auto n_hi = _mm_unpackhi_epi8(n, zero);
      auto res_lo = _mm_add_epi16(_mm_add_epi16(p_lo, _mm_slli_epi16(s_lo, 1)), n_lo);
      auto res_hi = _mm_add_epi16(_mm_add_epi16(p_hi, _mm_slli_epi16(s_hi, 1)), n_hi);
      res_lo = _mm_srli_epi16(_mm_add_epi16(res_lo, two), 2);
      res_hi = _mm_srli_epi16(_mm_add_epi16(res_hi, two), 2);
      auto result = _mm_packus_epi16(res_lo, res_hi);
      _mm_store_si128(reinterpret_cast<__m128i*>(dstp + x), result);
      // dstp[x] = (srcp[x - 1] + (srcp[x] << 1) + srcp[x + 1] + 2) >> 2;
    }

    srcp += stride;
    dstp += stride;
  }

}
#endif

void HorizontalBlur3_c(const uint8_t* srcp, uint8_t* dstp, int stride, int width, int height)
{
  for (int y = 0; y < height; ++y)
  {
    dstp[0] = (srcp[0] + srcp[1] + 1) >> 1;

    for (int x = 1; x < width - 1; ++x)
      dstp[x] = (srcp[x - 1] + (srcp[x] << 1) + srcp[x + 1] + 2) >> 2;

    dstp[width - 1] = (srcp[width - 2] + srcp[width - 1] + 1) >> 1;

    srcp += stride;
    dstp += stride;
  }
}

#ifdef INTEL_INTRINSICS
void HorizontalBlur6_SSE2_simd(const uint8_t* srcp, uint8_t* dstp, int stride, int width, int height)
{
  auto zero = _mm_setzero_si128();
  auto eight = _mm_set1_epi16(8);
  auto six = _mm_set1_epi16(6);

  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x += 16) {
      __m128i pp = _mm_loadu_si128(reinterpret_cast<const __m128i*>(srcp + x - 2));
      __m128i p = _mm_loadu_si128(reinterpret_cast<const __m128i*>(srcp + x - 1));
      __m128i s = _mm_load_si128(reinterpret_cast<const __m128i*>(srcp + x));
      __m128i n = _mm_loadu_si128(reinterpret_cast<const __m128i*>(srcp + x + 1));
      __m128i nn = _mm_loadu_si128(reinterpret_cast<const __m128i*>(srcp + x + 2));

      auto pp_lo = _mm_unpacklo_epi8(pp, zero);
      auto pp_hi = _mm_unpackhi_epi8(pp, zero);
      auto p_lo = _mm_unpacklo_epi8(p, zero);
      auto p_hi = _mm_unpackhi_epi8(p, zero);
      auto s_lo = _mm_unpacklo_epi8(s, zero);
      auto s_hi = _mm_unpackhi_epi8(s, zero);
      auto n_lo = _mm_unpacklo_epi8(n, zero);
      auto n_hi = _mm_unpackhi_epi8(n, zero);
      auto nn_lo = _mm_unpacklo_epi8(nn, zero);
      auto nn_hi = _mm_unpackhi_epi8(nn, zero);
      
      auto centermulsix_lo = _mm_mullo_epi16(s_lo, six);
      auto centermulsix_hi = _mm_mullo_epi16(s_hi, six);
      auto res_lo = _mm_add_epi16(centermulsix_lo, _mm_add_epi16(_mm_add_epi16(pp_lo, _mm_slli_epi16(_mm_add_epi16(p_lo, n_lo), 2)), nn_lo));
      auto res_hi = _mm_add_epi16(centermulsix_hi, _mm_add_epi16(_mm_add_epi16(pp_hi, _mm_slli_epi16(_mm_add_epi16(p_hi, n_hi), 2)), nn_hi));

      res_lo = _mm_srli_epi16(_mm_add_epi16(res_lo, eight), 4);
      res_hi = _mm_srli_epi16(_mm_add_epi16(res_hi, eight), 4);
      auto result = _mm_packus_epi16(res_lo, res_hi);
      _mm_store_si128(reinterpret_cast<__m128i*>(dstp + x), result);
      // dstp[x] = (srcp[x - 2] + ((srcp[x - 1] + srcp[x + 1]) << 2) + srcp[x] * 6 + srcp[x + 2] + 8) >> 4;
    }

    srcp += stride;
    dstp += stride;
  }
}
#endif

void HorizontalBlur6_c(const uint8_t* srcp, uint8_t* dstp, int stride, int width, int height)
{
  for (int y = 0; y < height; y++)
  {
    dstp[0] = (srcp[0] * 6 + (srcp[1] << 3) + (srcp[2] << 1) + 8) >> 4;
    dstp[1] = (((srcp[0] + srcp[2]) << 2) + srcp[1] * 6 + (srcp[3] << 1) + 8) >> 4;

    for (int x = 2; x < width - 2; ++x)
      dstp[x] = (srcp[x - 2] + ((srcp[x - 1] + srcp[x + 1]) << 2) + srcp[x] * 6 + srcp[x + 2] + 8) >> 4;

    dstp[width - 2] = ((srcp[width - 4] << 1) + ((srcp[width - 3] + srcp[width - 1]) << 2) + srcp[width - 2] * 6 + 8) >> 4;
    dstp[width - 1] = ((srcp[width - 3] << 1) + (srcp[width - 2] << 3) + srcp[width - 1] * 6 + 8) >> 4;

    srcp += stride;
    dstp += stride;
  }
}

#ifdef INTEL_INTRINSICS
void andNeighborsInPlace_SSE2_simd(uint8_t* srcp, int stride, int width, int height)
{
  uint8_t* srcpp = srcp - stride;
  uint8_t* srcpn = srcp + stride;

  for (int y = 0; y < height; y++)
  {
    for (int x = 0; x < width; x += 16) {
      __m128i src_0 = _mm_load_si128(reinterpret_cast<const __m128i*>(srcp + x));
      __m128i src_p_m1 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(srcpp + x - 1));
      __m128i src_p = _mm_loadu_si128(reinterpret_cast<const __m128i*>(srcpp + x));
      __m128i src_p_p1 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(srcpp + x + 1));
      __m128i src_n_m1 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(srcpn + x - 1));
      __m128i src_n = _mm_loadu_si128(reinterpret_cast<const __m128i*>(srcpn + x));
      __m128i src_n_p1 = _mm_loadu_si128(reinterpret_cast<const __m128i*>(srcpn + x + 1));
      auto result_p = _mm_or_si128(_mm_or_si128(src_p_m1, src_p), src_p_p1);
      auto result_n = _mm_or_si128(_mm_or_si128(src_n_m1, src_n), src_n_p1);
      auto result = _mm_and_si128(src_0, _mm_or_si128(result_p, result_n));
      _mm_store_si128(reinterpret_cast<__m128i*>(srcp + x), result);
      // srcp[x] &= (srcpp[x - 1] | srcpp[x] | srcpp[x + 1] | srcpn[x - 1] | srcpn[x] | srcpn[x + 1]);
    }

    srcpp += stride;
    srcp += stride;
    srcpn += stride;
  }
}
#endif
