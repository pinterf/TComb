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

#if defined(_WIN32) && !defined(INTEL_INTRINSICS)
#error Forgot to set INTEL_INTRINSICS? Comment out this line if not
#endif

#include "avisynth.h"
#include "common.h"
#include <stdint.h>
#include <stdio.h>
#include "PlanarFrame.h"

#define VERSION "v2.2"

//#define OLD_ASM

#define min3(a,b,c) std::min(std::min(a,b),c)
#define max3(a,b,c) std::max(std::max(a,b),c)
#define min4(a,b,c,d) std::min(std::min(a,b),std::min(c,d))
#define max4(a,b,c,d) std::max(std::max(a,b),std::max(c,d))

class TCombFrame
{
public:
  int fnum;
  bool sc;
  bool isValid[11];
  PlanarFrame* orig, * msk1, * msk2;
  PlanarFrame** b, * avg, * omsk;
  TCombFrame();
  TCombFrame(VideoInfo& vi, int cpuFlags);
  ~TCombFrame();
  void setFNum(int i);
};

class TCombCache
{
public:
  TCombFrame** frames;
  int start_pos, size;
  TCombCache();
  TCombCache(int _size, VideoInfo& vi, int cpuFlags);
  ~TCombCache();
  void resetCacheStart(int first, int last);
  int getCachePos(int n);
};

class TComb : public GenericVideoFilter
{
public:
  PVideoFrame __stdcall GetFrame(int n, IScriptEnvironment* env);
  TComb(PClip _child, int _mode, int _fthreshL, int _fthreshC, int _othreshL,
    int othreshC, bool _map, double _scthresh, bool _debug, int _opt, IScriptEnvironment* env);
  ~TComb();
private:
  bool map, debug;
  int fthreshL, fthreshC;
  int othreshL, othreshC;
  int mode, opt;
  unsigned long diffmaxsc;
  double scthresh;
  PlanarFrame* dstPF, * tmpPF;
  PlanarFrame* minPF, * maxPF;
  PlanarFrame* padPF;
  TCombCache* tdc;
  char buf[256];
  int mapn(int n);
  void getAverages(int lc, IScriptEnvironment* env);
  void buildOscillationMasks(int lc, IScriptEnvironment* env);
  void getFinalMasks(int lc, IScriptEnvironment* env);
  void insertFrame(PVideoFrame& src, int pos, int fnum, int lc, IScriptEnvironment* env);
  void buildDiffMask(TCombFrame* tf1, TCombFrame* tf2, int lc, IScriptEnvironment* env);
  void buildDiffMasks(int lc, IScriptEnvironment* env);
  void absDiff(PlanarFrame* src1, PlanarFrame* src2, PlanarFrame* dst,
    int lc, IScriptEnvironment* env);
  void absDiffAndMinMask(PlanarFrame* src1, PlanarFrame* src2, PlanarFrame* dst,
    int lc, IScriptEnvironment* env);
  void VerticalBlur3(PlanarFrame* src, PlanarFrame* dst, int lc, IScriptEnvironment* env);
  void HorizontalBlur3(PlanarFrame* src, PlanarFrame* dst, int lc, IScriptEnvironment* env);
  void getStartStop(int lc, int& start, int& stop);
  void buildFinalFrame(PlanarFrame* p2, PlanarFrame* p1, PlanarFrame* src,
    PlanarFrame* n1, PlanarFrame* n2, PlanarFrame* m1, PlanarFrame* m2, PlanarFrame* m3,
    PlanarFrame* dst, int lc, IScriptEnvironment* env);
  void copyPad(PlanarFrame* src, PlanarFrame* dst, int lc, IScriptEnvironment* env);
  void MinMax(PlanarFrame* src, PlanarFrame* dmin, PlanarFrame* dmax, int lc,
    IScriptEnvironment* env);
  void HorizontalBlur6(PlanarFrame* src, PlanarFrame* dst, int lc, IScriptEnvironment* env);
  void absDiffAndMinMaskThresh(PlanarFrame* src1, PlanarFrame* src2, PlanarFrame* dst,
    int lc, IScriptEnvironment* env);
  void buildFinalMask(PlanarFrame* s1, PlanarFrame* s2, PlanarFrame* m1,
    PlanarFrame* dst, int lc, IScriptEnvironment* env);
  void calcAverages(PlanarFrame* s1, PlanarFrame* s2, PlanarFrame* dst, int lc, IScriptEnvironment* env);
  void checkOscillation5(PlanarFrame* p2, PlanarFrame* p1, PlanarFrame* s1,
    PlanarFrame* n1, PlanarFrame* n2, PlanarFrame* dst, int lc, IScriptEnvironment* env);
  void checkAvgOscCorrelation(PlanarFrame* s1, PlanarFrame* s2, PlanarFrame* s3,
    PlanarFrame* s4, PlanarFrame* dst, int lc, IScriptEnvironment* env);
  void or3Masks(PlanarFrame* s1, PlanarFrame* s2, PlanarFrame* s3,
    PlanarFrame* dst, int lc, IScriptEnvironment* env);
  void orAndMasks(PlanarFrame* s1, PlanarFrame* s2, PlanarFrame* dst, int lc, IScriptEnvironment* env);
  void andMasks(PlanarFrame* s1, PlanarFrame* s2, PlanarFrame* dst, int lc, IScriptEnvironment* env);
  bool checkSceneChange(PlanarFrame* s1, PlanarFrame* s2, int n, IScriptEnvironment* env);
  void andNeighborsInPlace(PlanarFrame* src, int lc, IScriptEnvironment* env);
};

void checkSceneChangePlanar_1_SSE2_simd(const uint8_t* prvp, const uint8_t* srcp,
  int height, int width, int prv_pitch, int src_pitch, uint64_t& diffp);

template<typename pixel_t>
void checkSceneChangePlanar_1_c(const pixel_t* prvp, const pixel_t* srcp,
  int height, int width, int prv_pitch, int src_pitch, uint64_t& diffp);

void andMasks_SSE2_simd(const uint8_t* s1p, const uint8_t* s2p, uint8_t* dstp, int stride, int width, int height);
void andMasks_c(const uint8_t* s1p, const uint8_t* s2p, uint8_t* dstp, int stride, int width, int height);

void orAndMasks_SSE2_simd(const uint8_t* s1p, const uint8_t* s2p, uint8_t* dstp, int stride, int width, int height);
void orAndMasks_c(const uint8_t* s1p, const uint8_t* s2p, uint8_t* dstp, int stride, int width, int height);

void or3Masks_SSE2_simd(const uint8_t * s1p, const uint8_t * s2p, const uint8_t * s3p, uint8_t * dstp, int stride, int width, int height);
void or3Masks_c(const uint8_t* s1p, const uint8_t* s2p, const uint8_t* s3p, uint8_t* dstp, int stride, int width, int height);

void calcAverages_SSE2_simd(const uint8_t * s1p, const uint8_t * s2p, uint8_t * dstp, int stride, int width, int height);
void calcAverages_c(const uint8_t* s1p, const uint8_t* s2p, uint8_t* dstp, int stride, int width, int height);

void MinMax_SSE2_simd(const uint8_t * srcp, uint8_t * dstpMin, uint8_t * dstpMax, int src_stride, int dmin_stride, int width, int height, int thresh);
void MinMax_c(const uint8_t* srcp, uint8_t* dstpMin, uint8_t* dstpMax, int src_stride, int dmin_stride, int width, int height, int thresh);

void absDiff_SSE2_simd(const uint8_t * srcp1, const uint8_t * srcp2, uint8_t * dstp, int stride, int width, int height);
void absDiff_c(const uint8_t* srcp1, const uint8_t* srcp2, uint8_t* dstp, int stride, int width, int height);

void buildFinalMask_SSE2_simd(const uint8_t * s1p, const uint8_t * s2p, const uint8_t * m1p, uint8_t * dstp, int stride, int width, int height, int thresh);
void buildFinalMask_c(const uint8_t * s1p, const uint8_t * s2p, const uint8_t * m1p, uint8_t * dstp, int stride, int width, int height, int thresh);

void checkOscillation5_SSE2_simd(const uint8_t * p2p, const uint8_t * p1p, const uint8_t * s1p, const uint8_t * n1p, const uint8_t * n2p, uint8_t * dstp, int stride, int width, int height, int thresh);
void checkOscillation5_c(const uint8_t * p2p, const uint8_t * p1p, const uint8_t * s1p, const uint8_t * n1p, const uint8_t * n2p, uint8_t * dstp, int stride, int width, int height, int thresh);

void absDiffAndMinMaskThresh_SSE2_simd(const uint8_t * srcp1, const uint8_t * srcp2, uint8_t * dstp, int stride, int width, int height, int thresh);
void absDiffAndMinMaskThresh_c(const uint8_t * srcp1, const uint8_t * srcp2, uint8_t * dstp, int stride, int width, int height, int thresh);

void absDiffAndMinMask_SSE2_simd(const uint8_t * srcp1, const uint8_t * srcp2, uint8_t * dstp, int stride, int width, int height);
void absDiffAndMinMask_c(const uint8_t * srcp1, const uint8_t * srcp2, uint8_t * dstp, int stride, int width, int height);

void checkAvgOscCorrelation_SSE2_simd(const uint8_t * s1p, const uint8_t * s2p, const uint8_t * s3p, const uint8_t * s4p, uint8_t * dstp, int stride, int width, int height, int thresh);
void checkAvgOscCorrelation_c(const uint8_t * s1p, const uint8_t * s2p, const uint8_t * s3p, const uint8_t * s4p, uint8_t * dstp, int stride, int width, int height, int thresh);

void VerticalBlur3_SSE2_simd(const uint8_t * srcp, uint8_t * dstp, int stride, int width, int height);
void VerticalBlur3_c(const uint8_t* srcp, uint8_t* dstp, int stride, int width, int height);

void HorizontalBlur3_SSE2_simd(const uint8_t * srcp, uint8_t * dstp, int stride, int width, int height);
void HorizontalBlur3_c(const uint8_t* srcp, uint8_t* dstp, int stride, int width, int height);

void HorizontalBlur6_SSE2_simd(const uint8_t* srcp, uint8_t* dstp, int stride, int width, int height);
void HorizontalBlur6_c(const uint8_t * srcp, uint8_t * dstp, int stride, int width, int height);

void andNeighborsInPlace_SSE2_simd(uint8_t * srcp, int stride, int width, int height);
// no distinct C here

