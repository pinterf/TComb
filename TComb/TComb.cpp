/*
**                    TComb v2.0 for Avisynth 2.6.x
**
**   TComb is a temporal comb filter (it reduces cross-luminance (rainbowing)
**   and cross-chrominance (dot crawl) artifacts in static areas of the picture).
**   It will ONLY work with NTSC material, and WILL NOT work with telecined material
**   where the rainbowing/dotcrawl was introduced prior to the telecine process!
**   It must be used before ivtc or deinterlace.
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
#include <inttypes.h>

extern "C" void buildFinalMask_SSE2(const uint8_t *s1p, const uint8_t *s2p, const uint8_t *m1p, uint8_t *dstp, int stride, int width, int height, int thresh);
extern "C" void andNeighborsInPlace_SSE2(uint8_t *srcp, int stride, int width, int height);
extern "C" void absDiff_SSE2(const uint8_t *srcp1, const uint8_t *srcp2, uint8_t *dstp, int stride, int width, int height);
extern "C" void absDiffAndMinMask_SSE2(const uint8_t *srcp1, const uint8_t *srcp2, uint8_t *dstp, int stride, int width, int height);
extern "C" void absDiffAndMinMaskThresh_SSE2(const uint8_t *srcp1, const uint8_t *srcp2, uint8_t *dstp, int stride, int width, int height, int thresh);
extern "C" void MinMax_SSE2(const uint8_t *srcp, uint8_t *dstpMin, uint8_t *dstpMax, int src_stride, int dmin_stride, int width, int height, int thresh);
extern "C" void checkOscillation5_SSE2(const uint8_t *p2p, const uint8_t *p1p, const uint8_t *s1p, const uint8_t *n1p, const uint8_t *n2p, uint8_t *dstp, int stride, int width, int height, int thresh);
extern "C" void calcAverages_SSE2(const uint8_t *s1p, const uint8_t *s2p, uint8_t *dstp, int stride, int width, int height);
extern "C" void checkAvgOscCorrelation_SSE2(const uint8_t *s1p, const uint8_t *s2p, const uint8_t *s3p, const uint8_t *s4p, uint8_t *dstp, int stride, int width, int height, int thresh);
extern "C" void or3Masks_SSE2(const uint8_t *s1p, const uint8_t *s2p, const uint8_t *s3p, uint8_t *dstp, int stride, int width, int height);
extern "C" void orAndMasks_SSE2(const uint8_t *s1p, const uint8_t *s2p, uint8_t *dstp, int stride, int width, int height);
extern "C" void andMasks_SSE2(const uint8_t *s1p, const uint8_t *s2p, uint8_t *dstp, int stride, int width, int height);
extern "C" void checkSceneChange_SSE2(const uint8_t *s1p, const uint8_t *s2p, int stride, int width, int height, int64_t *diffp);
extern "C" void VerticalBlur3_SSE2(const uint8_t *srcp, uint8_t *dstp, int stride, int width, int height);
extern "C" void HorizontalBlur3_SSE2(const uint8_t *srcp, uint8_t *dstp, int stride, int width, int height);
extern "C" void HorizontalBlur6_SSE2(const uint8_t *srcp, uint8_t *dstp, int stride, int width, int height);

PVideoFrame __stdcall TComb::GetFrame(int n, IScriptEnvironment *env)
{
	tdc->resetCacheStart(n - 10, n + 10);
	int lc = mode == 2 ? 0x111 : (mode == 1 ? 0x110 : 0x1);
	for (int i = -10; i <= 10; ++i)
	{
		if (tdc->frames[tdc->getCachePos(10 + i)]->fnum != n + i)
			insertFrame(child->GetFrame(mapn(n + i), env), 10 + i, n + i, lc, env);
	}
	if (mode == 0 || mode == 2)
		buildDiffMasks(0x1, env);
	getAverages(lc, env);
	buildOscillationMasks(lc, env);
	getFinalMasks(lc, env);
	buildFinalFrame(
		tdc->frames[tdc->getCachePos(6)]->orig,
		tdc->frames[tdc->getCachePos(8)]->orig,
		tdc->frames[tdc->getCachePos(10)]->orig,
		tdc->frames[tdc->getCachePos(12)]->orig,
		tdc->frames[tdc->getCachePos(14)]->orig,
		tdc->frames[tdc->getCachePos(10)]->msk2,
		tdc->frames[tdc->getCachePos(12)]->msk2,
		tdc->frames[tdc->getCachePos(14)]->msk2,
		dstPF, lc, env);
	PVideoFrame dst = env->NewVideoFrame(vi);
	dstPF->copyTo(dst, vi);
	return dst;
}

void TComb::insertFrame(PVideoFrame &src, int pos, int fnum, int lc, IScriptEnvironment *env)
{
	TCombFrame *tf = tdc->frames[tdc->getCachePos(pos)];
	tf->setFNum(fnum);
	tf->orig->copyFrom(src, vi);
	if (pos > 1)
		tf->sc = checkSceneChange(tf->orig,
		tdc->frames[tdc->getCachePos(pos - 2)]->orig, fnum, env);
	if (lc & 0x1)
	{
		HorizontalBlur3(tf->orig, tf->b[0], 0x1, env);
		VerticalBlur3(tf->orig, tf->b[1], 0x1, env);
		HorizontalBlur3(tf->b[1], tf->b[2], 0x1, env);
		HorizontalBlur6(tf->orig, tf->b[3], 0x1, env);
		VerticalBlur3(tf->b[1], tf->b[4], 0x1, env);
		HorizontalBlur6(tf->b[4], tf->b[5], 0x1, env);
	}
	memset(tf->isValid, 1, 7 * sizeof(bool));
}

void TComb::buildDiffMasks(int lc, IScriptEnvironment *env)
{
	for (int i = 2; i <= 20; i += 2)
	{
		TCombFrame *tf = tdc->frames[tdc->getCachePos(i)];
		if (!tf->isValid[7])
		{
			buildDiffMask(tdc->frames[tdc->getCachePos(i - 2)], tf, lc, env);
			tf->isValid[7] = true;
		}
	}
}

void TComb::buildDiffMask(TCombFrame *tf1, TCombFrame *tf2, int lc, IScriptEnvironment *env)
{
	absDiff(tf1->orig, tf2->orig, tf2->msk1, lc, env);
	for (int i = 0; i<5; ++i)
		absDiffAndMinMask(tf1->b[i], tf2->b[i], tf2->msk1, lc, env);
	absDiffAndMinMaskThresh(tf1->b[5], tf2->b[5], tf2->msk1, lc, env);
}

void TComb::getAverages(int lc, IScriptEnvironment *env)
{
	for (int i = 2; i <= 20; i += 2)
	{
		TCombFrame *tf = tdc->frames[tdc->getCachePos(i)];
		if (!tf->isValid[9])
		{
			calcAverages(tf->orig, tdc->frames[tdc->getCachePos(i - 2)]->orig, tf->avg, lc, env);
			tf->isValid[9] = true;
		}
	}
}

void TComb::buildOscillationMasks(int lc, IScriptEnvironment *env)
{
	for (int i = 8; i <= 20; i += 2)
	{
		TCombFrame *tf = tdc->frames[tdc->getCachePos(i)];
		if (!tf->isValid[10])
		{
			checkOscillation5(tdc->frames[tdc->getCachePos(i - 8)]->orig,
				tdc->frames[tdc->getCachePos(i - 6)]->orig,
				tdc->frames[tdc->getCachePos(i - 4)]->orig,
				tdc->frames[tdc->getCachePos(i - 2)]->orig,
				tf->orig, tf->omsk, lc, env);
			checkAvgOscCorrelation(
				tdc->frames[tdc->getCachePos(i - 6)]->avg,
				tdc->frames[tdc->getCachePos(i - 4)]->avg,
				tdc->frames[tdc->getCachePos(i - 2)]->avg,
				tf->avg, tf->omsk, lc, env);
			tf->isValid[10] = true;
		}
	}
}

void TComb::getFinalMasks(int lc, IScriptEnvironment *env)
{
	for (int i = 10; i <= 14; i += 2)
	{
		TCombFrame *tf = tdc->frames[tdc->getCachePos(i)];
		if (!tf->isValid[8])
		{
			if (tdc->frames[tdc->getCachePos(i - 2)]->sc ||
				tdc->frames[tdc->getCachePos(i - 0)]->sc)
			{
				for (int j = 0; j<3; ++j)
					memset(tf->msk2->GetPtr(j), 0, tf->msk2->GetPitch(j)*tf->msk2->GetHeight(j));
			}
			else
			{
				if (lc & 0x1)
				{
					andMasks(tdc->frames[tdc->getCachePos(i - 2)]->omsk,
						tdc->frames[tdc->getCachePos(i)]->omsk,
						tmpPF, 0x1, env);
					for (int j = 0; j <= 4; j += 2)
						orAndMasks(tdc->frames[tdc->getCachePos(i + j)]->omsk,
						tdc->frames[tdc->getCachePos(i + j + 2)]->omsk,
						tmpPF, 0x1, env);
					andNeighborsInPlace(tmpPF, 0x1, env);
					orAndMasks(tdc->frames[tdc->getCachePos(i - 2)]->msk1,
						tf->msk1, tmpPF, 0x1, env);
				}
				if (lc & 0x110)
				{
					or3Masks(tf->omsk,
						tdc->frames[tdc->getCachePos(i + 2)]->omsk,
						tdc->frames[tdc->getCachePos(i + 4)]->omsk,
						tmpPF, (lc & 0x110), env);
				}
				buildFinalMask(tdc->frames[tdc->getCachePos(i - 4)]->orig,
					tf->orig, tmpPF, tf->msk2, lc, env);
			}
			tf->isValid[8] = true;
		}
	}
}

TComb::TComb(PClip _child, int _mode, int _fthreshL, int _fthreshC, int _othreshL, int _othreshC,
	bool _map, double _scthresh, bool _debug, int _opt, IScriptEnvironment *env) :
	GenericVideoFilter(_child), mode(_mode), fthreshL(_fthreshL), fthreshC(_fthreshC),
	othreshL(_othreshL), othreshC(_othreshC), map(_map), scthresh(_scthresh), debug(_debug), opt(_opt)
{
	dstPF = tmpPF = NULL;
	minPF = maxPF = NULL;
	padPF = NULL;
	tdc = NULL;

	if (!vi.IsYV12() && !vi.IsYUY2())
		env->ThrowError("TComb: Only YV12 and YUY2 colorspace are supported!");
	if (fthreshC < 1 || fthreshC > 255)
		env->ThrowError("TComb: threshC must be in the range 1 <= fthreshC <= 255!");
	if (fthreshL < 1 || fthreshL > 255)
		env->ThrowError("TComb: threshL must be in the range 1 <= fthreshL <= 255!");
	if (othreshC < 1 || othreshC > 255)
		env->ThrowError("TComb: threshC must be in the range 1 <= othreshC <= 255!");
	if (othreshL < 1 || othreshL > 255)
		env->ThrowError("TComb: threshL must be in the range 1 <= othreshL <= 255!");

	child->SetCacheHints(CACHE_NOTHING, 0);

  const int cpuFlags = env->GetCPUFlags();

	dstPF = new PlanarFrame(vi, cpuFlags);
	tmpPF = new PlanarFrame(vi, cpuFlags);
	minPF = new PlanarFrame(vi, cpuFlags);
	maxPF = new PlanarFrame(vi, cpuFlags);
	padPF = new PlanarFrame(cpuFlags);
	padPF->createPlanar(vi.height + 4, vi.width + 4, vi.IsYV12() ? 1 : 2);
	tdc = new TCombCache(21, vi, cpuFlags);

	if (scthresh < 0.0)
		diffmaxsc = 0xFFFFFFFF;
	else
		diffmaxsc = unsigned long(((vi.width >> 4) << 4)*vi.height*219.0*scthresh / 100.0);

	if (debug)
	{
		sprintf_s(buf, "TComb %s by tritical\n", VERSION);
		OutputDebugString(buf);
	}
}

TComb::~TComb()
{
	if (dstPF) delete dstPF;
	if (tmpPF) delete tmpPF;
	if (minPF) delete minPF;
	if (maxPF) delete maxPF;
	if (padPF) delete padPF;
	if (tdc) delete tdc;
}

void TComb::buildFinalFrame(PlanarFrame *p2, PlanarFrame *p1, PlanarFrame *src,
	PlanarFrame *n1, PlanarFrame *n2, PlanarFrame *m1, PlanarFrame *m2, PlanarFrame *m3,
	PlanarFrame *dst, int lc, IScriptEnvironment *env)
{
	if (!map)
		dst->copyFrom(*src);
	else
		for (int b = 0; b < 3; ++b)
			memset(dst->GetPtr(b), 0, dst->GetPitch(b)*dst->GetHeight(b));

	int start = 0, stop = 3;
	getStartStop(lc, start, stop);
	MinMax(src, minPF, maxPF, lc, env);
	for (int b = start; b < stop; ++b)
	{
		const uint8_t *p2p = p2->GetPtr(b);
		const int p2_pitch = p2->GetPitch(b);
		const uint8_t *p1p = p1->GetPtr(b);
		const int p1_pitch = p1->GetPitch(b);
		const uint8_t *srcp = src->GetPtr(b);
		const int src_pitch = src->GetPitch(b);
		const int height = src->GetHeight(b);
		const int width = src->GetWidth(b);
		const uint8_t *n1p = n1->GetPtr(b);
		const int n1_pitch = n1->GetPitch(b);
		const uint8_t *n2p = n2->GetPtr(b);
		const int n2_pitch = n2->GetPitch(b);
		const uint8_t *m1p = m1->GetPtr(b);
		const int m1_pitch = m1->GetPitch(b);
		const uint8_t *m2p = m2->GetPtr(b);
		const int m2_pitch = m2->GetPitch(b);
		const uint8_t *m3p = m3->GetPtr(b);
		const int m3_pitch = m3->GetPitch(b);
		const uint8_t *minp = minPF->GetPtr(b);
		const int min_pitch = minPF->GetPitch(b);
		const uint8_t *maxp = maxPF->GetPtr(b);
		const int max_pitch = maxPF->GetPitch(b);
		uint8_t *dstp = dst->GetPtr(b);
		const int dst_pitch = dst->GetPitch(b);

		if (!map)
		{
			for (int y = 0; y < height; ++y)
			{
				for (int x = 0; x < width; ++x)
				{
					if (m2p[x])
					{
						const int val = (p1p[x] + (srcp[x] << 1) + n1p[x] + 2) >> 2;
						if (val >= minp[x] && val <= maxp[x]) { dstp[x] = val; continue; }
					}

					if (m1p[x])
					{
						const int val = (p2p[x] + (p1p[x] << 1) + srcp[x] + 2) >> 2;
						if (val >= minp[x] && val <= maxp[x]) { dstp[x] = val; continue; }
					}

					if (m3p[x])
					{
						const int val = (srcp[x] + (n1p[x] << 1) + n2p[x] + 2) >> 2;
						if (val >= minp[x] && val <= maxp[x]) { dstp[x] = val; continue; }
					}
				}

				m1p += m1_pitch;
				m2p += m2_pitch;
				m3p += m3_pitch;
				p2p += p2_pitch;
				p1p += p1_pitch;
				srcp += src_pitch;
				n1p += n1_pitch;
				n2p += n2_pitch;
				dstp += dst_pitch;
				minp += min_pitch;
				maxp += max_pitch;
			}
		}
		else
		{
			for (int y = 0; y < height; ++y)
			{
				for (int x = 0; x < width; ++x)
				{
					if (m2p[x])
					{
						const int val = (p1p[x] + (srcp[x] << 1) + n1p[x] + 2) >> 2;
						if (val >= minp[x] && val <= maxp[x]) { dstp[x] = 255; continue; }
					}

					if (m1p[x])
					{
						const int val = (p2p[x] + (p1p[x] << 1) + srcp[x] + 2) >> 2;
						if (val >= minp[x] && val <= maxp[x]) { dstp[x] = 170; continue; }
					}

					if (m3p[x])
					{
						const int val = (srcp[x] + (n1p[x] << 1) + n2p[x] + 2) >> 2;
						if (val >= minp[x] && val <= maxp[x]) { dstp[x] = 85; continue; }
					}
				}

				m1p += m1_pitch;
				m2p += m2_pitch;
				m3p += m3_pitch;
				p2p += p2_pitch;
				p1p += p1_pitch;
				srcp += src_pitch;
				n1p += n1_pitch;
				n2p += n2_pitch;
				dstp += dst_pitch;
				minp += min_pitch;
				maxp += max_pitch;
			}
		}
	}
}

void TComb::buildFinalMask(PlanarFrame *s1, PlanarFrame *s2, PlanarFrame *m1,
	PlanarFrame *dst, int lc, IScriptEnvironment *env)
{
	const int cpu = opt == 0 ? 0 : env->GetCPUFlags();

	int start = 0, stop = 3;
	getStartStop(lc, start, stop);
	for (int b = start; b<stop; ++b)
	{
		const uint8_t *s1p = s1->GetPtr(b);
		const int stride = s1->GetPitch(b);
		const int width = s1->GetWidth(b);
		const int height = s1->GetHeight(b);
		const uint8_t *s2p = s2->GetPtr(b);
		const uint8_t *m1p = m1->GetPtr(b);
		uint8_t *dstp = dst->GetPtr(b);
		const int thresh = b == 0 ? othreshL : othreshC;

		if (cpu&CPUF_SSE2)
			buildFinalMask_SSE2(s1p, s2p, m1p, dstp, stride, width, height, thresh);
		else
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
	}
}

void TComb::andNeighborsInPlace(PlanarFrame *src, int lc, IScriptEnvironment *env)
{
	const int cpu = opt == 0 ? 0 : env->GetCPUFlags();

	uint8_t *srcp = src->GetPtr(0);
	const int height = src->GetHeight(0);
	const int width = src->GetWidth(0);
	const int stride = src->GetPitch(0);
	uint8_t *srcpp = srcp - stride;
	uint8_t *srcpn = srcp + stride;

	srcp[0] &= (srcpn[0] | srcpn[1]);
	for (int x = 1; x < width - 1; ++x)
		srcp[x] &= (srcpn[x - 1] | srcpn[x] | srcpn[x + 1]);
	srcp[width - 1] &= (srcpn[width - 2] | srcpn[width - 1]);

	srcpp += stride;
	srcp += stride;
	srcpn += stride;

	if (cpu&CPUF_SSE2)
	{
		const int widtha = (width % 16) ? ((width >> 4) << 4) : width - 16;

		andNeighborsInPlace_SSE2(srcp + 16, stride, widtha - 16, height - 2);

		for (int y = 1; y < height - 1; y++)
		{
			srcp[0] &= (srcpp[0] | srcpp[1] | srcpn[0] | srcpn[1]);
			for (int x = 1; x < 16; x++)
				srcp[x] &= (srcpp[x - 1] | srcpp[x] | srcpp[x + 1] | srcpn[x - 1] | srcpn[x] | srcpn[x + 1]);

			for (int x = widtha; x < width - 1; x++)
				srcp[x] &= (srcpp[x - 1] | srcpp[x] | srcpp[x + 1] | srcpn[x - 1] | srcpn[x] | srcpn[x + 1]);
			srcp[width - 1] &= (srcpp[width - 2] | srcpp[width - 1] | srcpn[width - 2] | srcpn[width - 1]);

			srcpp += stride;
			srcp += stride;
			srcpn += stride;
		}
	}
	else
	{
		for (int y = 1; y < height - 1; ++y)
		{
			srcp[0] &= (srcpp[0] | srcpp[1] | srcpn[0] | srcpn[1]);
			for (int x = 1; x < width - 1; ++x)
				srcp[x] &= (srcpp[x - 1] | srcpp[x] | srcpp[x + 1] | srcpn[x - 1] | srcpn[x] | srcpn[x + 1]);
			srcp[width - 1] &= (srcpp[width - 2] | srcpp[width - 1] | srcpn[width - 2] | srcpn[width - 1]);

			srcpp += stride;
			srcp += stride;
			srcpn += stride;
		}
	}

	srcp[0] &= (srcpp[0] | srcpp[1]);
	for (int x = 1; x < width - 1; ++x)
		srcp[x] &= (srcpp[x - 1] | srcpp[x] | srcpp[x + 1]);
	srcp[width - 1] &= (srcpp[width - 2] | srcpp[width - 1]);
}

void TComb::absDiff(PlanarFrame *src1, PlanarFrame *src2, PlanarFrame *dst, int lc, IScriptEnvironment *env)
{
	const int cpu = opt == 0 ? 0 : env->GetCPUFlags();

	const uint8_t *srcp1 = src1->GetPtr(0);
	const uint8_t *srcp2 = src2->GetPtr(0);
	uint8_t *dstp = dst->GetPtr(0);
	const int height = src1->GetHeight(0);
	const int width = src1->GetWidth(0);
	const int stride = src1->GetPitch(0);
	const int src2_pitch = src2->GetPitch(0);
	const int dst_pitch = dst->GetPitch(0);
  // sanity check; FIXME: I haved a feeling that strides are not necessarily equal?
  if (src2_pitch != stride)
    env->ThrowError("AbsDiff sources have difference pitches!");

  if (dst_pitch != stride)
    env->ThrowError("AbsDiff source and destination have difference pitches!");

	if (cpu&CPUF_SSE2)
		absDiff_SSE2(srcp1, srcp2, dstp, stride, width, height);
	else
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
}

void TComb::absDiffAndMinMask(PlanarFrame *src1, PlanarFrame *src2, PlanarFrame *dst, int lc, IScriptEnvironment *env)
{
	const int cpu = opt == 0 ? 0 : env->GetCPUFlags();

	const uint8_t *srcp1 = src1->GetPtr(0);
	const uint8_t *srcp2 = src2->GetPtr(0);
	uint8_t *dstp = dst->GetPtr(0);
	const int height = src1->GetHeight(0);
	const int width = src1->GetWidth(0);
	const int stride = src1->GetPitch(0);

	if (cpu&CPUF_SSE2)
		absDiffAndMinMask_SSE2(srcp1, srcp2, dstp, stride, width, height);
	else
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
}

void TComb::absDiffAndMinMaskThresh(PlanarFrame *src1, PlanarFrame *src2, PlanarFrame *dst,
	int lc, IScriptEnvironment *env)
{
	const int cpu = opt == 0 ? 0 : env->GetCPUFlags();

	const uint8_t *srcp1 = src1->GetPtr(0);
	const uint8_t *srcp2 = src2->GetPtr(0);
	uint8_t *dstp = dst->GetPtr(0);
	const int height = src1->GetHeight(0);
	const int width = src1->GetWidth(0);
	const int stride = src1->GetPitch(0);
	const int thresh = fthreshL;

	if (cpu&CPUF_SSE2)
		absDiffAndMinMaskThresh_SSE2(srcp1, srcp2, dstp, stride, width, height, thresh);
	else
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
}

void TComb::copyPad(PlanarFrame *src, PlanarFrame *dst, int lc, IScriptEnvironment *env)
{
	int start = 0, stop = 3;
	getStartStop(lc, start, stop);
	for (int b = start; b < stop; ++b)
	{
		env->BitBlt(dst->GetPtr(b) + dst->GetPitch(b) + 1, dst->GetPitch(b),
			src->GetPtr(b), src->GetPitch(b), src->GetWidth(b), src->GetHeight(b));

		int height = dst->GetHeight(b);
		int width = dst->GetWidth(b);
		if (b == 0)
		{
			height -= 2;
			width -= 2;
		}
		else if (b > 0 && vi.IsYUY2())
			height -= 2;

		int pitch = dst->GetPitch(b);
		uint8_t *dstp = dst->GetPtr(b) + pitch * 2;
		for (int y = 2; y < height - 2; ++y)
		{
			dstp[0] = dstp[1];
			dstp[width - 1] = dstp[width - 2];
			dstp += pitch;
		}

		env->BitBlt(dst->GetPtr(b), dst->GetPitch(b), dst->GetPtr(b) + dst->GetPitch(b),
			dst->GetPitch(b), dst->GetPitch(b), 1);
		env->BitBlt(dst->GetPtr(b) + dst->GetPitch(b) * (height - 1), dst->GetPitch(b),
			dst->GetPtr(b) + dst->GetPitch(b) * (height - 2), dst->GetPitch(b),
			dst->GetPitch(b), 1);
	}
}

void TComb::MinMax(PlanarFrame *src, PlanarFrame *dmin, PlanarFrame *dmax, int lc,
	IScriptEnvironment *env)
{
	const int cpu = opt == 0 ? 0 : env->GetCPUFlags();

	copyPad(src, padPF, lc, env);

	int start = 0, stop = 3;
	getStartStop(lc, start, stop);
	for (int b = start; b < stop; ++b)
	{
		const int src_stride = padPF->GetPitch(b);
		const uint8_t *srcp = padPF->GetPtr(b) + src_stride + 1;
		const int width = src->GetWidth(b);
		const int height = src->GetHeight(b);
		uint8_t *dminp = dmin->GetPtr(b);
		const int min_stride = dmin->GetPitch(b);
		uint8_t *dmaxp = dmax->GetPtr(b);
		const int thresh = b == 0 ? 2 : 8;

		if (cpu&CPUF_SSE2)
			MinMax_SSE2(srcp, dminp, dmaxp, src_stride, min_stride, width, height, thresh);
		else
		{
			const uint8_t *srcpp = srcp - src_stride;
			const uint8_t *srcpn = srcp + src_stride;

			for (int y = 0; y < height; ++y)
			{
				for (int x = 0; x < width; ++x)
				{
					dminp[x] = max(min(min(min(min(srcpp[x - 1], srcpp[x]),
						min(srcpp[x + 1], srcp[x - 1])),
						min(min(srcp[x], srcp[x + 1]),
						min(srcpn[x - 1], srcpn[x]))), srcpn[x + 1]) - thresh, 0);
					dmaxp[x] = min(max(max(max(max(srcpp[x - 1], srcpp[x]),
						max(srcpp[x + 1], srcp[x - 1])),
						max(max(srcp[x], srcp[x + 1]),
						max(srcpn[x - 1], srcpn[x]))), srcpn[x + 1]) + thresh, 255);
				}

				srcpp += src_stride;
				srcp += src_stride;
				srcpn += src_stride;
				dminp += min_stride;
				dmaxp += min_stride;
			}
		}
	}
}

void TComb::checkOscillation5(PlanarFrame *p2, PlanarFrame *p1, PlanarFrame *s1,
	PlanarFrame *n1, PlanarFrame *n2, PlanarFrame *dst, int lc, IScriptEnvironment *env)
{
	const int cpu = opt == 0 ? 0 : env->GetCPUFlags();

	int start = 0, stop = 3;
	getStartStop(lc, start, stop);
	for (int b = start; b < stop; ++b)
	{
		const uint8_t *p1p = p1->GetPtr(b);
		const int stride = p1->GetPitch(b);
		const int width = p1->GetWidth(b);
		const int height = p1->GetHeight(b);
		const uint8_t *p2p = p2->GetPtr(b);
		const uint8_t *s1p = s1->GetPtr(b);
		const uint8_t *n1p = n1->GetPtr(b);
		const uint8_t *n2p = n2->GetPtr(b);
		uint8_t *dstp = dst->GetPtr(b);
		const int thresh = b == 0 ? othreshL : othreshC;

		if (cpu&CPUF_SSE2)
			checkOscillation5_SSE2(p2p, p1p, s1p, n1p, n2p, dstp, stride, width, height, thresh);
		else
		{
			for (int y = 0; y < height; ++y)
			{
				for (int x = 0; x < width; ++x)
				{
					const int min31 = min3(p2p[x], s1p[x], n2p[x]);
					const int max31 = max3(p2p[x], s1p[x], n2p[x]);
					const int min22 = min(p1p[x], n1p[x]);
					const int max22 = max(p1p[x], n1p[x]);
					if (((min31 > max22) || max22 == 0 || (max31 < min22) || max31 == 0) &&
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
	}
}

void TComb::calcAverages(PlanarFrame *s1, PlanarFrame *s2, PlanarFrame *dst, int lc,
	IScriptEnvironment *env)
{
	const int cpu = opt == 0 ? 0 : env->GetCPUFlags();

	int start = 0, stop = 3;
	getStartStop(lc, start, stop);
	for (int b = start; b < stop; ++b)
	{
		const uint8_t *s1p = s1->GetPtr(b);
		const int stride = s1->GetPitch(b);
		const int height = s1->GetHeight(b);
		const int width = s1->GetWidth(b);
		const uint8_t *s2p = s2->GetPtr(b);
		uint8_t *dstp = dst->GetPtr(b);

		if (cpu&CPUF_SSE2)
			calcAverages_SSE2(s1p, s2p, dstp, stride, width, height);
		else
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
	}
}

void TComb::checkAvgOscCorrelation(PlanarFrame *s1, PlanarFrame *s2, PlanarFrame *s3,
	PlanarFrame *s4, PlanarFrame *dst, int lc, IScriptEnvironment *env)
{
	const int cpu = opt == 0 ? 0 : env->GetCPUFlags();

	int start = 0, stop = 3;
	getStartStop(lc, start, stop);
	for (int b = start; b < stop; ++b)
	{
		const uint8_t *s1p = s1->GetPtr(b);
		const int stride = s1->GetPitch(b);
		const int width = s1->GetWidth(b);
		const int height = s1->GetHeight(b);
		const uint8_t *s2p = s2->GetPtr(b);
		const uint8_t *s3p = s3->GetPtr(b);
		const uint8_t *s4p = s4->GetPtr(b);
		uint8_t *dstp = dst->GetPtr(b);
		const int thresh = b == 0 ? fthreshL : fthreshC;

		if (cpu&CPUF_SSE2)
			checkAvgOscCorrelation_SSE2(s1p, s2p, s3p, s4p, dstp, stride, width, height, thresh);
		else
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
	}
}

void TComb::or3Masks(PlanarFrame *s1, PlanarFrame *s2, PlanarFrame *s3,
	PlanarFrame *dst, int lc, IScriptEnvironment *env)
{
	const int cpu = opt == 0 ? 0 : env->GetCPUFlags();

	int start = 0, stop = 3;
	getStartStop(lc, start, stop);
	for (int b = start; b < stop; ++b)
	{
		const uint8_t *s1p = s1->GetPtr(b);
		const int stride = s1->GetPitch(b);
		const int width = s1->GetWidth(b);
		const int height = s1->GetHeight(b);
		const uint8_t *s2p = s2->GetPtr(b);
		const uint8_t *s3p = s3->GetPtr(b);
		uint8_t *dstp = dst->GetPtr(b);

		if (cpu&CPUF_SSE2)
			or3Masks_SSE2(s1p, s2p, s3p, dstp, stride, width, height);
		else
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
	}
}

void TComb::orAndMasks(PlanarFrame *s1, PlanarFrame *s2, PlanarFrame *dst, int lc,
	IScriptEnvironment *env)
{
	const int cpu = opt == 0 ? 0 : env->GetCPUFlags();

	int start = 0, stop = 3;
	getStartStop(lc, start, stop);
	for (int b = start; b < stop; ++b)
	{
		const uint8_t *s1p = s1->GetPtr(b);
		const int stride = s1->GetPitch(b);
		const int width = s1->GetWidth(b);
		const int height = s1->GetHeight(b);
		const uint8_t *s2p = s2->GetPtr(b);
		uint8_t *dstp = dst->GetPtr(b);

		if (cpu&CPUF_SSE2)
			orAndMasks_SSE2(s1p, s2p, dstp, stride, width, height);
		else
		{
			for (int y = 0; y<height; ++y)
			{
				for (int x = 0; x<width; ++x)
					dstp[x] |= (s1p[x] & s2p[x]);

				s1p += stride;
				s2p += stride;
				dstp += stride;
			}
		}
	}
}

void TComb::andMasks(PlanarFrame *s1, PlanarFrame *s2, PlanarFrame *dst, int lc,
	IScriptEnvironment *env)
{
	const int cpu = opt == 0 ? 0 : env->GetCPUFlags();

	int start = 0, stop = 3;
	getStartStop(lc, start, stop);
	for (int b = start; b < stop; ++b)
	{
		const uint8_t *s1p = s1->GetPtr(b);
		const int stride = s1->GetPitch(b);
		const int width = s1->GetWidth(b);
		const int height = s1->GetHeight(b);
		const uint8_t *s2p = s2->GetPtr(b);
		uint8_t *dstp = dst->GetPtr(b);

		if (cpu&CPUF_SSE2)
			andMasks_SSE2(s1p, s2p, dstp, stride, width, height);
		else
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
	}
}

bool TComb::checkSceneChange(PlanarFrame *s1, PlanarFrame *s2, int n, IScriptEnvironment *env)
{
	if (scthresh < 0.0)
	{
		if (debug)
		{
			sprintf_s(buf, "TComb: Frame %d: Not a SceneChange\n", n);
			OutputDebugString(buf);
		}
		return false;
	}

	const int cpu = opt == 0 ? 0 : env->GetCPUFlags();

	const uint8_t *s1p = s1->GetPtr(0);
	const uint8_t *s2p = s2->GetPtr(0);
	const int height = s1->GetHeight(0);
	const int width = (s1->GetWidth(0) >> 4) << 4;
	const int stride = s1->GetPitch(0);
	int64_t diff = 0;

	if (cpu&CPUF_SSE2)
		checkSceneChange_SSE2(s1p, s2p, stride, width, height, &diff);
	else
	{
		for (int y = 0; y < height; ++y)
		{
			for (int x = 0; x < width; x += 4)
			{
				diff += abs(s1p[x + 0] - s2p[x + 0]);
				diff += abs(s1p[x + 1] - s2p[x + 1]);
				diff += abs(s1p[x + 2] - s2p[x + 2]);
				diff += abs(s1p[x + 3] - s2p[x + 3]);
			}

			s1p += stride;
			s2p += stride;
		}
	}

	if (diff > diffmaxsc)
	{
		if (debug)
		{
			sprintf_s(buf, "TComb: Frame %d: SceneChange detected! (%" PRId64 ",%u)\n", n, diff, diffmaxsc);
			OutputDebugString(buf);
		}

		return true;
	}

	if (debug)
	{
		sprintf_s(buf, "TComb: Frame %d: Not a SceneChange (%" PRId64 ",%u)\n", n, diff, diffmaxsc);
		OutputDebugString(buf);
	}

	return false;
}

void TComb::VerticalBlur3(PlanarFrame *src, PlanarFrame *dst, int lc, IScriptEnvironment *env)
{
	const int cpu = opt == 0 ? 0 : env->GetCPUFlags();

	const uint8_t *srcp = src->GetPtr(0);
	uint8_t *dstp = dst->GetPtr(0);
	const int stride = src->GetPitch(0);
	const int width = src->GetWidth(0);
	const int height = src->GetHeight(0);

	if (cpu&CPUF_SSE2)
		VerticalBlur3_SSE2(srcp, dstp, stride, width, height);
	else
	{
		const uint8_t *srcpp = srcp - stride;
		const uint8_t *srcpn = srcp + stride;

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
}

void TComb::HorizontalBlur3(PlanarFrame *src, PlanarFrame *dst, int lc, IScriptEnvironment *env)
{
	const int cpu = opt == 0 ? 0 : env->GetCPUFlags();

	const uint8_t *srcp = src->GetPtr(0);
	uint8_t *dstp = dst->GetPtr(0);
	const int stride = src->GetPitch(0);
	const int width = src->GetWidth(0);
	const int height = src->GetHeight(0);

	if ((cpu&CPUF_SSE2) && width >= 16)
	{
		const int widtha = (width >> 4) << 4;

		HorizontalBlur3_SSE2(srcp + 16, dstp + 16, stride, widtha - 32, height);

		for (int y = 0; y < height; y++)
		{
			dstp[0] = (srcp[0] + srcp[1] + 1) >> 1;

			for (int x = 1; x < 16; x++)
				dstp[x] = (srcp[x - 1] + (srcp[x] << 1) + srcp[x + 1] + 2) >> 2;

			for (int x = widtha - 16; x < width - 1; x++)
				dstp[x] = (srcp[x - 1] + (srcp[x] << 1) + srcp[x + 1] + 2) >> 2;

			dstp[width - 1] = (srcp[width - 2] + srcp[width - 1] + 1) >> 1;

			srcp += stride;
			dstp += stride;
		}
	}
	else
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
}

void TComb::HorizontalBlur6(PlanarFrame *src, PlanarFrame *dst, int lc, IScriptEnvironment *env)
{
	const int cpu = opt == 0 ? 0 : env->GetCPUFlags();

	const uint8_t *srcp = src->GetPtr(0);
	uint8_t *dstp = dst->GetPtr(0);
	const int stride = src->GetPitch(0);
	const int width = src->GetWidth(0);
	const int height = src->GetHeight(0);

	if ((cpu&CPUF_SSE2) && width >= 16)
	{
		const int widtha = (width >> 4) << 4;

		HorizontalBlur6_SSE2(srcp + 16, dstp + 16, stride, widtha - 32, height);

		for (int y = 0; y < height; y++)
		{
			dstp[0] = (srcp[0] * 6 + (srcp[1] << 4) + (srcp[2] << 1) + 8) >> 4;
			dstp[1] = (((srcp[0] + srcp[2]) >> 2) + srcp[1] * 6 + (srcp[3] << 1) + 8) >> 4;

			for (int x = 2; x < 16; x++)
				dstp[x] = (srcp[x - 2] + ((srcp[x - 1] + srcp[x + 1]) << 2) + srcp[x] * 6 + srcp[x + 2] + 8) >> 4;

			for (int x = widtha - 16; x < width - 2; x++)
				dstp[x] = (srcp[x - 2] + ((srcp[x - 1] + srcp[x + 1]) << 2) + srcp[x] * 6 + srcp[x + 2] + 8) >> 4;

			dstp[width - 2] = ((srcp[width - 4] << 1) + ((srcp[width - 3] + srcp[width - 1]) << 2) + srcp[width - 2] * 6 + 8) >> 4;
			dstp[width - 1] = ((srcp[width - 3] << 1) + (srcp[width - 2] << 3) + srcp[width - 1] * 6 + 8) >> 4;

			srcp += stride;
			dstp += stride;
		}
	}
	else
	{
		dstp[0] = (srcp[0] * 6 + (srcp[1] << 4) + (srcp[2] << 1) + 8) >> 4;
		dstp[1] = (((srcp[0] + srcp[2]) << 2) + srcp[1] * 6 + (srcp[3] << 1) + 8) >> 4;

		for (int x = 2; x < width - 2; ++x)
			dstp[x] = (srcp[x - 2] + ((srcp[x - 1] + srcp[x + 1]) << 2) + srcp[x] * 6 + srcp[x + 2] + 8) >> 4;

		dstp[width - 2] = ((srcp[width - 4] << 1) + ((srcp[width - 3] + srcp[width - 1]) << 2) + srcp[width - 2] * 6 + 8) >> 4;
		dstp[width - 1] = ((srcp[width - 3] << 1) + (srcp[width - 2] << 3) + srcp[width - 1] * 6 + 8) >> 4;

		srcp += stride;
		dstp += stride;
	}
}

void TComb::getStartStop(int lc, int &start, int &stop)
{
	switch (lc)
	{
	case 0x0:
		stop = 0;
		return;
	case 0x1:
		stop = 1;
		return;
	case 0x10:
		start = 1;
		stop = 2;
		return;
	case 0x11:
		stop = 2;
		return;
	case 0x100:
		start = 2;
		return;
	case 0x110:
		start = 1;
		return;
	case 0x101:
	case 0x111:
		return;
	default:
		return;
	}
}

int TComb::mapn(int n)
{
	if (n < 0) return 0;
	if (n >= vi.num_frames) return vi.num_frames - 1;
	return n;
}

TCombFrame::TCombFrame()
{
	fnum = -20;
	sc = false;
	memset(isValid, 0, 11 * sizeof(bool));
	orig = msk1 = msk2 = omsk = avg = NULL;
	b = NULL;
}

TCombFrame::TCombFrame(VideoInfo &vi, int cpuFlags)
{
	fnum = -20;
	sc = false;
	memset(isValid, 0, 11 * sizeof(bool));
	orig = new PlanarFrame(vi, cpuFlags);
	msk1 = new PlanarFrame(vi, cpuFlags);
	msk2 = new PlanarFrame(vi, cpuFlags);
	omsk = new PlanarFrame(vi, cpuFlags);
	avg = new PlanarFrame(vi, cpuFlags);
	b = (PlanarFrame **)malloc(6 * sizeof(PlanarFrame*));
	for (int i = 0; i<6; ++i)
	{
		b[i] = new PlanarFrame(vi, cpuFlags);
	}
}

TCombFrame::~TCombFrame()
{
	if (orig) delete orig;
	if (msk1) delete msk1;
	if (msk2) delete msk2;
	if (avg) delete avg;
	if (omsk) delete omsk;
	if (b)
	{
		for (int i = 0; i<6; ++i)
		{
			if (b[i]) delete b[i];
		}
		free(b);
	}
}

void TCombFrame::setFNum(int i)
{
	memset(isValid, 0, 11 * sizeof(bool));
	sc = false;
	fnum = i;
}

TCombCache::TCombCache()
{
	frames = NULL;
	start_pos = size = -20;
}

TCombCache::TCombCache(int _size, VideoInfo &vi, int cpuFlags)
{
	frames = NULL;
	start_pos = size = -20;
	if (_size > 0)
	{
		start_pos = 0;
		size = _size;
		frames = (TCombFrame**)malloc(size*sizeof(TCombFrame*));
		memset(frames, 0, size*sizeof(TCombFrame*));
		for (int i = 0; i<size; ++i)
			frames[i] = new TCombFrame(vi, cpuFlags);
	}
}

TCombCache::~TCombCache()
{
	if (frames)
	{
		for (int i = 0; i<size; ++i)
		{
			if (frames[i]) delete frames[i];
		}
		free(frames);
	}
}

void TCombCache::resetCacheStart(int first, int last)
{
	for (int j = first; j <= last; ++j)
	{
		for (int i = 0; i<size; ++i)
		{
			if (frames[i]->fnum == j)
			{
				start_pos = i - j + first;
				if (start_pos < 0) start_pos += size;
				else if (start_pos >= size) start_pos -= size;
				return;
			}
		}
	}
}

int TCombCache::getCachePos(int n)
{
	return (start_pos + n) % size;
}

AVSValue __cdecl Create_TComb(AVSValue args, void* user_data, IScriptEnvironment* env)
{
	if (!args[0].IsClip())
		env->ThrowError("TComb: First argument must be a clip!");
	if (args[0].AsClip()->GetVideoInfo().IsFieldBased())
		env->ThrowError("TComb: Input clip is already field-based! Use AssumeFrameBased().");
	PClip sepf;
	try
	{
		sepf = env->Invoke("SeparateFields", args[0].AsClip()).AsClip();
	}
	catch (IScriptEnvironment::NotFound)
	{
		env->ThrowError("TComb: Error on invoke; SeparateFields was not found!");
	}
	catch (AvisynthError e)
	{
		env->ThrowError("TComb: Caught an AvisynthError on invoke (%s)!", e.msg);
	}
	sepf = new TComb(sepf, args[1].AsInt(2), args[2].AsInt(4), args[3].AsInt(5), args[4].AsInt(5),
		args[5].AsInt(6), args[6].AsBool(false), args[7].AsFloat(12.0), args[8].AsBool(false), args[9].AsInt(-1), env);
	try
	{
		sepf = env->Invoke("Weave", sepf).AsClip();
	}
	catch (IScriptEnvironment::NotFound)
	{
		env->ThrowError("TComb: Error on invoke; Weave was not found!");
	}
	catch (AvisynthError e)
	{
		env->ThrowError("TComb: Caught an AvisynthError on invoke (%s)!", e.msg);
	}
	return sepf;
}

const AVS_Linkage *AVS_linkage = nullptr;

extern "C" __declspec(dllexport) const char* __stdcall AvisynthPluginInit3(IScriptEnvironment* env, const AVS_Linkage* const vectors)
{
	AVS_linkage = vectors;

	env->AddFunction("TComb", "c[mode]i[fthreshL]i[fthreshC]i[othreshL]i[othreshC]i[map]b" \
		"[scthresh]f[debug]b[opt]i", Create_TComb, 0);

	return "TComb plugin";
}