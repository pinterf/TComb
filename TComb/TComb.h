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

#include <windows.h>
#include <stdint.h>
#include <stdio.h>
#include "avisynth.h"
#include "PlanarFrame.h"

#define VERSION "v2.0"


#define min3(a,b,c) min(min(a,b),c)
#define max3(a,b,c) max(max(a,b),c)
#define min4(a,b,c,d) min(min(a,b),min(c,d))
#define max4(a,b,c,d) max(max(a,b),max(c,d))

class TCombFrame
{
public:
	int fnum;
	bool sc;
	bool isValid[11];
	PlanarFrame *orig, *msk1, *msk2;
	PlanarFrame **b, *avg, *omsk;
	TCombFrame::TCombFrame();
	TCombFrame::TCombFrame(VideoInfo &vi, int cpuFlags);
	TCombFrame::~TCombFrame();
	void TCombFrame::setFNum(int i);
};

class TCombCache
{
public:
	TCombFrame **frames;
	int start_pos, size;
	TCombCache::TCombCache();
	TCombCache::TCombCache(int _size, VideoInfo &vi, int cpuFlags);
	TCombCache::~TCombCache();
	void TCombCache::resetCacheStart(int first, int last);
	int TCombCache::getCachePos(int n);
};

class TComb : public GenericVideoFilter
{
public:
	PVideoFrame __stdcall TComb::GetFrame(int n, IScriptEnvironment *env);
	TComb::TComb(PClip _child, int _mode, int _fthreshL, int _fthreshC, int _othreshL,
		int othreshC, bool _map, double _scthresh, bool _debug, int _opt, IScriptEnvironment *env);
	TComb::~TComb();
private:
	bool map, debug;
	int fthreshL, fthreshC;
	int othreshL, othreshC;
	int mode, opt;
	unsigned long diffmaxsc;
	double scthresh;
	PlanarFrame *dstPF, *tmpPF;
	PlanarFrame *minPF, *maxPF;
	PlanarFrame *padPF;
	TCombCache *tdc;
	char buf[256];
	int TComb::mapn(int n);
	void TComb::getAverages(int lc, IScriptEnvironment *env);
	void TComb::buildOscillationMasks(int lc, IScriptEnvironment *env);
	void TComb::getFinalMasks(int lc, IScriptEnvironment *env);
	void TComb::insertFrame(PVideoFrame &src, int pos, int fnum, int lc, IScriptEnvironment *env);
	void TComb::buildDiffMask(TCombFrame *tf1, TCombFrame *tf2, int lc, IScriptEnvironment *env);
	void TComb::buildDiffMasks(int lc, IScriptEnvironment *env);
	void TComb::absDiff(PlanarFrame *src1, PlanarFrame *src2, PlanarFrame *dst, 
		int lc, IScriptEnvironment *env);
	void TComb::absDiffAndMinMask(PlanarFrame *src1, PlanarFrame *src2, PlanarFrame *dst, 
		int lc, IScriptEnvironment *env);
	void TComb::VerticalBlur3(PlanarFrame *src, PlanarFrame *dst, int lc, IScriptEnvironment *env);
	void TComb::HorizontalBlur3(PlanarFrame *src, PlanarFrame *dst, int lc, IScriptEnvironment *env);
	void TComb::getStartStop(int lc, int &start, int &stop);
	void TComb::buildFinalFrame(PlanarFrame *p2, PlanarFrame *p1, PlanarFrame *src, 
		PlanarFrame *n1, PlanarFrame *n2, PlanarFrame *m1, PlanarFrame *m2, PlanarFrame *m3,
		PlanarFrame *dst, int lc, IScriptEnvironment *env);
	void TComb::copyPad(PlanarFrame *src, PlanarFrame *dst, int lc, IScriptEnvironment *env);
	void TComb::MinMax(PlanarFrame *src, PlanarFrame *dmin, PlanarFrame *dmax, int lc, 
		IScriptEnvironment *env);
	void TComb::HorizontalBlur6(PlanarFrame *src, PlanarFrame *dst, int lc, IScriptEnvironment *env);
	void TComb::absDiffAndMinMaskThresh(PlanarFrame *src1, PlanarFrame *src2, PlanarFrame *dst, 
		int lc, IScriptEnvironment *env);
	void TComb::buildFinalMask(PlanarFrame *s1, PlanarFrame *s2, PlanarFrame *m1,
		PlanarFrame *dst, int lc, IScriptEnvironment *env);
	void TComb::calcAverages(PlanarFrame *s1, PlanarFrame *s2, PlanarFrame *dst, int lc,IScriptEnvironment *env);
	void TComb::checkOscillation5(PlanarFrame *p2, PlanarFrame *p1, PlanarFrame *s1,
		PlanarFrame *n1, PlanarFrame *n2, PlanarFrame *dst, int lc, IScriptEnvironment *env);
	void TComb::checkAvgOscCorrelation(PlanarFrame *s1, PlanarFrame *s2, PlanarFrame *s3,
		PlanarFrame *s4, PlanarFrame *dst, int lc, IScriptEnvironment *env);
	void TComb::or3Masks(PlanarFrame *s1, PlanarFrame *s2, PlanarFrame *s3,
		PlanarFrame *dst, int lc, IScriptEnvironment *env);
	void TComb::orAndMasks(PlanarFrame *s1, PlanarFrame *s2, PlanarFrame *dst, int lc, IScriptEnvironment *env);
	void TComb::andMasks(PlanarFrame *s1, PlanarFrame *s2, PlanarFrame *dst, int lc, IScriptEnvironment *env);
	bool TComb::checkSceneChange(PlanarFrame *s1, PlanarFrame *s2, int n, IScriptEnvironment *env);
	void TComb::andNeighborsInPlace(PlanarFrame *src, int lc, IScriptEnvironment *env);
};