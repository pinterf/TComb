                                                                                                      |
                                    TComb for AviSynth v2.6.x                                         |
                                       v2.1 (22 February 2021)                                        |
                                           by tritical                                                |
                                       modified by Elegant (v2.0; 17 July 2015)                       |
                                       additional work by pinterf                                     |
                                                                                                      |
                                            HELP FILE                                                 |
-------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------


INFO:


      TComb is a temporal comb filter (it reduces cross-luminance (rainbowing) and cross-chrominance
   (dot crawl) artifacts in static areas of the picture).  It will ONLY work with NTSC material, and
   WILL NOT work with telecined material where the rainbowing/dotcrawl was introduced prior to the
   telecine process!  It must be used before ivtc or deinterlace in order to work.  In terms of what
   it does it is similar to guavacomb/dedot. TComb currently supports YV12 and YUY2 colorspaces.

      TComb does support seeking... that is, jumping to a random frame will produce the same result
   as if you had linearly run up to that frame.  For dot crawl removal tcomb requires at least 3
   static fields of the same parity and for rainbow removal tcomb requires at least 5 static fields
   of the same parity.


   Syntax =>

      TComb(int mode, int fthreshL, int fthreshC, int othreshL, int othreshC, bool map,
              double scthresh, bool debug, int opt)



PARAMETERS:


   mode - (limit processing to luma or chroma only)

       Controls whether both luma/chroma are processed or only one or the other.  Possible settings:

           0 - process luma only    (dot crawl removal)
           1 - process chroma only  (rainbow removal)
           2 - process both

       default:  2  (int)


   fthreshL/fthreshC - (filtered pixel correlation thresholds)

       One of the things TComb checks for is correlation between filtered values over the length
       of the filtering window.  If all values differ by less than fthreshL (for luma) or fthreshC
       (for chroma) then the filtered values are considered to be correlated.  Larger values will
       allow more filtering (will be more effective at removing rainbowing/dot crawl), but will also
       create more artifacts.  Smaller values will produce less artifacts, but will be less effective
       in removing rainbowing/dot crawl. A good range of values is between 4 and 7.

       default:  fthreshL -> 4  (int)
                 fthreshC -> 5


   othreshL/othreshC - (original pixel correlation thresholds)

       One of the things TComb checks for is correlation between original pixel values from every
       other field of the same parity.  Due to the oscillation period, these values should be equal
       or very similar in static areas containing dot crawl or rainbowing.  If the pixel values
       differ by less than othreshL (for luma) or othreshC (for chroma) then the pixels are considered
       to be correlated.  Larger values will allow more filtering (will be more effective at removing
       rainbowing/dotcrawl), but will also create more artifacts.  Smaller values will produce less
       artifacts, but will be less effective in removing rainbowing/dotcrawl. A good range of values
       is between 4 and 8.

       default:  othreshL -> 5  (int)
                 othreshC -> 6


   map -

       Identifies pixels that are being replaced with filtered values.  Each pixel in the output
       frame will have one of the following values indicating how it is being filtered:

            0 - not being filtered
           85 - [1 2 1] average of (n,n+1,n+2)
          170 - [1 2 1] average of (n-2,n-1,n)
          255 - [1 2 1] average of (n-1,n,n+1)

           ** n = current frame

       default:  false  (bool)


   scthresh - (scenechange threshold)

       Sets the scenechange detection threshold as a percentage of maximum change on the luma
       plane.  Use the debug output to see which frames are detected as scenechanges and the
       scenechange statistics.

       default:  12.0  (float)


   debug -

       Will enable debug output.  The only thing it shows are the scenechange stats.  The info
       is output via OutputDebugString().  You can use the utility "DebugView" from sysinternals
       to view the output.  The frame numbers in the debug output correspond to the input clip
       after a separatefields() call.  TComb internally invokes separatefields() before itself
       and weave() after itself.

       default:  false  (bool)

   opt - (another debug parameter: CPU)

       0: C only (no assembly at all)
       other: automatically choose SSE2 or C

       For development use: opt parameters can appear/disappear/change their meaning between versions

       default:  -1 (int)


BASIC SETUP/USAGE:


   Setting up TComb is pretty simple.  The only values that would ever really need adjusting
   are fthreshL/fthreshC, othreshL/othreshC, and mode.

   Set mode to 0 if you want to do dot crawl removal only, set it to 1 if you want to
   do rainbow removal only, or set it to 2 to do both.

   Dot Crawl Removal Tweaking (fthreshL/othreshL):

      To find good values for fthreshL/othreshL, start with the following line:

             tcomb(mode=0,fthreshL=255,othreshL=255)

      Now, keep othreshL at 255 but set fthreshL down to 1.  Keep increasing fthreshL
      in steps of 1 to 2 until you find the point at which all dot crawl is removed.
      Remember that value.  Next, set fthreshL back to 255, and set othreshL to 1.
      Now, increase othreshL in steps of 1 or 2 until you find the point at which all
      dot crawl is removed.  You've now got values for fthreshL/othreshL.

   Rainbowing Removal Tweaking (fthreshC/othreshC):

      To find good values for fthreshC/othreshC, start with the following line:

             tcomb(mode=1,fthreshC=255,othreshC=255)

      Now, keep othreshC at 255 but set fthreshC down to 1.  Keep increasing fthreshC
      in steps of 1 to 2 until you find the point at which all (or most) rainbowing is
      removed.  Remember that value.  Next, set fthreshC back to 255, and set othreshC
      to 1. Now, increase othreshC in steps of 1 or 2 until you find the point at which
      all (or most) rainbowing is removed.  You've now got values for fthreshC/othreshC.

   Once you've got values for mode, fthreshL/fthreshC, and othreshL/othreshC, add the
   necessary tcomb() line into your script and run through part of it.  If you see any
   artifacts try lowering your fthresh/othresh values.



CHANGE LIST:

   02/22/2021  v2.1 (pinterf)
       + project forked to https://github.com/pinterf/TComb/ 
       + param 'opt' is back for debug. 0 means pure C code
       + Fix bug in x64 assembler buildFinalMask_SSE2
       + Fix crash in 32bit version of VerticalBlur3_SSE2
       + Fix: scenechange SSE2 did not work
       + Fix: x64 assembler HorizontalBlur6_SSE2
       + Fix: HorizontalBlur6: C only did top 2 lines. SSE2 bad top 2 lines
       + Fix: HorizontalBlur3_SSE2 artifacts (both x86 and x64)
       + Fix: HorizontalBlur3_SSE2 missing rounder (both x86 and x64)
       (now C and SSE2 is giving identical results)
       + Code:
           + Update to Visual Studio 2019
           + update to actual Avisynth+ headers
           + clang-friendly code
           + removed memcpy and bitblt variants
           + replaced planarframes module with the one I updated in tivtc project for 
             avisynth+ and hbd preparation    
           + Fix debug build configuration in VS project settings

   07/26/2015  v2.0.0.1 (Elegant)
       + Corrected the masks used in HorizontalBlur6 for x64.

	   
   07/17/2015  v2.0 (Elegant)

       + Removed buffering of frames/info that weren't actually used (was there for
         development/testing purposes). Should save a lot of RAM usage.
       + Switched to AVS 2.6 API since AviSynth 2.6.0 was released.
	   + Added x64 support which also utilizes SSE2. This also includes some missing
	     SSE2 functions (andNeighborsInPlace_SSE2).
	   + Restructured debug and error messages so that it was apparent that TComb was
	     responsible.
	   - Removed MMX/ISSE support as times have changed and the support was not going to
	     be carried over to x64.
	   - Removed "opt" parameter. TComb will now use SSE2 if available and will fallback
	     on C++ if it is not supported.


End of tritical version history
------------------------------------------------------------------------
		 

   05/16/2006  v2.0 Beta 2

       + Stricter checking of othreshL/othreshC when looking for oscillation
       + For dot crawl detection require at least one vertical neighbor (y-1/y+1, x-1/x/x+1)
       - fixed possible crash with yuy2 input (sse2 planar<->packed conversions)


   03/31/2006  v2.0 Beta 1

       - complete rewrite


   06/24/2005  v0.9.0

       - Initial Release
	   


contact:    GitHub (@Elegant996)
