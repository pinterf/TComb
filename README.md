# TComb

This is an update to tritical's TComb v2.0 Beta 2 moving it from beta to release as it encompasses all the changes in tritical's To-Do-List.

### Requirements

This filter requires AviSynth 2.6.0 or AviSynth+ as well as the Visual C++ Redistributable Package for Visual Studio 2015-19.

### Syntax and Parameters

The syntax and parameters are identical to the original TComb with the exception of the "opt" parameter. To see a list refer to this [link](http://avisynth.nl/index.php/TComb).

### Changes

In 2015 many changes were made when updating TComb in order to improve speed (see full changelog for more details):

* Removed buffering of frames/info that weren't actually used
* Switched to AVS 2.6 API
* Added x64 support which also utilizes SSE2
* Restructured debug and error messages
* Removed MMX/ISSE support
* Removed/changed "opt" parameter

In 2021 came a general bugfix release by pinterf.
TO-DO: linux port, external assembler to SIMD intrinsics.

### Programmer Notes

This program was compiled using Visual Studio 2019 and falls under the GNU General Public License.

I (Elegant) would like to thank jpsdr and dubhater for their work on nnedi3 and the VapourSynth version of TComb (respectively). Their work led to the port of this project.
I'd also like to thank the masm32 community who were very helpful as I explored assembly.