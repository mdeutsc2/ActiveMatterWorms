# ActiveMatterWorms
Contains all simulation source code for nematic worms project (forked from mvarga6/Nematic-Worms)

## Chapel build instructons
Use chapel-1.32.0 and cuda-11.8.0 (some dependencies required)
CHPL_LLVM=bundled CHPL_CUDA_PATH=/usr/lib/cuda CHPL_LOCALE_MODEL=gpu make -j 6 OPTIMIZE=1 PROFILE=1



## Notes
 - For converting images to .avi files that can be read by imageJ for tracking, first run the command:
 ```ffmpeg -framerate 1 -i frame%05d.png -c:v libx264 -r 1 pix_fmt yuv420p anim.avi
 - Then to convert to raw video run
 ``` ffmpeg -i anim.avi -pix_fmt nv12 -f avi -vcodec rawvideo converted.avi