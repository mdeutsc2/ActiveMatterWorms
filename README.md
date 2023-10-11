# ActiveMatterWorms
Contains all simulation source code for nematic worms project (forked from mvarga6/Nematic-Worms)

## Chapel build instructons
Use chapel-1.32.0 and cuda-11.8.0 (some dependencies required)
CHPL_LLVM=bundled CHPL_CUDA_PATH=/usr/lib/cuda CHPL_LOCALE_MODEL=gpu make -j 6 OPTIMIZE=1 PROFILE=1
