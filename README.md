Walk-on-Boundary on PBRTv3
--------------------------

This is an extension of PBRT version 3 that implements our SIGGRAPH NA 2023 paper,
[*A Practical Walk-on-Boundary Method for Boundary Value Problems*](https://rsugimoto.net/WoBforBVPsProject/) by Sugimoto, Chen, Jiang, Batty, and Hachisuka.
Specifically, this is the codebase used to generate the figure featured in section 5.3, *WoB within MC rendering*.
This codebase only implements a solver for the interior Dirichlet problem.

## Directory structure

All of the relevant extension code is contained within the `src/wob` subfolder. In particular:

- `wob.cpp` and `wob.h`: Main implementation of the extension
- `scenes`: Example scenes used in rendering. The figure in the paper used the bunny scene.
- `thirdparty`: Dependency libraries

Additional modifications were made to `src/core/api.cpp` to integrate the new `WoBIntegrator` class,
and to `src/core/geometry.h` to fix zero vector normalization.

## Dependencies

As with base PBRT, dependencies are included as git submodules.

    git clone --recurse-submodules https://github.com/tyxchen/WoBPBRT

There is only one external dependency:

- [tinycolormap](https://github.com/yuki-koyama/tinycolormap): a header-only library for colormaps

## Building

```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make pbrt_exe
```

Due to how PBRT is architected, the solution must be rendered separately from the scene;
the separate renders must then be composited to produce the final image.
To render the bunny example:

    python render.py src/wob/scenes/bunny.pbrt --pbrt build/pbrt 

This will produce the following files in `src/wob/scenes`:

- `bunny.png`: the interior solution as rendered by PBRT
- `bunny.png.mask`: a mask file for compositing
- `bunny.png.txt`: the raw numerical solution as rendered by PBRT
- `bunny-scene.png`: the rendered scene
- `bunny-plt.png`: the final composited image

## Credits

PBRT version 3 is released by Matt Pharr, Wenzel Jacob, and Greg Humphreys under the BSD license.

The watertight Stanford Bunny model included in this repository is courtesy of Christopher Batty.
