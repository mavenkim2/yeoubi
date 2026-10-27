A physically based CUDA and OptiX path tracer for USD.

## Features:

### Path Guiding:
Path guiding is a variance reduction technique that improves the efficiency of Monte Carlo rendering.
In a standard path tracer, samples are often generated using a combination of BSDF sampling 
and direct light sampling, also known as next event estimation (NEE).
While these techniques work well in simple scenes, they can struggle to sample complex effects 
such as caustics and long indirect illumination. 
Path guiding helps address this issue by learning an approximation of the scene's radiance distribution, which is then sampled 
to guide paths towards directions with higher expected contributions.

Yeoubi's implementation of path guiding has the following components:
- A breadth-first kd-tree builder amenable to GPU execution
- A CUDA port of OpenPGL's[^openpgl] fitting code for parallax von Mises-Fisher distributions.
- Storing samples in a wavefront renderer.
- Illumination aware spatial subdivision[^subdivision].

Enable path guiding by compiling with the following CMake option:
```bash
cmake -S . -B build -DYBI_ENABLE_PATH_GUIDING=ON
```

YBI_ENABLE_PATH_GUIDING

## References
[^openpgl]: Intel. Open Path Guiding Library (OpenPGL).
[Source code](https://github.com/OpenPathGuidingLibrary/openpgl).

[^ruppert]: Ruppert, Herholz, and Lensch.
*Robust Fitting of Parallax-Aware Mixtures for Path Guiding*.
ACM TOG, 2020.
[Paper](https://doi.org/10.1145/3386569.3392421).

[^subdivision]: *Illumination-Aware Spatial Subdivision for Path Guiding*.
Disney Research, 2026.
[Project page](https://studios.disneyresearch.com/2026/07/01/illumination-aware-spatial-subdivision-for-path-guiding/).
