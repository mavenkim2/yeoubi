A physically based CUDA and OptiX path tracer for USD.

## Features:

### 

### Path Guiding
Path guiding is a variance reduction technique that improves the efficiency of Monte Carlo rendering.
In a standard path tracer, samples are often generated using a combination of BSDF sampling 
and direct light sampling, also known as next event estimation (NEE).
While these techniques work well in simple scenes, they can struggle to sample complex effects 
such as caustics and long indirect illumination. 

Path guiding helps address this issue by learning an approximation of the scene's radiance distribution, which is then sampled 
to guide paths towards directions with higher expected contributions.

Yeoubi's implementation of path guiding has the following components:
- A breadth-first kd-tree builder designed for GPU execution
- A CUDA port of OpenPGL's[^openpgl] mixture fitting code[^ruppert]
- Recording training samples in a wavefront renderer[^wavefront]
- Illumination aware spatial subdivision[^subdivision]
- RIS-based guided sampling[^ris]

**Results**

Enable path guiding when configuring the build:
```bash
cmake -S . -B build -DYBI_ENABLE_PATH_GUIDING=ON
```

### Roadmap:

[^openpgl]: Sebastian Herholz and Addis Dittebrandt. 2022.
    [Intel® Open Path Guiding Library](https://github.com/OpenPathGuidingLibrary/openpgl).

[^wavefront]: Lea Reichardt, Brian Green, Yining Karl Li, and Marco Manzi. 2025.
    [Path Guiding Surfaces and Volumes in Disney's Hyperion Renderer: A Case Study](https://www.yiningkarlli.com/projects/pathguidingcourse2025.html).
    In *ACM SIGGRAPH 2025 Course Notes: Path Guiding in Production and Recent Advancements*. 30–66.

[^ruppert]: Lukas Ruppert, Sebastian Herholz, and Hendrik P. A. Lensch. 2020.
    [Robust Fitting of Parallax-Aware Mixtures for Path Guiding](https://doi.org/10.1145/3386569.3392421).
    *ACM Transactions on Graphics* 39, 4 (Jul. 2020), Article 147.

[^subdivision]: Fengshi Zheng, Christoph Peters, Sebastian Herholz, Marco Manzi, and Elmar Eisemann. 2026.
    [Illumination-Aware Spatial Subdivision for Path Guiding](https://studios.disneyresearch.com/2026/07/01/illumination-aware-spatial-subdivision-for-path-guiding/).
    In *Eurographics Symposium on Rendering (EGSR 2026)*.

[^ris]: Sebastian Herholz. 2025.
    [Integrating Path Guiding into a Production Render: The Nitty Gritty Details](https://sherholz.github.io/siggraph2025-path-guiding-course/).
    In *ACM SIGGRAPH 2025 Course Notes: Path Guiding in Production and Recent Advancements*. 4–29, Section 3.3.
