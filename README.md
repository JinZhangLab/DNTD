# Direct non-trilinear decomposition (DNTD)
Source code of DNTD.

**Abstract**： Trilinear decomposition has been employed in the analysis of three-way analytical data. However, it is difficult to keep a perfect trilinear structure for the data in real applications. A direct non-trilinear decomposition (DNTD) algorithm is proposed in this study for analyzing the three-way data with imperfect trilinearity caused by the shift or the variation of the chromatographic or spectral profiles in repeated measurements. In the method, two alternating steps are contained to resolve the constant and shifting profiles for each component, respectively. The first step is based on a pseudo-trilinear decomposition, and the second step involves a slice-based matrix decomposition. Furthermore, for reducing the ambiguity caused by the relaxation of trilinear constraint, an average profile regularization is introduced to enhance the convergence. Three datasets were used to validate the proposed method, including a simulated spectral dataset, a gas chromatography-mass spectrometry (GC-MS) dataset and a flow injection analysis (FIA) dataset. Furthermore, a comparison with parallel factor analysis (PARAFAC), alternating trilinear decomposition (ATLD) and PARAFAC2 was performed. The profiles with variant shifts were correctly extracted and a further understanding of the FIA dataset was obtained.
**Key words**: High-dimension data analysis; Direct non-trilinear decomposition; Parallel factor analysis; Alternating trilinear decomposition; Parallel factor analysis 2; Shift decomposition.


# Use DNTD for a simulated tree-way NIR dataset with imperfect trilinearity
Run test_simul.m script directly in Matlab environment
```matlab
    test_simul.m
```

# If some code is helpful for you, please cites the following reference:
J Zhang, C Guo, WS Cai, XG Shao. Direct non-trilinear decomposition for analyzing high-dimensional data with imperfect trilinearity [J]. Chemometrics and Intelligent Laboratory Systems, 2021.

# File Description
N_way_toolbox:    A matlab toolbox containing useful functions for multi-way data analysis.    
dntd.m:                 Source code of the proposed method （Direct non-trilinear decomposition）.    
plot2way.m:          Function for ploting the resolved profiles by trilinear or non-trilinear decomposition method.    
test_simul.m:         Demonstration of non-trilinear decomposition for the simulated dataset by the proposed method.    
