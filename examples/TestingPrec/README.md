# Testing of Regularization and BFGS preconditioning

This is the place to test the effect different regularizations 


### H2, regularization_w1=1e1, regularization_length_scale = None

    Inversion.LBFGS: ********** iteration  11 **********
    Inversion.LBFGS:        F(m) = 3.64408
    Inversion.ERT-H2: search direction component 0 = Summary: inf=-1.56269e-05 sup=2.01867e-05 data points=14674.
    Inversion.ERT-H2: search direction component 1 = Summary: inf=-1.01924e-05 sup=1.02442e-05 data points=14674.
    Inversion.ERT-H2: search direction component 2 = Summary: inf=-5.48148e-06 sup=8.21463e-06 data points=14674.
    Inversion.LBFGS: Search direction scaled : p=5.37963e-06 (scale factor = 0.266493)
    Inversion.LBFGS.linesearch: Setting options: {'alphaMin': 2.136932660031764e-06}
    Inversion.LBFGS: Starting line search with alphaMin, alpha  = 2.13693e-06, 0.408414
    Inversion.LBFGS.linesearch: Initial values: phi(0)=3.64408, phi'(0)=-0.000110374, alpha= 0.408414.
    Inversion.ERT-H2: M = Summary: inf=-0.0329614 sup=0.0198522 data points=14674
    Inversion.ERT-H2: m = Summary: inf=-0.159087 sup=0.254463 data points=315556
    Inversion.ERT-H2: sigma_0 = Summary: inf=0.000485731 sup=0.000734511 data points=315556
    Inversion.ERT-H2: sigma_0_face = Summary: inf=0.00056949 sup=0.00056949 data points=2742
    Inversion.ERT-H2: sigma_0_stations = 0.0004915589556201047 - 0.0005418927961528474 
    Inversion.ERT-H2: 8 additive DC potentials calculated.
    Inversion.ERT-H2: misfit ERT, reg, curl - total         =  3.167377e+00, 4.766704e-01, 0.000000e+00 -  3.644047e+00
    Inversion.ERT-H2: ratios ERT, reg, curl  [%]    =  86.9192, 13.0808, 0
    Inversion.LBFGS.linesearch: Iteration 1, alpha=0.408414, phi(alpha)=3.64405
    Inversion.ERT-H2: 8 adjoint potentials calculated.
    Inversion.LBFGS.linesearch: phi'(alpha)=-6.53397e-05
    Inversion.LBFGS.linesearch: Strong Wolfe condition is fulfilled. We are done.
    Inversion.LBFGS.linesearch: Line search completed after 1 steps (alpha=0.408414, phi(alpha)=3.64405).
    Inversion.LBFGS: Search direction scaling found as alpha=0.408414
    Inversion.LBFGS: Solution checked: |m-m_old|=6.29465e-05, |m|*m_tol=3.29614e-06
    Inversion.LBFGS: F(m) = 3.64405
    Inversion.LBFGS: Gradient has converged: |F-Fold|=3.60476e-05 < g_tol*max(|F|,|Fold|)=3.64408e-05
    Inversion.LBFGS: Success after 11 iterations!

### H2, regularization_w1=1e1, regularization_length_scale = None
    Inversion.LBFGS: ********** iteration  21 **********
    Inversion.LBFGS:        F(m) = 1.84977
    Inversion.ERT-H2: search direction component 0 = Summary: inf=-0.000374912 sup=0.000576077 data points=14674.
    Inversion.ERT-H2: search direction component 1 = Summary: inf=-0.000211613 sup=0.000213644 data points=14674.
    Inversion.ERT-H2: search direction component 2 = Summary: inf=-0.00056838 sup=0 data points=14674.
    Inversion.LBFGS: Search direction scaled : p=0.000587238 (scale factor = 1.01937)
    Inversion.LBFGS.linesearch: Setting options: {'alphaMin': 2.1321835502483634e-07}
    Inversion.LBFGS: Starting line search with alphaMin, alpha  = 2.13218e-07, 0.229212
    Inversion.LBFGS.linesearch: Initial values: phi(0)=1.84977, phi'(0)=-0.000336269, alpha= 0.229212.
    Inversion.ERT-H2: M = Summary: inf=-0.148494 sup=0.0945324 data points=14674
    Inversion.ERT-H2: m = Summary: inf=-0.679288 sup=1.11518 data points=315556
    Inversion.ERT-H2: sigma_0 = Summary: inf=0.000288719 sup=0.00173701 data points=315556
    Inversion.ERT-H2: sigma_0_face = Summary: inf=0.00056949 sup=0.00056949 data points=2742
    Inversion.ERT-H2: sigma_0_stations = 0.0003025157536409129 - 0.000478273363845199 
    Inversion.ERT-H2: 8 additive DC potentials calculated.
    Inversion.ERT-H2: misfit ERT, reg, curl - total         =  9.273667e-01, 9.251450e-01, 0.000000e+00 -  1.852512e+00
    Inversion.ERT-H2: ratios ERT, reg, curl  [%]    =  50.06, 49.94, 0
    Inversion.LBFGS.linesearch: Iteration 1, alpha=0.229212, phi(alpha)=1.85251
    Inversion.LBFGS.linesearch: Sufficient decrease condition or decrease criterium are violated -> start zoom.
    Inversion.LBFGS.linesearch.zoom: Iteration 1: alpha range =[0, 0.229212] (width= 0.229212)
    Inversion.LBFGS.linesearch.zoom: Interpolation nodes: [0.22921176861410367, 0.0]
    Inversion.LBFGS.linesearch.zoom: Quadratic interpolation is applied.
    Inversion.LBFGS.linesearch.zoom: Insufficient search interval reduction.
    Inversion.LBFGS.linesearch.zoom: Bisection is applied.
    Inversion.ERT-H2: M = Summary: inf=-0.147896 sup=0.0939751 data points=14674
    Inversion.ERT-H2: m = Summary: inf=-0.67668 sup=1.10935 data points=315556
    Inversion.ERT-H2: sigma_0 = Summary: inf=0.000289473 sup=0.00172691 data points=315556
    Inversion.ERT-H2: sigma_0_face = Summary: inf=0.00056949 sup=0.00056949 data points=2742
    Inversion.ERT-H2: sigma_0_stations = 0.0003034939112334561 - 0.0004767107581247185 
    Inversion.ERT-H2: 8 additive DC potentials calculated.
    Inversion.ERT-H2: misfit ERT, reg, curl - total         =  9.346304e-01, 9.158102e-01, 0.000000e+00 -  1.850441e+00
    Inversion.ERT-H2: ratios ERT, reg, curl  [%]    =  50.5085, 49.4915, 0
    Inversion.LBFGS.linesearch.zoom: Iteration 1, alpha=0.114606, phi(alpha)=1.85044 (interpolation order = 1)
    Inversion.LBFGS.linesearch.zoom: Iteration 1: alpha_hi updated, now = 0.114606
    Inversion.LBFGS.linesearch.zoom: Iteration 2: alpha range =[0, 0.114606] (width= 0.114606)
    Inversion.LBFGS.linesearch.zoom: Interpolation nodes: [0.11460599091622935, 0.22921176861410367, 0.0]
    Inversion.LBFGS.linesearch.zoom: Cubic interpolation is applied.
    Inversion.LBFGS.linesearch.zoom: Insufficient search interval reduction.
    Inversion.LBFGS.linesearch.zoom: Quadratic interpolation is applied.
    Inversion.LBFGS.linesearch.zoom: Insufficient search interval reduction.
    Inversion.LBFGS.linesearch.zoom: Bisection is applied.
    Inversion.ERT-H2: M = Summary: inf=-0.147598 sup=0.0936964 data points=14674
    Inversion.ERT-H2: m = Summary: inf=-0.67542 sup=1.10643 data points=315556
    Inversion.ERT-H2: sigma_0 = Summary: inf=0.000289838 sup=0.00172188 data points=315556
    Inversion.ERT-H2: sigma_0_face = Summary: inf=0.00056949 sup=0.00056949 data points=2742
    Inversion.ERT-H2: sigma_0_stations = 0.00030398417543570105 - 0.00047593137080174404 
    Inversion.ERT-H2: 8 additive DC potentials calculated.
    Inversion.ERT-H2: misfit ERT, reg, curl - total         =  9.387639e-01, 9.111674e-01, 0.000000e+00 -  1.849931e+00
    Inversion.ERT-H2: ratios ERT, reg, curl  [%]    =  50.7459, 49.2541, 0
    Inversion.LBFGS.linesearch.zoom: Iteration 2, alpha=0.0573031, phi(alpha)=1.84993 (interpolation order = 1)
    Inversion.LBFGS.linesearch.zoom: Iteration 2: alpha_hi updated, now = 0.0573031
    Inversion.LBFGS.linesearch.zoom: Iteration 3: alpha range =[0, 0.0573031] (width= 0.0573031)
    Inversion.LBFGS.linesearch.zoom: Interpolation nodes: [0.05730310206729219, 0.11460599091622935, 0.0]
    Inversion.LBFGS.linesearch.zoom: Cubic interpolation is applied.
    Inversion.LBFGS.linesearch.zoom: Insufficient search interval reduction.
    Inversion.LBFGS.linesearch.zoom: Quadratic interpolation is applied.
    Inversion.LBFGS.linesearch.zoom: Insufficient search interval reduction.
    Inversion.LBFGS.linesearch.zoom: Bisection is applied.
    Inversion.ERT-H2: M = Summary: inf=-0.147448 sup=0.0935571 data points=14674
    Inversion.ERT-H2: m = Summary: inf=-0.67479 sup=1.10497 data points=315556
    Inversion.ERT-H2: sigma_0 = Summary: inf=0.000290021 sup=0.00171938 data points=315556
    Inversion.ERT-H2: sigma_0_face = Summary: inf=0.00056949 sup=0.00056949 data points=2742
    Inversion.ERT-H2: sigma_0_stations = 0.0003042296044467934 - 0.0004755421551110801 
    Inversion.ERT-H2: 8 additive DC potentials calculated.
    Inversion.ERT-H2: misfit ERT, reg, curl - total         =  9.409567e-01, 9.088522e-01, 0.000000e+00 -  1.849809e+00
    Inversion.ERT-H2: ratios ERT, reg, curl  [%]    =  50.8678, 49.1322, 0
    Inversion.LBFGS.linesearch.zoom: Iteration 3, alpha=0.0286517, phi(alpha)=1.84981 (interpolation order = 1)
    Inversion.LBFGS.linesearch.zoom: Iteration 3: alpha_hi updated, now = 0.0286517
    Inversion.LBFGS.linesearch.zoom: Iteration 4: alpha range =[0, 0.0286517] (width= 0.0286517)
    Inversion.LBFGS.linesearch.zoom: Interpolation nodes: [0.028651657642823607, 0.05730310206729219, 0.0]
    Inversion.LBFGS.linesearch.zoom: Cubic interpolation is applied.
    Inversion.LBFGS.linesearch.zoom: Insufficient search interval reduction.
    Inversion.LBFGS.linesearch.zoom: Quadratic interpolation is applied.
    Inversion.LBFGS.linesearch.zoom: Insufficient search interval reduction.
    Inversion.LBFGS.linesearch.zoom: Bisection is applied.
    Inversion.ERT-H2: M = Summary: inf=-0.147374 sup=0.0934874 data points=14674
    Inversion.ERT-H2: m = Summary: inf=-0.674474 sup=1.10425 data points=315556
    Inversion.ERT-H2: sigma_0 = Summary: inf=0.000290112 sup=0.00171812 data points=315556
    Inversion.ERT-H2: sigma_0_face = Summary: inf=0.00056949 sup=0.00056949 data points=2742
    Inversion.ERT-H2: sigma_0_stations = 0.00030435239324974714 - 0.00047534766664444195 
    Inversion.ERT-H2: 8 additive DC potentials calculated.
    Inversion.ERT-H2: misfit ERT, reg, curl - total         =  9.420847e-01, 9.076961e-01, 0.000000e+00 -  1.849781e+00
    Inversion.ERT-H2: ratios ERT, reg, curl  [%]    =  50.9295, 49.0705, 0
    Inversion.LBFGS.linesearch.zoom: Iteration 4, alpha=0.0143259, phi(alpha)=1.84978 (interpolation order = 1)
    Inversion.LBFGS.linesearch.zoom: Iteration 4: alpha_hi updated, now = 0.0143259
    Inversion.LBFGS.linesearch.zoom: Iteration 5: alpha range =[0, 0.0143259] (width= 0.0143259)
    Inversion.LBFGS.linesearch.zoom: Interpolation nodes: [0.014325935430589316, 0.028651657642823607, 0.0]
    Inversion.LBFGS.linesearch.zoom: Cubic interpolation is applied.
    Inversion.LBFGS.linesearch.zoom: Insufficient search interval reduction.
    Inversion.LBFGS.linesearch.zoom: Quadratic interpolation is applied.
    Inversion.LBFGS.linesearch.zoom: Insufficient search interval reduction.
    Inversion.LBFGS.linesearch.zoom: Bisection is applied.
    Inversion.ERT-H2: M = Summary: inf=-0.147337 sup=0.0934526 data points=14674
    Inversion.ERT-H2: m = Summary: inf=-0.674317 sup=1.10388 data points=315556
    Inversion.ERT-H2: sigma_0 = Summary: inf=0.000290158 sup=0.0017175 data points=315556
    Inversion.ERT-H2: sigma_0_face = Summary: inf=0.00056949 sup=0.00056949 data points=2742
    Inversion.ERT-H2: sigma_0_stations = 0.00030441380623432167 - 0.0004752504522415557 
    Inversion.ERT-H2: 8 additive DC potentials calculated.
    Inversion.ERT-H2: misfit ERT, reg, curl - total         =  9.426566e-01, 9.071185e-01, 0.000000e+00 -  1.849775e+00
    Inversion.ERT-H2: ratios ERT, reg, curl  [%]    =  50.9606, 49.0394, 0
    Inversion.LBFGS.linesearch.zoom: Iteration 5, alpha=0.00716307, phi(alpha)=1.84978 (interpolation order = 1)
    Inversion.LBFGS.linesearch.zoom: Iteration 5: alpha_hi updated, now = 0.00716307
    Inversion.LBFGS.linesearch.zoom: Iteration 6: alpha range =[0, 0.00716307] (width= 0.00716307)
    Inversion.LBFGS.linesearch.zoom: Interpolation nodes: [0.007163074324472171, 0.014325935430589316, 0.0]
    Inversion.LBFGS.linesearch.zoom: Cubic interpolation is applied.
    Inversion.ERT-H2: M = Summary: inf=-0.147316 sup=0.0934339 data points=14674
    Inversion.ERT-H2: m = Summary: inf=-0.674232 sup=1.10369 data points=315556
    Inversion.ERT-H2: sigma_0 = Summary: inf=0.000290182 sup=0.00171716 data points=315556
    Inversion.ERT-H2: sigma_0_face = Summary: inf=0.00056949 sup=0.00056949 data points=2742
    Inversion.ERT-H2: sigma_0_stations = 0.0003044467268715124 - 0.0004751983563914578 
    Inversion.ERT-H2: 8 additive DC potentials calculated.
    Inversion.ERT-H2: misfit ERT, reg, curl - total         =  9.429653e-01, 9.068090e-01, 0.000000e+00 -  1.849774e+00
    Inversion.ERT-H2: ratios ERT, reg, curl  [%]    =  50.9773, 49.0227, 0
    Inversion.LBFGS.linesearch.zoom: Iteration 6, alpha=0.00332399, phi(alpha)=1.84977 (interpolation order = 3)
    Inversion.ERT-H2: 8 adjoint potentials calculated.
    Inversion.LBFGS.linesearch.zoom: phi'(alpha)=2.22295e-05
    Inversion.LBFGS.linesearch.zoom: Zoom completed after 6 steps: alpha=0.00332399, phi(alpha)=1.84977, phi'(alpha)=2.22295e-05.
    Inversion.LBFGS.linesearch.zoom: Strong Wolfe condition is fulfilled. We are done.
    Inversion.LBFGS.linesearch: Line search completed after 1 steps (alpha=0.00332399, phi(alpha)=1.84977).
    Inversion.LBFGS: Search direction scaling found as alpha=0.00332399
    Inversion.LBFGS: Solution checked: |m-m_old|=2.29634e-05, |m|*m_tol=1.47316e-05
    Inversion.LBFGS: F(m) = 1.84977
    Inversion.LBFGS: Gradient has converged: |F-Fold|=5.91035e-07 < g_tol*max(|F|,|Fold|)=1.84977e-05
    Inversion.LBFGS: Success after 21 iterations!

### H2, regularization_w1=1e-1, regularization_length_scale = None
    Inversion.LBFGS: ********** iteration  42 **********
    Inversion.LBFGS:        F(m) = 0.388102
    Inversion.ERT-H2: search direction component 0 = Summary: inf=-0.00661094 sup=0.00721155 data points=14674.
    Inversion.ERT-H2: search direction component 1 = Summary: inf=-0.00241269 sup=0.00252678 data points=14674.
    Inversion.ERT-H2: search direction component 2 = Summary: inf=-0.00318268 sup=0.00390604 data points=14674.
    Inversion.LBFGS: Search direction scaled : p=0.000617117 (scale factor = 0.0855734)
    Inversion.LBFGS.linesearch: Setting options: {'alphaMin': 1.7195489472745922e-06}
    Inversion.LBFGS: Starting line search with alphaMin, alpha  = 1.71955e-06, 0.100933
    Inversion.LBFGS.linesearch: Initial values: phi(0)=0.388102, phi'(0)=-5.79004e-05, alpha= 0.100933.
    Inversion.ERT-H2: M = Summary: inf=-0.290292 sup=0.165606 data points=14674
    Inversion.ERT-H2: m = Summary: inf=-1.25288 sup=2.16263 data points=315556
    Inversion.ERT-H2: sigma_0 = Summary: inf=0.000162693 sup=0.00495113 data points=315556
    Inversion.ERT-H2: sigma_0_face = Summary: inf=0.00056949 sup=0.00056949 data points=2742
    Inversion.ERT-H2: sigma_0_stations = 0.00017875588477634426 - 0.0003757215898754952 
    Inversion.ERT-H2: 8 additive DC potentials calculated.
    Inversion.ERT-H2: misfit ERT, reg, curl - total         =  6.394112e-02, 3.242382e-01, 0.000000e+00 -  3.881793e-01
    Inversion.ERT-H2: ratios ERT, reg, curl  [%]    =  16.4721, 83.5279, 0
    Inversion.LBFGS.linesearch: Iteration 1, alpha=0.100933, phi(alpha)=0.388179
    Inversion.LBFGS.linesearch: Sufficient decrease condition or decrease criterium are violated -> start zoom.
    Inversion.LBFGS.linesearch.zoom: Iteration 1: alpha range =[0, 0.100933] (width= 0.100933)
    Inversion.LBFGS.linesearch.zoom: Interpolation nodes: [0.10093268377142908, 0.0]
    Inversion.LBFGS.linesearch.zoom: Quadratic interpolation is applied.
    Inversion.LBFGS.linesearch.zoom: Insufficient search interval reduction.
    Inversion.LBFGS.linesearch.zoom: Bisection is applied.
    Inversion.ERT-H2: M = Summary: inf=-0.290289 sup=0.165598 data points=14674
    Inversion.ERT-H2: m = Summary: inf=-1.25207 sup=2.16243 data points=315556
    Inversion.ERT-H2: sigma_0 = Summary: inf=0.000162824 sup=0.00495013 data points=315556
    Inversion.ERT-H2: sigma_0_face = Summary: inf=0.00056949 sup=0.00056949 data points=2742
    Inversion.ERT-H2: sigma_0_stations = 0.00017888168031570482 - 0.0003758747238180299 
    Inversion.ERT-H2: 8 additive DC potentials calculated.
    Inversion.ERT-H2: misfit ERT, reg, curl - total         =  6.402725e-02, 3.240924e-01, 0.000000e+00 -  3.881196e-01
    Inversion.ERT-H2: ratios ERT, reg, curl  [%]    =  16.4968, 83.5032, 0
    Inversion.LBFGS.linesearch.zoom: Iteration 1, alpha=0.0504672, phi(alpha)=0.38812 (interpolation order = 1)
    Inversion.LBFGS.linesearch.zoom: Iteration 1: alpha_hi updated, now = 0.0504672
    Inversion.LBFGS.linesearch.zoom: Iteration 2: alpha range =[0, 0.0504672] (width= 0.0504672)
    Inversion.LBFGS.linesearch.zoom: Interpolation nodes: [0.05046720166018818, 0.10093268377142908, 0.0]
    Inversion.LBFGS.linesearch.zoom: Cubic interpolation is applied.
    Inversion.LBFGS.linesearch.zoom: Insufficient search interval reduction.
    Inversion.LBFGS.linesearch.zoom: Quadratic interpolation is applied.
    Inversion.LBFGS.linesearch.zoom: Insufficient search interval reduction.
    Inversion.LBFGS.linesearch.zoom: Bisection is applied.
    Inversion.ERT-H2: M = Summary: inf=-0.290288 sup=0.165593 data points=14674
    Inversion.ERT-H2: m = Summary: inf=-1.25167 sup=2.16233 data points=315556
    Inversion.ERT-H2: sigma_0 = Summary: inf=0.000162889 sup=0.00494964 data points=315556
    Inversion.ERT-H2: sigma_0_face = Summary: inf=0.00056949 sup=0.00056949 data points=2742
    Inversion.ERT-H2: sigma_0_stations = 0.0001789446112786883 - 0.0003759513141926746 
    Inversion.ERT-H2: 8 additive DC potentials calculated.
    Inversion.ERT-H2: misfit ERT, reg, curl - total         =  6.408594e-02, 3.240195e-01, 0.000000e+00 -  3.881054e-01
    Inversion.ERT-H2: ratios ERT, reg, curl  [%]    =  16.5125, 83.4875, 0
    Inversion.LBFGS.linesearch.zoom: Iteration 2, alpha=0.0252345, phi(alpha)=0.388105 (interpolation order = 1)
    Inversion.LBFGS.linesearch.zoom: Iteration 2: alpha_hi updated, now = 0.0252345
    Inversion.LBFGS.linesearch.zoom: Iteration 3: alpha range =[0, 0.0252345] (width= 0.0252345)
    Inversion.LBFGS.linesearch.zoom: Interpolation nodes: [0.025234460604567726, 0.05046720166018818, 0.0]
    Inversion.LBFGS.linesearch.zoom: Cubic interpolation is applied.
    Inversion.LBFGS.linesearch.zoom: Insufficient search interval reduction.
    Inversion.LBFGS.linesearch.zoom: Quadratic interpolation is applied.
    Inversion.LBFGS.linesearch.zoom: Insufficient search interval reduction.
    Inversion.LBFGS.linesearch.zoom: Bisection is applied.
    Inversion.ERT-H2: M = Summary: inf=-0.290288 sup=0.165591 data points=14674
    Inversion.ERT-H2: m = Summary: inf=-1.25147 sup=2.16228 data points=315556
    Inversion.ERT-H2: sigma_0 = Summary: inf=0.000162922 sup=0.00494939 data points=315556
    Inversion.ERT-H2: sigma_0_face = Summary: inf=0.00056949 sup=0.00056949 data points=2742
    Inversion.ERT-H2: sigma_0_stations = 0.00017897608506191193 - 0.0003759896152322331 
    Inversion.ERT-H2: 8 additive DC potentials calculated.
    Inversion.ERT-H2: misfit ERT, reg, curl - total         =  6.411919e-02, 3.239830e-01, 0.000000e+00 -  3.881022e-01
    Inversion.ERT-H2: ratios ERT, reg, curl  [%]    =  16.5212, 83.4788, 0
    Inversion.LBFGS.linesearch.zoom: Iteration 3, alpha=0.0126181, phi(alpha)=0.388102 (interpolation order = 1)
    Inversion.LBFGS.linesearch.zoom: Iteration 3: alpha_hi updated, now = 0.0126181
    Inversion.LBFGS.linesearch.zoom: Iteration 4: alpha range =[0, 0.0126181] (width= 0.0126181)
    Inversion.LBFGS.linesearch.zoom: Interpolation nodes: [0.0126180900767575, 0.025234460604567726, 0.0]
    Inversion.LBFGS.linesearch.zoom: Cubic interpolation is applied.
    Inversion.LBFGS.linesearch.zoom: Insufficient search interval reduction.
    Inversion.LBFGS.linesearch.zoom: Quadratic interpolation is applied.
    Inversion.LBFGS.linesearch.zoom: Insufficient search interval reduction.
    Inversion.LBFGS.linesearch.zoom: Bisection is applied.
    Inversion.ERT-H2: M = Summary: inf=-0.290287 sup=0.16559 data points=14674
    Inversion.ERT-H2: m = Summary: inf=-1.25137 sup=2.16225 data points=315556
    Inversion.ERT-H2: sigma_0 = Summary: inf=0.000162938 sup=0.00494926 data points=315556
    Inversion.ERT-H2: sigma_0_face = Summary: inf=0.00056949 sup=0.00056949 data points=2742
    Inversion.ERT-H2: sigma_0_stations = 0.00017899182402938263 - 0.00037600876721524505 
    Inversion.ERT-H2: 8 additive DC potentials calculated.
    Inversion.ERT-H2: misfit ERT, reg, curl - total         =  6.413678e-02, 3.239648e-01, 0.000000e+00 -  3.881016e-01
    Inversion.ERT-H2: ratios ERT, reg, curl  [%]    =  16.5258, 83.4742, 0
    Inversion.LBFGS.linesearch.zoom: Iteration 4, alpha=0.0063099, phi(alpha)=0.388102 (interpolation order = 1)
    Inversion.ERT-H2: 8 adjoint potentials calculated.
    Inversion.LBFGS.linesearch.zoom: phi'(alpha)=4.54537e-05
    Inversion.LBFGS.linesearch.zoom: Zoom completed after 4 steps: alpha=0.0063099, phi(alpha)=0.388102, phi'(alpha)=4.54537e-05.
    Inversion.LBFGS.linesearch.zoom: Strong Wolfe condition is fulfilled. We are done.
    Inversion.LBFGS.linesearch: Line search completed after 1 steps (alpha=0.0063099, phi(alpha)=0.388102).
    Inversion.LBFGS: Search direction scaling found as alpha=0.0063099
    Inversion.LBFGS: F(m) = 0.388102
    Inversion.LBFGS: Solution has converged: |m-m_old|=1.06521e-05 < |m|*m_tol=2.90287e-05
    Inversion.LBFGS: Success after 42 iterations!

### PseudoGauss, regularization_w1=1e-3, regularization_length_scale = 5
    Inversion.LBFGS: ********** iteration  16 **********
    Inversion.LBFGS:        F(m) = 1.13251
    Inversion.ERT-Gauss: search direction component 0 = Summary: inf=-0.0015567 sup=0.00219379 data points=14674.
    Inversion.ERT-Gauss: search direction component 1 = Summary: inf=-0.00119268 sup=0.00060357 data points=14674.
    Inversion.ERT-Gauss: search direction component 2 = Summary: inf=-0.000774598 sup=0.000520821 data points=14674.
    Inversion.ERT-Gauss: search direction component 3 = Summary: inf=-0.00154158 sup=0.00192125 data points=14674.
    Inversion.LBFGS: Search direction scaled : p=0.00178228 (scale factor = 0.812419)
    Inversion.LBFGS.linesearch: Setting options: {'alphaMin': 3.2916804327598444e-06}
    Inversion.LBFGS: Starting line search with alphaMin, alpha  = 3.29168e-06, 0.487156
    Inversion.LBFGS.linesearch: Initial values: phi(0)=1.13251, phi'(0)=-1.39453e-05, alpha= 0.487156.
    Inversion.ERT-Gauss: m = Summary: inf=-1.63258 sup=0.784265 data points=315556
    Inversion.ERT-Gauss: sigma_0 = Summary: inf=0.000111292 sup=0.00124764 data points=315556
    Inversion.ERT-Gauss: sigma_0_face = Summary: inf=0.00056949 sup=0.00056949 data points=2742
    Inversion.ERT-Gauss: sigma_0_stations = 0.0001853517076381842 - 0.0007710954367138727 
    Inversion.ERT-Gauss: 8 additive DC potentials calculated.
    Inversion.ERT-Gauss: misfit ERT, reg, div, curl; total  =  5.233104e-01, 3.025722e-01, 6.474094e-01, 2.684098e-01; 1.132506e+00
    Inversion.ERT-Gauss: ratios ERT, reg  [%]       =  46.2082, 53.7918
    Inversion.LBFGS.linesearch: Iteration 1, alpha=0.487156, phi(alpha)=1.13251
    Inversion.ERT-Gauss: 8 adjoint potentials calculated.
    Inversion.LBFGS.linesearch: phi'(alpha)=-6.20139e-06
    Inversion.LBFGS.linesearch: Strong Wolfe condition is fulfilled. We are done.
    Inversion.LBFGS.linesearch: Line search completed after 1 steps (alpha=0.487156, phi(alpha)=1.13251).
    Inversion.LBFGS: Search direction scaling found as alpha=0.487156
    Inversion.LBFGS: Solution checked: |m-m_old|=0.00342209, |m|*m_tol=0.000231062
    Inversion.LBFGS: F(m) = 1.13251
    Inversion.LBFGS: Gradient has converged: |F-Fold|=4.95567e-06 < g_tol*max(|F|,|Fold|)=1.13251e-05
    Inversion.LBFGS: Success after 16 iterations!

### PseudoGauss, regularization_w1=1e-4, regularization_length_scale = 5
    Inversion.LBFGS: ********** iteration  21 **********
    Inversion.LBFGS:        F(m) = 0.214754
    Inversion.ERT-Gauss: search direction component 0 = Summary: inf=-0.0981558 sup=0.00805389 data points=14674.
    Inversion.ERT-Gauss: search direction component 1 = Summary: inf=-0.021271 sup=0.0498985 data points=14674.
    Inversion.ERT-Gauss: search direction component 2 = Summary: inf=-0.02351 sup=0.0188912 data points=14674.
    Inversion.ERT-Gauss: search direction component 3 = Summary: inf=-0.103768 sup=0.00656513 data points=14674.
    Inversion.LBFGS: Search direction scaled : p=0.0289573 (scale factor = 0.279058)
    Inversion.LBFGS.linesearch: Setting options: {'alphaMin': 9.440527570718719e-07}
    Inversion.LBFGS: Starting line search with alphaMin, alpha  = 9.44053e-07, 0.388771
    Inversion.LBFGS.linesearch: Initial values: phi(0)=0.214754, phi'(0)=-0.000191937, alpha= 0.388771.
    Inversion.ERT-Gauss: m = Summary: inf=-2.60068 sup=1.44855 data points=315556
    Inversion.ERT-Gauss: sigma_0 = Summary: inf=4.22691e-05 sup=0.00242428 data points=315556
    Inversion.ERT-Gauss: sigma_0_face = Summary: inf=0.00056949 sup=0.00056949 data points=2742
    Inversion.ERT-Gauss: sigma_0_stations = 9.031042016492083e-05 - 0.0012725735949518765 
    Inversion.ERT-Gauss: 8 additive DC potentials calculated.
    Inversion.ERT-Gauss: misfit ERT, reg, div, curl; total  =  2.602809e-02, 9.713133e-02, 2.012707e-01, 8.078286e-02; 2.156205e-01
    Inversion.ERT-Gauss: ratios ERT, reg  [%]       =  12.0712, 87.9288
    Inversion.LBFGS.linesearch: Iteration 1, alpha=0.388771, phi(alpha)=0.215621
    Inversion.LBFGS.linesearch: Sufficient decrease condition or decrease criterium are violated -> start zoom.
    Inversion.LBFGS.linesearch.zoom: Iteration 1: alpha range =[0, 0.388771] (width= 0.388771)
    Inversion.LBFGS.linesearch.zoom: Interpolation nodes: [0.38877078562640727, 0.0]
    Inversion.LBFGS.linesearch.zoom: Quadratic interpolation is applied.
    Inversion.LBFGS.linesearch.zoom: Insufficient search interval reduction.
    Inversion.LBFGS.linesearch.zoom: Bisection is applied.
    Inversion.ERT-Gauss: m = Summary: inf=-2.60491 sup=1.44346 data points=315556
    Inversion.ERT-Gauss: sigma_0 = Summary: inf=4.2091e-05 sup=0.00241199 data points=315556
    Inversion.ERT-Gauss: sigma_0_face = Summary: inf=0.00056949 sup=0.00056949 data points=2742
    Inversion.ERT-Gauss: sigma_0_stations = 8.997573584392243e-05 - 0.001267025250888924 
    Inversion.ERT-Gauss: 8 additive DC potentials calculated.
    Inversion.ERT-Gauss: misfit ERT, reg, div, curl; total  =  2.573047e-02, 9.693962e-02, 2.007944e-01, 8.070957e-02; 2.149523e-01
    Inversion.ERT-Gauss: ratios ERT, reg  [%]       =  11.9703, 88.0297
    Inversion.LBFGS.linesearch.zoom: Iteration 1, alpha=0.194386, phi(alpha)=0.214952 (interpolation order = 1)
    Inversion.LBFGS.linesearch.zoom: Iteration 1: alpha_hi updated, now = 0.194386
    Inversion.LBFGS.linesearch.zoom: Iteration 2: alpha range =[0, 0.194386] (width= 0.194386)
    Inversion.LBFGS.linesearch.zoom: Interpolation nodes: [0.19438586483958217, 0.38877078562640727, 0.0]
    Inversion.LBFGS.linesearch.zoom: Cubic interpolation is applied.
    Inversion.LBFGS.linesearch.zoom: Insufficient search interval reduction.
    Inversion.LBFGS.linesearch.zoom: Quadratic interpolation is applied.
    Inversion.LBFGS.linesearch.zoom: Insufficient search interval reduction.
    Inversion.LBFGS.linesearch.zoom: Bisection is applied.
    Inversion.ERT-Gauss: m = Summary: inf=-2.60702 sup=1.44092 data points=315556
    Inversion.ERT-Gauss: sigma_0 = Summary: inf=4.20022e-05 sup=0.00240587 data points=315556
    Inversion.ERT-Gauss: sigma_0_face = Summary: inf=0.00056949 sup=0.00056949 data points=2742
    Inversion.ERT-Gauss: sigma_0_stations = 8.980885909016971e-05 - 0.0012642601568770124 
    Inversion.ERT-Gauss: 8 additive DC potentials calculated.
    Inversion.ERT-Gauss: misfit ERT, reg, div, curl; total  =  2.575663e-02, 9.684423e-02, 2.005574e-01, 8.067329e-02; 2.147941e-01
    Inversion.ERT-Gauss: ratios ERT, reg  [%]       =  11.9913, 88.0087
    Inversion.LBFGS.linesearch.zoom: Iteration 2, alpha=0.0971934, phi(alpha)=0.214794 (interpolation order = 1)
    Inversion.LBFGS.linesearch.zoom: Iteration 2: alpha_hi updated, now = 0.0971934
    Inversion.LBFGS.linesearch.zoom: Iteration 3: alpha range =[0, 0.0971934] (width= 0.0971934)
    Inversion.LBFGS.linesearch.zoom: Interpolation nodes: [0.09719340444616963, 0.19438586483958217, 0.0]
    Inversion.LBFGS.linesearch.zoom: Cubic interpolation is applied.
    Inversion.LBFGS.linesearch.zoom: Insufficient search interval reduction.
    Inversion.LBFGS.linesearch.zoom: Quadratic interpolation is applied.
    Inversion.LBFGS.linesearch.zoom: Insufficient search interval reduction.
    Inversion.LBFGS.linesearch.zoom: Bisection is applied.
    Inversion.ERT-Gauss: m = Summary: inf=-2.60808 sup=1.43965 data points=315556
    Inversion.ERT-Gauss: sigma_0 = Summary: inf=4.19579e-05 sup=0.00240281 data points=315556
    Inversion.ERT-Gauss: sigma_0_face = Summary: inf=0.00056949 sup=0.00056949 data points=2742
    Inversion.ERT-Gauss: sigma_0_stations = 8.972553681319868e-05 - 0.0012628798735970905 
    Inversion.ERT-Gauss: 8 additive DC potentials calculated.
    Inversion.ERT-Gauss: misfit ERT, reg, div, curl; total  =  2.581362e-02, 9.679665e-02, 2.004393e-01, 8.065524e-02; 2.147592e-01
    Inversion.ERT-Gauss: ratios ERT, reg  [%]       =  12.0198, 87.9802
    Inversion.LBFGS.linesearch.zoom: Iteration 3, alpha=0.0485972, phi(alpha)=0.214759 (interpolation order = 1)
    Inversion.LBFGS.linesearch.zoom: Iteration 3: alpha_hi updated, now = 0.0485972
    Inversion.LBFGS.linesearch.zoom: Iteration 4: alpha range =[0, 0.0485972] (width= 0.0485972)
    Inversion.LBFGS.linesearch.zoom: Interpolation nodes: [0.048597174249463346, 0.09719340444616963, 0.0]
    Inversion.LBFGS.linesearch.zoom: Cubic interpolation is applied.
    Inversion.LBFGS.linesearch.zoom: Insufficient search interval reduction.
    Inversion.LBFGS.linesearch.zoom: Quadratic interpolation is applied.
    Inversion.LBFGS.linesearch.zoom: Insufficient search interval reduction.
    Inversion.LBFGS.linesearch.zoom: Bisection is applied.
    Inversion.ERT-Gauss: m = Summary: inf=-2.6086 sup=1.43902 data points=315556
    Inversion.ERT-Gauss: sigma_0 = Summary: inf=4.19357e-05 sup=0.00240129 data points=315556
    Inversion.ERT-Gauss: sigma_0_face = Summary: inf=0.00056949 sup=0.00056949 data points=2742
    Inversion.ERT-Gauss: sigma_0_stations = 8.968390466827157e-05 - 0.0012621902971677594 
    Inversion.ERT-Gauss: 8 additive DC potentials calculated.
    Inversion.ERT-Gauss: misfit ERT, reg, div, curl; total  =  2.585311e-02, 9.677289e-02, 2.003802e-01, 8.064624e-02; 2.147528e-01
    Inversion.ERT-Gauss: ratios ERT, reg  [%]       =  12.0385, 87.9615
    Inversion.LBFGS.linesearch.zoom: Iteration 4, alpha=0.0242991, phi(alpha)=0.214753 (interpolation order = 1)
    Inversion.ERT-Gauss: 8 adjoint potentials calculated.
    Inversion.LBFGS.linesearch.zoom: phi'(alpha)=0.000111966
    Inversion.LBFGS.linesearch.zoom: Zoom completed after 4 steps: alpha=0.0242991, phi(alpha)=0.214753, phi'(alpha)=0.000111966.
    Inversion.LBFGS.linesearch.zoom: Strong Wolfe condition is fulfilled. We are done.
    Inversion.LBFGS.linesearch: Line search completed after 1 steps (alpha=0.0242991, phi(alpha)=0.214753).
    Inversion.LBFGS: Search direction scaling found as alpha=0.0242991
    Inversion.LBFGS: Solution checked: |m-m_old|=0.000949099, |m|*m_tol=0.000368669
    Inversion.LBFGS: F(m) = 0.214753
    Inversion.LBFGS: Gradient has converged: |F-Fold|=9.80848e-07 < g_tol*max(|F|,|Fold|)=2.14754e-06
    Inversion.LBFGS: Success after 21 iterations!

### Preparation

Prespare mesh:

    ./prepare.py

This uses `brick.geo`. If you want to change the line set-up you need to change both files
`prepare.py` and `brick.geo`.

Then create synthetic data:

    runSynthetic.py -d config

Then inversion can run:

    runERTinversion -d config