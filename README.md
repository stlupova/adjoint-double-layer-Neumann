# adjoint-double-layer-Neumann
Code for https://arxiv.org/abs/2310.00188

For questions/comments contact Svetlana Tlupova (svetlana.tlupova@farmingdale.edu) 

Reference:
J. T. Beale, M. Storm, and S. Tlupova. The adjoint double layer potential on smooth surfaces in R^3 and the Neumann problem. arXiv; Cornell University Library, 2023. https://arxiv.org/abs/2310.00188

This code was developed as part of work supported by the National Science Foundation grant DMS-2012371.

## Testing instructions

### Running the code:
1.	cd into one of the subdirectories
2.	make updates to parameters as described in tests below
3.	compile using  `make`
4.	run using `./main.out`


### Test 1: Direct evaluation of the modified ADL on unit sphere.	Exact solution based on spherical harmonic (equation (67) in referenced paper).

1.	Code subdirectory: `ADL_sphere_ellipsoid`
2.	In utilities.h, set the following:
    * grid size (h=1/N_gridlines): N_gridlines=16, 32, 64, 128
    * order of regularization in constant “regularization”:
        * regularization=0: no regularization
        * regularization=3: 3rd order regularization 
        * regularization=5 (or anything else): 5th order regularization
    * example = 1
    * direct_or_IE = 1
    * ellipse_a = ellipse_b = ellipse_c = 1.0
3.	In main.cpp, set the desired delta:
    * DEL = 1.5 * pow(h,4.0/5.0)
    * DEL = 3 * h
    * DEL = 2 * h


### Test 2: Integral equation solution on unit sphere. Exact solution based on spherical harmonic (equation (67) in referenced paper).

1.	Code subdirectory: `ADL_sphere_ellipsoid`
2.	Set same values as Test 1, with the exception of in utilities.h
    * direct_or_IE = 2
3.	In addition, in utilities.h, set the following:
    * beta = 0.7
    * IE_tol = 1.e-08
    * IE_iter_max = 100


### Test 3: Integral equation solution on ellipsoid. Exact solution is a harmonic function given in equation (71) in referenced paper.

1.	Code subdirectory: `ADL_sphere_ellipsoid`
2.	Set same values as Test 2, with the exception of in utilities.h
    * example = 2
    * ellipse_a = 1.0, ellipse_b = ellipse_c = 0.5


### Test 4: Integral equation solution on molecular surface. Exact solution is a harmonic function given in equation (71) in referenced paper.

1.	Code subdirectory: `ADL_molecule`
2.	The difference between the codes is only in utilities.cpp, in the way the surface is defined.
3.	Set same values as Test 3, with the exception of in utilities.h
    * example = 2
    * ellipse_a = ellipse_b = ellipse_c = 0.0 (these can be anything really; they are not taken into account)
