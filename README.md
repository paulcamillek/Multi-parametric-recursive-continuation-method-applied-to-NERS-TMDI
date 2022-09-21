# Multi-parametric-recursive-continuation-method-applied-to-NERS-TMDI
This set of codes can help continue the Frequency Response of a System into multiple dimensions. 
The main file is named main_function and all other functions all called directly through it. 
The code first generates Frequency response of the system (FRF-function). Limit Point Cycles (LPC), Neimark Sacker (NS), and other bifurcation points are tracked 
within branch follow. 
The newton function computes the solution to the nonlinear equations in nondim_temp (nondimensional temporary files). 
Once a bifurcation point is found, the user is asked if he wants to continue along the FRF or track the bifurcation point, for example a (LPC). 
There is also an option to compute the bifurcation diagram in F-xn space. L1 is for a level 1 continuation and has its related files. 
L2 is a level 2 continuation and has its related files. More levels can be computed following the format of L1 and L2 and referring to 
[https://doi.org/10.1016/j.ymssp.2019.03.011].
Please beware that some function will need to be significantly modified for ones application
For further questions about this code, please email the author @ paulc97@vt.edu
