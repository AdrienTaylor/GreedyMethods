Authors: Y. Drori 		Google Inc.,
         A. Taylor		INRIA/Ecole Normale Supérieure (Paris, PSL University).

Date:   March, 2018

Version: March 13, 2018

In this repository, you will find resources for the work entitled 
         "Efficient First-order Methods for Convex Minimization: a Constructive Approach";
it contents:

(1) [Folder PESTO_files] Matlab files for generating the worst-case comparisons between the different methods 
(fast gradient, triple momentum and optimized gradient) for miminimizing smooth (strongly) 
convex functions, using the Performance Estimation Toolbox
(available from https://github.com/AdrienTaylor/Performance-Estimation-Toolbox).

(2) [Folder SSEP_files] Files to generate the step sizes for the SSEP-based gradient method for smooth strongly
convex minimization.

----- Setup

In order to use the code, the user should have YALMIP installed, along
with some semidefinite solver (e.g. Mosek, Sedumi, SDPT3, ...).

Once YALMIP and the SDP solver installed, the user can put the code in the
directory of his choice and execute the demo file step by step.

Note that figures from the work were generated using Mosek (other solvers should provide
similar results; usually not with the same accuracy).

----- Links

Link to YALMIP: http://users.isy.liu.se/johanl/yalmip/pmwiki.php?n=Main.Download

Link to Mosek: https://mosek.com/
Link to SeDuMi: http://sedumi.ie.lehigh.edu/
Link to SDPT3: http://www.math.cmu.edu/~reha/sdpt3.html

