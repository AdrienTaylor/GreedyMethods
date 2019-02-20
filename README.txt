Authors: Y. Drori 		Google Inc.,
         A. Taylor		INRIA/Ecole Normale Supérieure (Paris, PSL University).

Date:   Feb, 2019

Version: Feb 20, 2019

In this repository, you will find resources for the work entitled 
         "Efficient First-order Methods for Convex Minimization: a Constructive Approach";
it contents:

(1) [Folder SSEP_files] Files to generate the step sizes for the SSEP-based gradient method for smooth strongly
convex minimization. Those steps could be generated via PESTO, but the resulting code would be more difficult
to read.

(2) [Folder PESTO_files] Matlab files for generating the worst-case comparisons between the different methods 
(subgradient, fast gradient, triple momentum, optimized gradient, and SSEP-based methods) for miminimizing 
the different classes of functions treated in the work, using the Performance Estimation Toolbox
(available from https://github.com/AdrienTaylor/Performance-Estimation-Toolbox).

(3) [Folder Data] This folder contains the step-sizes generated via the SSEP. Those files can be used for
computing the worst-case guarantees for the resulting fixed-step procedures (e.g., via the PESTO toolbox
and the content of the PESTO_files folder.

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

