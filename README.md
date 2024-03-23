# MaxVol_Dual
 
This code solves the maxvol-dual problem
       $max_{Z,Theta,Delta} det(Z)^2 - lambda ||Delta||^2$
                 $such that Z = [Theta; e'] and Y' Theta <= 1 + Delta$
****** Input ******
- X      :  the input matrix
- r      :  the rank of the sought approximation
- lambda :  the regularization parameter
---Options---
- .maxiter    : the maximum number of iterations performed by the algorithm
               -default = 20.
- .epsilon    : the tolerance level for convergence
               -default = 1e-2.
- .num_workers: number of parallelized solutions used
****** Output ******

- v1          :    estimated center vector
- West        :    estimated W
- best_theta  :    final Theta at the convergence with maximum dual volume
- iter        :    number of iterations till convergence
- Y           :    projected points in r-1 dimensions
- C           :    projection matrix
