# MaxVol_Dual
 
This code solves the maxvol-dual problem

       $$ 
       \begin{align*} \max_{Z,\Theta,\Delta} & det(Z)^2 - \lambda ||\Delta||^2 \\
                 & \text{such that } Z = [\Theta; e^\top] \text{ and } Y^\top \Theta <= 1 + \Delta 
       \end{align*}
       $$

                 
### ****** Input ******
- X      :  the input matrix
- r      :  the rank of the sought approximation
- lambda :  the regularization parameter

  
### ****** Output ******
- v1          :    estimated center vector
- West        :    estimated W
- best_theta  :    final Theta at the convergence with maximum dual volume
- iter        :    number of iterations till convergence
- Y           :    projected points in r-1 dimensions
- C           :    projection matrix
