# Maxium Volume Simplex Dual (MV-Dual) 
 
This MATLAB code allows one to solve the simplex-structured matrix factorization (SSMF) problem: Given an $m$-by-$n$ data matrix $X$ and a factorization rank $r$, it looks for a matrix $W$ with $r$ columns and a column-stochastic matrix $H$ with $r$ (entries in each column are nonnegative and sum to one). 

This code solves SSMF via a dual approach, by solving the dual problem

```math
\begin{aligned} \max_{Z,\Theta,\Delta} & \quad det(Z)^2 - \lambda ||\Delta||^2 \\
& \text{such that } Z = [\Theta; e^\top] \text{ and } Y^\top \Theta <= 1 + \Delta , 
\end{aligned}
```

where 
- $Y = U^T (X-v e^T)$ with $U$ containing the first $r$ singular vectors of $X-v e^T$,
- $v$ is a translation vector contained in the interior of the convex hull of the columns of $X$ (e.g., $v = Xe/n$),
- $e$ is the vector of all ones, and 
- $\lambda$ is a penalty parameter balancing volume and noise.

The columns of the variable $\Theta$ contains an approximation of the facets of the convex hull of the volumns of $W$, and $\Delta$ models the noise.  

See the paper "Dual Simplex Volume Maximization for Simplex-Structured Matrix Factorization", by M. Abdolali, G. Barbarino and N. Gillis, 2024. 
