# Maxium Volume Simplex Dual (MV-Dual) 
 
This MATLAB code allows one to solve the simplex-structured matrix factorization (SSMF) problem: Given a data matrix $X \in \mathbb{R}^{m \times n}$ and a factorization rank $r$, it looks for a matrix $W \in \mathbb{R}^{m \times r}$ with $r$ columns and a column-stochastic matrix $H \in \mathbb{R}^{r \times n}_+$ with $r$ (entries in each column are nonnegative and sum to one) such that $X \approx WH$. 

This code solves SSMF via a dual approach, by solving the dual problem

```math
\begin{aligned} \max_{Z,\Theta,\Delta} & \quad \text{det}(Z)^2 - \lambda ||\Delta||^2 \\
 \text{such that } &  \quad Z = [\Theta; e^\top] \text{ and } Y^\top \Theta \leq 1 + \Delta , 
\end{aligned}
```

where 
- $Y = U^T (X-v e^T)$ with $U$ containing the first $r$ singular vectors of $X-v e^T$,
- $v$ is a translation vector contained in the interior of the convex hull of the columns of $X$ (e.g., $v = Xe/n$),
- $e$ is the vector of all ones, and 
- $\lambda$ is a penalty parameter balancing volume and noise.

The matix $Y$ is a low-dimensional projection of $X$, after translation around the origin. The columns of the variable $\Theta$ contains an approximation of the facets of the convex hull of the columns of $W$; in other words, the convex hull of the columns of $\Theta$ represents the polar of the convex hull of the columns of $W$. The matrix $\Delta$ models the noise.  

See the paper "Dual Simplex Volume Maximization for Simplex-Structured Matrix Factorization", by M. Abdolali, G. Barbarino and N. Gillis, 2024. 
