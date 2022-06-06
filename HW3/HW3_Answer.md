## MSA HW2

### Problem 1
**a.**
$||\mathbf{w}||_1=|2|+|3|+|-1|+|0|=6$

**b.**
$||\mathbf{w}||_2=\sqrt{2^2+3^2+(-1)^2+0^2}=\sqrt{14}$

**c.**
$||\mathbf{w}||_\infty=\max\{|2|,|3|,|-1|,|0|\}=3$

**d.**
$\text{supp}(\mathbf{w})=\{1,2,3\}$

**e.**
$\text{card}(\text{supp}(\mathbf{w}))=3$

### Problem 2
According to [FR2013], a non-negative function is called a norm if 
(a) $||\mathbf{x}||=0$ if and only if $\mathbf{x}=0$
(b) $||\lambda\mathbf{x}||=|\lambda|\times{}||\mathbf{x}||$ for all scaler $\lambda$ and $\mathbf{x}$
(c) $||\mathbf{x}+\mathbf{y}||\leq{}||\mathbf{x}||+||\mathbf{y}||$
However, (b) does not hold for all pairs of $\mathbf{x}$ and $\mathbf{y}$. For instance, let $\lambda=2$ and $\mathbf{y}=[3]$, then $||\lambda\mathbf{y}||=||[6]||=1\neq{}2=|2|\times{}||[3]||$

### Problem 3
**a.**
The Young's Inequality states that 
$\gamma\delta\leq{}\frac{\gamma^p}{p}+\frac{\delta^q}{q}$
Let $\gamma_i=\frac{|u_i^*|}{||\mathbf{u}||_p}$ and $\delta_i=\frac{|v_i|}{||\mathbf{v}||_q}$, we have 
$\frac{|u_i^*|}{||\mathbf{u}||_p}\frac{|v_i|}{||\mathbf{v}||_q}\leq{}\frac{1}{p}(\frac{|u_i^*|}{||\mathbf{u}||_p})^p+\frac{1}{q}(\frac{|v_i|}{||\mathbf{v}||_q})^q$
summing by $i$ for both sides, we get
$\frac{\sum_i|u_i^*v_i|}{||\mathbf{u}||_p||\mathbf{v}||_q}\leq{}\frac{1}{p}\sum_i(\frac{u_i^*}{||\mathbf{u}||_p})^p+\frac{1}{q}\sum_i(\frac{v_i}{||\mathbf{v}||_q})^q=\frac{1}{p}+\frac{1}{q}=1$
$\rightarrow{}\sum_i|u_i^*v_i|\leq{}||\mathbf{u}||_p||\mathbf{v}||_q\ -(1)$
Also, by the triangular inequality, $\sum_i|u_i^*v_i|\geq{}|\sum_iu_i^*v_i|=|\mathbf{u}^H\mathbf{v}|\ -(2)$
The inequality $|\mathbf{u}^H\mathbf{v}|\leq{}||\mathbf{u}||_p||\mathbf{v}||_q$ then holds due to $(1)$ and $(2)$

**b.**
<!-- The equality holds for the Young's Inequality iff $\gamma^p=\delta^q$, that is $(\frac{|u_i^*|}{||\mathbf{u}||_p})^p=(\frac{|v_i|}{||\mathbf{v}||_q})^q$
$\leftrightarrow{}|u_i|^p(||\mathbf{v}||_q)^q=|u_i^*|^p(||\mathbf{v}||_q)^q=|v_i|^q(||\mathbf{u}||_p)^p$
Let $\alpha=(||\mathbf{v}||_q)^q$ and $\beta=(||\mathbf{u}||_p)^p$
$\leftrightarrow{}\alpha|u_i|^p=\beta|v_i|^q$ -->
$|\mathbf{u}^H\mathbf{v}|=|\sum_iu_i^*v_i|=|\sum_i\frac{|v_i|^q}{v_i||\mathbf{v}||_q^{q-1}}v_i|=|\sum_i\frac{|v_i|^q}{||\mathbf{v}||_q^{q-1}}|=|\frac{\sum_i|v_i|^q}{||\mathbf{v}||_q^{q-1}}|=||\mathbf{v}||_q$
Also, 
$||\mathbf{u}||_p=\sqrt[p]{\sum_i|u_i|^p}=\sqrt[p]{\sum_i[\frac{|v_i|^q}{|v_i|||\mathbf{v}||_q^{q-1}}]^p}=\sqrt[p]{\sum_i\frac{|v_i|^{pq}}{|v_i|^p||\mathbf{v}||_q^{pq-p}}}=\sqrt[p]{\sum_i\frac{|v_i|^{p+q}}{|v_i|^p||\mathbf{v}||_q^{q}}}=\sqrt[p]{\sum_i\frac{|v_i|^{q}}{||\mathbf{v}||_q^{q}}}=\sqrt[p]{\frac{||\mathbf{v}||_q^{q}}{||\mathbf{v}||_q^{q}}}=1$
Therefore,
$|\mathbf{u}^H\mathbf{v}|=||\mathbf{u}||_p||\mathbf{v}||_q$ when the condition holds.

### Problem 4
**a.**
Assume that $\mathbf{A}=[\mathbf{a}_1\ \mathbf{a}_2\ \dots{}\mathbf{a}_6]$
Since $\mathbf{a}_i\neq{}c\mathbf{a}_j$ for all pairs of $i,j$ and scaling factor $c$, the spark is not 2.
Since $\mathbf{a}_3-2\mathbf{a}_5-\mathbf{a}_6=\mathbf{0}$, $\mathbf{a}_3,\mathbf{a}_5,\mathbf{a}_6$ is linearly dependent.
Therefore, the spark of the matrix is 3.

**b.**
The kruskal rank is therefore 2 due to **a.**.
The linearly dependent sets are:
$\mathcal{S}=\{\{\mathbf{a}_i,\mathbf{a}_j\}|i\neq{}j\}$

### Problem 5
For $||\mathbf{z}||_0=0$, it is obviouly impossible.
For $||\mathbf{z}||_0=1$, it is clear that no column vectors in $\mathbf{A}$ is a scaled from $\mathbf{y}$.
For $||\mathbf{z}||_0=2$, we check if all pairs of column vectors in $\mathbf{A}$ are linear independent with $\mathbf{y}$.
MATLAB code:
```
A = [-2 -4 4 9 8; ...
-7 9 4 -4 5; ...
5 -3 0 6 3];
z = [-22; 17; -15];
for i = 1 : 4
    for j = i+1 : 5
        if rank([A(:,i) A(:,j) z]) < 3
            fprintf("a_%d & a_%d\n",i,j)
        end
    end
end
```
The results find that only $\mathbf{a}_2,\mathbf{a}_4,\mathbf{y}$ are linearly dependent and this forms a unique solution.
Therefore, the solution set is:
$\mathcal{S}_\mathbf{z}=\{[0, 1, 0, -2, 0]^T\}$

### Problem 6
**a.**
<img src="./figure_jpg/HW3_6a.jpg" width=70%>
$Var(f): [9.38\times{}10^{-6}, 2.06\times{}10^{-6}, 6.87\times{}10^{-6}, 3.14\times{}10^{-6}, 2.35\times{}10^{-6}, 3.57\times{}10^{-5}]$
$Mean(f): [-2.50\times{}10^{-1}, -1.50\times{}10^{-1}, 9.96\times{}10^{-3}, 1.00\times{}10^{-1}, 2.00\times{}10^{-1}, 3.01\times{}10^{-1}]$
$MSE(f): [9.40\times{}10^{-6}, 2.06\times{}10^{-6}, 6.87\times{}10^{-6}, 3.14\times{}10^{-6}, 2.35\times{}10^{-6}, 3.65\times{}10^{-5}]$

**b.**
<img src="./figure_jpg/HW3_6b.jpg" width=70%>
$Var(f): [9.38\times{}10^{-6}, 2.06\times{}10^{-6}, 6.87\times{}10^{-6}, 3.14\times{}10^{-6}, 2.36\times{}10^{-6}, 3.56\times{}10^{-5}]$
$Mean(f): [-2.50\times{}10^{-1}, -1.50\times{}10^{-1}, 9.95\times{}10^{-3}, 1.00\times{}10^{-1}, 2.00\times{}10^{-1}, 3.01\times{}10^{-1}]$
$MSE(f): [9.40\times{}10^{-6}, 2.06\times{}10^{-6}, 6.87\times{}10^{-6}, 3.14\times{}10^{-6}, 2.36\times{}10^{-6}, 3.64\times{}10^{-5}]$