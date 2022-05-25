### Modern Spectral Analysis

#### Ch. 1 Introduction
* The Core question: given time samples such as $Ae^{j(2\pi{}fn+\phi)}$, we need to estimate $A,f,\phi$
* Convetional vs. Modern methods
    * Convetional methods: Using DFT 
    * Modern Approach: considering the structures of the data model

#### Ch. 2 Signal Models and Adaptive Signal Processing
* Autocorrelation $r_x(k)=\mathbb{E}[x(n)x^*(n-k)]$
* Stochastic Models:
    * $v(n),u(n)$: input and output of System
    * $r_u(k)=h(k)*r_v(k)*h^*(-k)$
    * $S_u(z)=H(z)H^*(1/z^*)S_v(z)$
    * $S_u(e^{j2\pi{}f})=|H(e^{j2\pi{}f})|^2S_v(e^{j2\pi{}f})$

#### Ch. 3 Nonparametric Methods
* DFT: $X(k)=\sum_{n=0}^{N-1}x(n)e^{-j2\pi{}nk/N}$
* Power Spectral Density: $S_x(e^{j2\pi{}f})=lim_{N\rightarrow{}\infty}\mathbb{E}[\frac{1}{N}|\sum_{n=0}^{N-1}x(n)e^{-j2\pi{}nf}|^2]$
* Periodogram: $\hat{S}_{x,per}(e^{j2\pi{}f})=\frac{1}{N}|\sum_{n=0}^{N-1}\tilde{x}(n)e^{-j2\pi{}nf}|^2$
    * For $f_k=\frac{k}{N}$, we have $\hat{S}_{x,per}(e^{j2\pi{}f_k})=\frac{1}{N}|\tilde{X}(k)|^2$
    * $\hat{S}_{x,per}(e^{j2\pi{}f_k})=(N|\tilde{x}(n)*h_k(n)|^2|_{n=N-1}$ where $h_k(n)=\frac{1}{N}e^{j2\pi{}f_kn}$ if $n=0\sim{}(N-1)$
    * $H_k(e^{j2\pi{}f})=\frac{1}{N}diric(f-f_k)e^{-j\pi{}(f-f_k)(N-1)}$ where $diric(x)=\frac{sin(\pi{}xN)}{sin(\pi{}x)}$
* Bias-Variance Decomposition: MSE = Bias$^2$ + Var
    * Bias of $\hat{S}_{x,per}(e^{j2\pi{}f})$ is $\mu-S_{x}(e^{j2\pi{}f}$) where $\mu=\mathbb{E}[\hat{S}_{x,per}(e^{j2\pi{}f})]$
    * Var of $\hat{S}_{x,per}(e^{j2\pi{}f})$ is $\mathbb{E}[|\hat{S}_{x,per}(e^{j2\pi{}f})-\mu|^2]$
* Periodogram Mean: $\mathbb{E}[\hat{S}_{x,per}(e^{j2\pi{}f})]=W_B(e^{j2\pi{}f})*S_{x}(e^{j2\pi{}f})$ where $W_B(e^{j2\pi{}f})=\frac{1}{N}diric(f)^2$
    * Smoothed PSD
    * It is biased
* Periodogram Covariance: $\sim{}S_x(e^{j2\pi{}f_1})S_x(e^{j2\pi{}f_2})(\frac{1}{N}diric(f_1-f_2))^2$
    * When $N$ is large, var $\sim{}S_x^2(e^{j2\pi{}f})$ in $f=0\sim{}0.5$ and var $\sim{}2S_x^2(e^{j2\pi{}f})$ in $f=0,0.5$
    * Since the variance is related to square of PSD, the variance is larger when the PSD of that freq. is large
* Blackman-Tukey (Periodogram Smoothing) method
    * $\hat{S}_{x,BT}(e^{j2\pi{}f})=\sum_{k=-N+1}^{N-1}\hat{r}_x(k)w_a(k)e^{-j2\pi{}fk}$ where $w_a(k)$ is a window function
    * $\hat{S}_{x,BT}(e^{j2\pi{}f})=\hat{S}_{x,per}(e^{j2\pi{}f})*W_a(e^{j2\pi{}f})$
    * Mean of $\hat{S}_{x,BT}(e^{j2\pi{}f})$ is $S_{x}(e^{j2\pi{}f})*W_B(e^{j2\pi{}f})*W_a(e^{j2\pi{}f})$
    * Var of $\hat{S}_{x,BT}(e^{j2\pi{}f})$ is $\sim{}\frac{E_w}{N}S^2_x(e^{j2\pi{}f})$ where $E_w=\sum_{k=-L+1}^{L-1}w_a^2(k)$ is the energy of the filter. The $\frac{E_w}{N}$ is called the variance reduction factor
* Barlett-Welch (Periodogram Averaging) method
    * $\hat{S}_{x,per,i}(e^{j2\pi{}f})=\frac{1}{L}|\sum_{n=0}^{L-1}\tilde{x}_i(n)e^{-j2\pi{}fn}|^2$ where $\tilde{x}_i(n)=\tilde{x}(iD+n)w(n)$ and $D$ is the offset distance
    * $\hat{S}_{x,PA}(e^{j2\pi{}f})=\frac{1}{K}\sum_{i=0}^{K-1}\hat{S}_{x,per,i}(e^{j2\pi{}f})=\frac{1}{KL}\sum_{i=0}^{K-1}|\sum_{n=0}^{L-1}\tilde{x}_i(n)e^{-j2\pi{}fn}|^2$
    * Mean of $\hat{S}_{x,PA}(e^{j2\pi{}f})$ is $\frac{1}{L}\int{}_{-0.5}^{0.5}S_x(e^{j2\pi{}\nu})|W(e^{j2\pi{}(f-\nu)})|^2d\nu$
    * Variance of $\hat{S}_{x,PA}(e^{j2\pi{}f})$ is $\sim\frac{1}{K}S_x^2(e^{j2\pi{}f})$, assuming the data segments are independent


#### Ch. 4 Harmonic Models and MVDR Spectrum Estimation
* Harmonic model: $x(n) = \sum_{p=1}^P\alpha_pe^{j2\pi{}f_pn}+w(n)$
    * Assumption
        * Magnitude $|\alpha_p|$ is deterministic
        * Phase $\angle{\alpha_p}$ is random (uniform, independent)
        * Freq. $f_p$ is deterministic
        * Noise $w(n)$ is WSS
    * Mean-value function $\mu_x(n)=0$
    * Autocorrelation function $r_x(k)=\sum_{p=1}^P\sigma{}_p^2e^{j2\pi{}f_pk}+\sigma{}_w^2\delta{}(k)$
    * PSD $S_x(e^{j2\pi{}f})=\sum_{p=1}^P\sigma{}_p^2\delta{}(f-f_p)+\sigma{}_w^2$
* Vector-form harmonic model:
    * $\mathbf{x}(x) = [x(n)\ x(n+1)\ \dots{}\ x(n+M-1)]^T=\mathbf{V}\mathbf{s}(n)+\mathbf{w}(n)$ 
    * $\mathbf{V}=[\mathbf{v}(f_1)\ \mathbf{v}(f_2)\ \dots{}\ \mathbf{v}(f_P)]$ and $\mathbf{v}(f)=[1\ e^{j2\pi{}f}\ \dots{}\ e^{j2\pi{}f(M-1)}]^T$
    * $\mathbf{s}(n)=[\alpha_1e^{j2\pi{}f_1n}\ \dots{}\ \alpha_Pe^{j2\pi{}f_Pn}]^T$
    * $\mathbf{R}=\mathbb{E}[\mathbf{x}(n)\mathbf{x}^H(n)]=\mathbf{V}\mathbf{R}_s\mathbf{V}^H+\sigma_w^2\mathbf{I}$ where $\mathbf{R}_s=\mathbb{E}[\mathbf{s}(n)\mathbf{s}^H(n)]=diag(\sigma_1^2,\dots{},\sigma_P^2)$
* Minimum Variance Distortionless Response (MVDR) Filter Design
    * $\mathbf{c}_{MVDR}=\argmin_{\mathbf{c}}\mathbb{E}[|y(n)|^2]=\argmin_{\mathbf{c}}\mathbf{c}^H\mathbf{R}\mathbf{c}$ subject to $\mathbf{c}^H\mathbf{v}(f)=1$ where $y(n)=\mathbf{c}^H\mathbf{x}(n)$
    * $\mathbf{c}_{MVDR}=\frac{\mathbf{R}^{-1}\mathbf{v}(f)}{\mathbf{v}^H(f)\mathbf{R}^{-1}\mathbf{v}(f)}$
    * MVDR spectrum: $\hat{S}_{x,MVDR}(e^{j2\pi{}f})=M\mathbb{E}[|y_{MVDR}(n)|^2]=\frac{M}{\mathbf{v}^H(f)\mathbf{R}^{-1}\mathbf{v}(f)}$
    * $\mathbf{R}$ needs to be estimate first
        * sample correlation matrix: $\hat{\mathbf{R}}=\frac{1}{K}\sum_{k=1}^K\tilde{\mathbf{x}}(k)\tilde{\mathbf{x}}^H(k)$
    * MVDR on harmonic model ($P=1$)
        * $\mathbf{R}=\sigma_1^2\mathbf{v}(f_1)\mathbf{v}^H(f_1)+\sigma_w^2\mathbf{I}$
        * $\mathbf{R}^{-1}=\frac{1}{\sigma_w^2}(\mathbf{I}-\frac{\sigma_1^2}{\sigma_w^2+M\sigma_1^2}\mathbf{v}(f_1)\mathbf{v}^H(f_1))$

#### Ch. 5 Parametric Methods for Line Spectra
* Signal and noise subspace
    * For $\mathbf{R}=\mathbf{V}\mathbf{R}_s\mathbf{V}^H+\sigma_w^2\mathbf{I}$
    * If $\mathbf{u}$ is in null space of $\mathbf{V}\mathbf{R}_s\mathbf{V}^H$, then $\mathbf{u}$ is the eigenvector of $\mathbf{R}$ with eigen value $\sigma_w^2$
    * If $\mathbf{u}$ is in eigen space of $\mathbf{V}\mathbf{R}_s\mathbf{V}^H$ with eigen value $\mu$, then $\mathbf{u}$ is the eigenvector of $\mathbf{R}$ with eigen value $\sigma_w^2+\mu$
    * If $M>P$, $\mathbf{R}$ has $P$ positive eigenvalue
    * Eigen decomposition: $\mathbf{V}\mathbf{R}_s\mathbf{V}^H=\mathbf{U}_s\mathbf{D}\mathbf{U}_s^H$
    * Let $\mathbf{U}=[\mathbf{U}_s|\mathbf{U}_n]$, then $\mathbf{R}=\mathbf{U}\begin{bmatrix}\mathbf{D}+\sigma_w^2\mathbf{I}_P & \mathbf{0}\\ \mathbf{0} & \sigma_w^2\mathbf{I}_{M-P}\end{bmatrix}\mathbf{U}^H$
* Multiple Signal Classification (MUSIC)
    * For $\mathbf{u}_p$ in noise subspace $\mathbf{V}\mathbf{R}_s\mathbf{V}^H\mathbf{u}_p=\mathbf{U}_s\mathbf{D}\mathbf{U}_s^H\mathbf{u}_p=0$
    * If $rank(\mathbf{V})=P$, $\mathbf{R}_s\mathbf{V}^H\mathbf{u}_p=0$
    * If $\sigma_i^2>0$, $\mathbf{V}^H\mathbf{u}_p=0$, that is $\mathbf{u}_p^H\mathbf{v}(f_i)=0$ for $i=1\sim{}P$ and $p=P+1\sim{}M$
    * MUSIC pseudo-specturm $P_{x,MUSIC}(f)=\frac{1}{|\mathbf{U}^H_n\mathbf{v}(f)|^2}=\frac{1}{\mathbf{v}^H(f)\mathbf{U}_n\mathbf{U}^H_n\mathbf{v}(f)}$
    * Implementation 
        * Use estimated signal $\hat{\mathbf{R}},\hat{\mathbf{U}}$ for $\mathbf{R},\mathbf{U}$
        * Find peaks by dense "grid"
* Root MUSIC
    * Gridless
    * Let $z=e^{j2\pi{}f}$, we have $\mathbf{v}(f)=[1\ z\ \dots{}\ z^{M-1}]^T=\mathbf{a}(z)$ and $\mathbf{v}(f)^H=\mathbf{a}^T(z^{-1})$
    * We want to solve $g(z)=\mathbf{a}^T(z^{-1})\mathbf{U}_n\mathbf{U}_n^H\mathbf{a}(z)=0$
    * [Self-Reciprocal roots] If $z=re^{j\phi}$ is a root of $g(z)$, then $z=\frac{1}{r}e^{j\phi}$ is also a root of $g(z)$, since $g(\frac{1}{r}e^{j\phi})=g^H(re^{j\phi})=0$.
    * [Degree of polynomial is $2(M-1)$] Degree of $z$ of $g(z)$ is ranging from $-M+1$ to $M-1$
    * Implementation
        * Calculate estimated noise subspace $\mathbf{U}_n$
        * Construct polynomial $P(z)=z^{M-1}\mathbf{a}^T(z^{-1})\hat{\mathbf{U}}_n\hat{\mathbf{U}}_n^H\mathbf{a}(z)$
        * Find roots of $P(z)$
        * Find roots that are inside the unit circle that are closest to the circle
        * Frequency estimation by $\hat{f}=\frac{1}{2\pi}\angle{\hat{z}_p}$
* Estimation of Signal Parameters via Rotational Invariance Techniques (ESPRIT)
    * Let $\mathbf{x}_M(n)=[x(n)\ x(n+1)\ \dots{}\ x(n+M-1)]^T$
    * $\mathbf{x}_M(n)=[\mathbf{x}_{M-1}(n)|x(n+M-1)]^T=[x(n)|\mathbf{x}_{M-1}(n+1)]^T$
    * Let $\mathbf{J}_1=[\mathbf{I}|\mathbf{0}]$ (remove last) and $\mathbf{J}_2=[\mathbf{0}|\mathbf{I}]$ (remove front)
    * $\mathbf{J}_1\mathbf{v}_M(f)=\mathbf{v}_{M-1}(f)$ and $\mathbf{J}_2\mathbf{v}_M(f)=\mathbf{v}_{M-1}(f)e^{j2\pi{}f}$
    * $\mathbf{J}_2\mathbf{V}_M=\mathbf{J}_1\mathbf{V}_M\boldsymbol{\Phi}$ where $\boldsymbol{\Phi}=diag([e^{j2\pi{}f_1}\dots{}e^{j2\pi{}f_P}])$
    * There exists an invertible transformation matrix $\mathbf{T}$ such that $\mathbf{V}_M=\mathbf{U}_s\mathbf{T}$
    * $\mathbf{J}_2\mathbf{U}_s\mathbf{T}=\mathbf{J}_1\mathbf{U}_s\mathbf{T}\boldsymbol{\Phi}$
    * Let $\mathbf{U}_{s,1}=\mathbf{J}_1\mathbf{U}_s$ and $\mathbf{U}_{s,2}=\mathbf{J}_2\mathbf{U}_s$, then $\mathbf{U}_{s,2}=\mathbf{U}_{s,1}(\mathbf{T}\boldsymbol{\Phi}\mathbf{T}^{-1})=\mathbf{U}_{s,1}\boldsymbol{\Psi}$
    * diag of $\boldsymbol{\Phi}$ is the eigenvalue of $\boldsymbol{\Psi}$
    * Implementation
        * $\hat{\boldsymbol{\Psi}}$ by estimated signal subspace $\hat{\mathbf{U}}_{s}$ using LS or TLS
        * $f_p=\frac{\angle{\phi_p}}{2\pi}$
    * $\hat{\mathbf{U}}_{s,2}+\hat{\mathbf{E}}_2=\hat{\mathbf{U}}_{s,1}\hat{\boldsymbol{\Psi}}$, we have $\hat{\boldsymbol{\Psi}}_{LS}=(\hat{\mathbf{U}}_{s,1}^H\hat{\mathbf{U}}_{s,1})^{-1}\hat{\mathbf{U}}_{s,1}^H\hat{\mathbf{U}}_{s,2}$
    * $\hat{\mathbf{U}}_{s,2}+\hat{\mathbf{E}}_2=(\hat{\mathbf{U}}_{s,1}+\hat{\mathbf{E}}_1)\hat{\boldsymbol{\Psi}}$, the goal is to find $\min|\hat{\mathbf{E}}_1|_F^2+|\hat{\mathbf{E}}_1|_F^2$
    
#### Ch. 6 Parameter Estimation and Cramer-Rao Bounds
* Parameter estimation
    * Available information 
        * prior information (ex. probabiltiy distribution, data model, known parameters)
        * observartions
    * Goal
        * Based on the available information, estimate the unknown parameter $\theta$
        * Design an estimater $\hat\theta$ for $\theta$
* Scalar parameter estimation
    * Observation data vectcor $\mathbf{x}=[x(0)\ x(1)\ \dots{}\ x(N-1)]^T$
    * Parameter of interest: real and scalar $\theta$ (Frequentist: deterministic; Bayesian: random)
    * Prior information: pdf of $\mathbf{x}$
    * Estimater $\hat\theta$ (random)
    * $MSE(\hat\theta)=Var(\hat\theta)+|Bias(\hat\theta)|^2$
    * Cramer-Rao Bound (CRB or CRLB)
        * Assume $\hat\theta$ is an unbiased estimator of $\theta$. The pdf is assumed to satified the regularity condition: $\mathbb{E}_x[\frac{\partial{}logp(\theta,\mathbf{x})}{\partial{}\theta}]=0$
        * The variance satisfies: $Var(\hat\theta)\geq{}\mathcal{I}(\theta)^{-1}$ where $\mathcal{I}(\theta)=-\mathbb{E}_x[\frac{\partial{}^2logp(\theta,\mathbf{x})}{\partial{}^2\theta}]$ is called the Fisher information
    * Example $x(n)=Acos(2\pi{}f_1n+\phi)+w(n)$, the variable $A,f_1,\phi,\sigma^2$ are deterministic but only $f_1$ is unknown
        * $s(n,f_1) = Acos(2\pi{}f_1n+\phi)$
        * PDF: $p(\mathbf{x},\theta)=\frac{1}{(2\pi{}\sigma^2)^{N/2}}exp(-\frac{1}{2\sigma^2}\sum_{n=0}^{N-1}(x(n)-s(n,f_1))^2)$
        * $\mathcal{I}(\theta)=\frac{1}{\sigma^2}\sum_{n=0}^{N-1}(\frac{\partial^2s(n,f_1)}{\partial^2f_1})^2$
        * Therfore, $Var(\hat{f}_1)\geq{}\frac{\sigma^2}{\sum_{n=0}^{N-1}(\frac{\partial^2s(n,f_1)}{\partial^2f_1})^2}$
        * Analysis:
        <img src='./CRB.jpg'>
* Vector parameter 
    * Assume $\hat{\boldsymbol\theta}(\mathbf{x})$ is an unbiased estimator of $\boldsymbol\theta$. The pdf is assumed to satified the regularity condition: $\mathbb{E}_x[\frac{\partial{}logp(\boldsymbol\theta,\mathbf{x})}{\partial{}\boldsymbol\theta}]=\mathbf{0}$
    * The variance satisfies: $Cov(\hat{\boldsymbol\theta})\geq{}\mathcal{I}^{-1}(\theta)$ where $[\mathcal{I}(\boldsymbol\theta)]_{i,j}=-\mathbb{E}_x[\frac{\partial{}^2logp(\theta,\mathbf{x})}{\partial{}\theta_i\partial\theta_j}]$ is called the Fisher information matrix
    * $Var(\hat{\theta}_l)\geq{}[\mathcal{I}^{-1}(\boldsymbol\theta)]_{l,l}$
    * Example $x(n)=\sum_{p=1}^P\alpha_pe^{2\pi{}f_pn}+w(n)$