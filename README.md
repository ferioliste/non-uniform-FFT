# NUFFT-II Implementation

This repository contains the implementation of the Non-uniform discrete Fast Fourier Transform of Type II (NUFFT-II) presented in [[Ruiz–Antol and Townsend, 2017]](https://arxiv.org/pdf/1701.04492) using Chebyshev polynomials and Taylor series expansions. 
We consider a sequence of complex numbers $c_1,...,c_{N-1}$, the sample points $x_0,...,x_{N-1}$ such that
```math
\left|x_j-\frac{j}{N}\right|\leq\frac{\gamma}{N}\qquad 0\leq j\leq N-1,
```
with $0\leq\gamma\leq 1/2$, and integer frequencies $\omega_j = j$ for $j=0,...,N-1$. We compute
```math
f_j = \sum_{k=0}^{N-1} c_k e^{-2\pi i x_k \omega_j} \qquad 0\leq j\leq N-1,
```
by contructing a low rank approximation of the transformation matrix making use of the FFT algorithm. Given a working precision $0<\varepsilon<1$, this allows to achieve a cost of $O(N\ \log N \ \log(1/\varepsilon) / \log\log(1/\varepsilon))$ in the Chebyshev case and of $O(N\ \log N\ \log(1/\varepsilon))$ when using a Taylor expansion. In general, this is better that $O(N^2)$ of the naive implementation.

## Files in the Repository

- `NUFFT_II_cheb`: Implementation using Chebyshev polynomials.
- `NUFFT_II_tayl`: Implementation using Taylor series expansion.
