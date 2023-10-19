# System Realization by Information Matrix

## SRIM
$$
\begin{aligned}
\mathbf{x}(k+1) &= A\mathbf{x}(k) + B\mathbf{u}(k) \\
\mathbf{y}(k) &= C\mathbf{x}(k) + D\mathbf{u}(k)
\end{aligned}
$$

$$
\begin{aligned}
\mathbf{y}_{p}(k) &= \mathcal{O}_{p}\mathbf{x}(k) + \mathcal{T}_{p}\mathbf{u}(k) \\
\begin{bmatrix}
\mathbf{y}(k) \\
\mathbf{y}(k+1) \\
\vdots \\
\mathbf{y}(k+p-1)
\end{bmatrix}
&=
\begin{bmatrix}
C \\
CA \\ 
CA^{2} \\ 
\vdots \\
CA^{p-1}
\end{bmatrix}
\mathbf{x}(k)
+
\begin{bmatrix}
D \\ 
CB & D \\
CAB & CB & D \\
\vdots & \vdots & \vdots & \ddots \\
CA^{p-2}B & CA^{p-3} & CA^{p-4} & \cdots & D
\end{bmatrix}
\begin{bmatrix}
\mathbf{u}(k) \\
\mathbf{u}(k+1) \\
\vdots \\
\mathbf{u}(k+p-1)
\end{bmatrix}
\end{aligned}
$$

$$
\mathbf{y}_{p}(k) = 
\begin{bmatrix}
\mathbf{y}(k) \\
\mathbf{y}(k+1) \\
\vdots \\
\mathbf{y}(k+p-1)
\end{bmatrix}
$$

$$
\mathbf{u}_{p}(k) = 
\begin{bmatrix}
\mathbf{u}(k) \\
\mathbf{u}(k+1) \\
\vdots \\
\mathbf{u}(k+p-1)
\end{bmatrix}
$$

Stack them horizontally as columns in a matrix, to get the matrix equation:

$$
\begin{aligned}
Y_{p}(k) &= \mathcal{O}_{p}X(k) + \mathcal{T}_{p}U(k) \\
\end{aligned}
$$

$$
\begin{aligned}
Y_{p}(k) &= \begin{bmatrix} \mathbf{y}_{p}(k) & \mathbf{y}_{p}(k+1) & \cdots & \mathbf{y}_{p}(k+N-1) \end{bmatrix} \\
&= \begin{bmatrix}
\mathbf{y}(k)     & \mathbf{y}(k+1) & \cdots & \mathbf{y}(k+N-1)\\
\mathbf{y}(k+1)   & \mathbf{y}(k+2) & \cdots & \mathbf{y}(k+N)  \\
\vdots            & \vdots          & \ddots & \vdots \\
\mathbf{y}(k+p-1) & \mathbf{y}(k+p) & \cdots & \mathbf{y}(k+N+p-2)
\end{bmatrix}
\end{aligned}
$$

$$
X(k) = \begin{bmatrix} \mathbf{x}(k) & \mathbf{x}(k+1) & \cdots & \mathbf{x}(k+N-1) \end{bmatrix} 
$$

$$
\begin{aligned}
U_{p}(k) &= \begin{bmatrix} \mathbf{u}_{p}(k) & \mathbf{u}_{p}(k+1) & \cdots & \mathbf{u}_{p}(k+N-1) \end{bmatrix} \\
&= \begin{bmatrix}
\mathbf{u}(k)     & \mathbf{u}(k+1) & \cdots & \mathbf{u}(k+N-1)\\
\mathbf{u}(k+1)   & \mathbf{u}(k+2) & \cdots & \mathbf{u}(k+N)  \\
\vdots            & \vdots          & \ddots & \vdots \\
\mathbf{u}(k+p-1) & \mathbf{u}(k+p) & \cdots & \mathbf{u}(k+N+p-2)
\end{bmatrix}
\end{aligned}
$$

Then, multiply the entire thing by $\frac{1}{N}Y_{p}^{*}(k)$ on the right:

$$
\begin{aligned}
Y_{p}(k) &= \mathcal{O}_{p}X(k) + \mathcal{T}_{p}U(k) \\
\end{aligned}
$$

$$
R_{yy} = \mathcal{} 
$$

$$
R_{yy} = \frac{1}{N}Y_{p}(k)Y_{p}^{*}(k) 
$$

## Power Spectral Density
$$
R_{yy} = \frac{1}{N}Y_{p}(k)Y_{p}^{*}(k)
$$

$$
\begin{aligned}
R_{yy}(\tau) &= \mathbb{E}\left\{\mathbf{y}(t+\tau)\mathbf{y}^{*}(t)\right\} \\
&= \int_{-\infty}^{\infty}{\mathbf{y}(t+\tau)\mathbf{y}^{*}(t)f_{y}(t)dt}
\end{aligned}
$$
