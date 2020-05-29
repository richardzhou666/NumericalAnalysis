# NumericalAnalysis
This project implements a number of famous numerical methods in python.
* ivp.py
Numerical methods for solving first-order Initial Value Problems, including linear one step methods and Linear Multi-step Methods (LMM).
	* Runge Kutta Method
	* Euler Method
	* Trapzoidal Method
	* Leap Frog Method (Midpoint)
	* Adams Bashforth
	
* ivpSystem.py
Numerical methods for solving second-order Initial Value Problems in matrix form using Forward Euler algorithm as following:
$$
\begin{cases}
\vec{x}(t) = A \vec{x}(t), &t\in [0,3] \\
\vec{x}(0) = \vec{x_0}
\end{cases}$$
where$$
A = \begin{pmatrix}
-33.4 & 66.6 \\
33.3 & -66.7
\end{pmatrix}, \vec{x_0}=\begin{pmatrix}3 \\ 0 \end{pmatrix}$$
with the exact solution
$$\vec{x}(t)=\begin{cases} e^{-100t}+2e^{-t/10} \\ -e^{-100t}+e^{-t/10}\end{cases}$$

* heat.py
Numerical methods for solving classic one-dimensional heat equation with Dirichlet Boundary Conditions:
$$\begin{cases}
u_t(t,x) - u_{xx}(t,x) = 0 \\
u(0,x) = e^{-100x^2} \\
u(-1,t) = u(1,t) = 0
\end{cases}$$
	* Forward Euler
	* Backward Euler
	* Crank Nicolson
