# 18.303: Linear Partial Differential Equations: Analysis and Numerics

Spring 2021, Dr. Vili Heinonen, Dept. of Mathematics.

Overview
--------

This is the home page for the 18.303 course at MIT in Spring 2021, where the syllabus, lecture materials, problem sets, and other miscellanea are posted.

> **Course description**
>
> _Provides students with the basic analytical and computational tools of linear partial differential equations (PDEs) for practical applications in science engineering, including heat/diffusion, wave, and Poisson equations. Analytics emphasize the viewpoint of linear algebra and the analogy with finite matrix problems. Studies operator adjoints and eigenproblems, series solutions, Green's functions, and separation of variables. Numerics focus on finite-difference and finite-element techniques to reduce PDEs to matrix problems, including stability and convergence analysis and implicit/explicit timestepping. [Julia](http://julialang.org/) (a Matlab-like environment) is introduced and used in homework for simple examples._
>
> Prerequisite: linear algebra ([18.06](http://web.mit.edu/18.06), 18.700, or equivalent).

Syllabus
--------

**Lectures**: TR 9:30-11am [https://mit.zoom.us/j/92763506665](https://mit.zoom.us/j/92763506665). 
**Office Hours:** Wednesday 9-10am [https://mit.zoom.us/j/92763506665](https://mit.zoom.us/j/92763506665).

**Grading**: 50% homework, 15% mid-term, 35% final project
(due the last day of class). Problem sets are due in class on the due date.
The lowest problem set score will be dropped at the end of the term. Missed
midterms require a letter from Student Support Services or Student Disabilities
Services to justify accommodations. Legitimate excuses include sports,
professional obligations, or illness. In the event of a justified absence, an
alternative make-up project will be assigned.

**Collaboration policy**: Make an effort to solve the problem on your own before
discussing with any classmates. When collaborating, write up the solution on
your own and acknowledge your collaborators.

**Books**: Introduction to Partial Differential Equations by Olver.

**Final project**: There is a final project instead of a final exam. In your project,
you should consider a PDE or possibly a numerical
method not treated in class, and write a 5â€“10 page academic-style paper that
includes:

*Review*: why is this PDE/method important, what is its history, and what are
the important publications and references? (A comprehensive bibliography is
expected: not just the sources you happened to consult, but a complete set of
sources you would recommend that a reader consult to learn a fuller picture.)
*Analysis*: what are the important general analytical properties? e.g.
conservation laws, algebraic structure, nature of solutions (oscillatory,
decaying, etcetera). Analytical solution of a simple problem.
*Numerics*: what numerical method do you use, and what are its convergence
properties (and stability, for timestepping)? Implement the method (e.g. in
Julia) and demonstrate results for some test problems. Validate your solution
(show that it converges in some known case).

You must submit a one-page proposal of your intended final-project topic,
summarizing what you intend to do. Some suggestions of
possible projects will be given before then.

Tentative Schedule
--------------------

- Why PDEs are interesting.
- The Fourier series and eigenfunction expansions for the Poisson equation
- Optional: **Julia Tutorial**
- Spectral methods for numerically solving PDEs
- Finite difference discretizations
- Properties of Hermitian operators
- (Semilinear) Heat : Equation
- Basic time stepping methods
- Method of Lines (MOL) Solutions
- Lax equivalence, stability, Von Neumann Analysis
- Higher dimensional PDEs
- Generalized boundary conditions
- Separation of Variables
- Wave Equation
- Traveling waves and D'Alembert's solution
- Numerical Dispersion
- Sturm-Liouville Operators
- Distributions
- Green's Functions
- Weak form and Galerkin expansions
- Finite Element Methods

Lecture Summary
-------------------

## Vector spaces and linear operators
[Lecture 1](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture01/lecture01.pdf) | [Sine series (Julia)](https://github.com/mitmath/18303/blob/master/julia/Sine-series.ipynb)

During the first week we covered basic properties of vector spaces and linear operators including norms and inner products. We defined the notion of linear operators and introduced adjoints of linear operators. We went through some examples of smooth functions on the interval \[0,1\]. We used these tools and notions to solve the Poisson equation with Dirichlet boundary conditions on this interval. We also talked about bases for vector spaces and introduced the notion of an orthogonal and orthonormal bases. We introduced the basic properties of Fourier transform on a finite interval.

## Finite differences
[Lecture 2](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture02/lecture02.pdf) | [Finite differences (Julia)](https://github.com/mitmath/18303/blob/master/julia/finite_differences.ipynb)

We covered the basic idea of discretizing functions and writing down finite difference approximations of differential operators. We introduced backward, forward, and center difference methods and used these to write a simple discretization for the Laplacian. We talked about matrix representations for the difference operators and the importance of boundary conditions. We briefly discussed how the matrices are invertible if a linear system has a unique solution. As an example we talked about Poisson equation. We also covered deriving finite difference operators using polynomial fitting. 

## Heat and wave equations
[Lecture 3](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture03/lecture03.pdf) | [Pset 1](https://github.com/mitmath/18303/blob/master/problem_sets/pset1.ipynb)

We showed that Laplacian operator is self-adjoint with Dirichlet boundaries. We introduced the notion of positive and negative (semi)definite operators. We talked about the _superposition principle_ and used it to solve the heat equation and the wave equation with Dirichlet boundaries. Important theme during this lecture was the ability to separate partial differential equations in sufficiently symmetric domains. 

## Boundary conditions
[Lecture 4](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture04/lecture04.pdf)

We discussed some general properties of different boundary conditions for partial differential equations. We showed that the general solution is the solution to the homogeneous problem with the desired boundary condition + the solution to the inhomogeneous problem with zero boundaries. We revisited boundary conditions within the framework of finite difference approximation. We saw how the uniqueness of the solution to a linear PDE given by the boundary condition translates in manifested in the finite difference approximation. 

## Problems in higher dimensions
[Lecture 5](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture05/lecture05.pdf)

We talked about how some of the ideas we used earlier extend to problems in higher dimensions. The main method here was separation of variables, which is a powerful technique to solve PDEs when the problem is nicely symmetric. This is true especially for time-dependent problems -- time is usually independent of the spatial dimensions so solving for the time evolution once you have solved the spatial part is often not too hard. We also discussed calculating Fourier coefficients in greater detail and defined the finite Fourier transform as a linear map from one vector space to another. 

## On stability of numerical schemes
[Lecture 6](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture06/lecture06.pdf)

In this lecture we learned some techniques for analyzing the stability of time integration schemes. We considered a wave equation with damping controlled by the coefficient µ. The time evolution of the system can change drastically depending on µ. We briefly talked about what we mean by saying that a time integrator is stable. We saw how we can use the eigenvalues of the Laplacian to analyze the stability of the heat equation with different time discretization schemes. 

## Inhomogeneous PDEs revisited
[Lecture 7](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture07/lecture07.pdf)

We talked about general linear PDEs with various boundary conditions. One of the main points was that the solution of an inhomogeneous PDE can be seen as the solution to the homogeneous problem with the given boundary conditions + the solution to the inhomogeneous problem with zero boundaries. We solved the heat equation in n+1 dimensions. If the eigenproblem can be solved for the spatial part, it is always easy to solve the time evolution of the coefficients of the field in the eigenbasis. We showed how to do this.

We also solved the n+1 dimensional wave equation in the eigenbasis of the Laplacian. We introduced the _Fourier transform_ to solve for the time evolution of the coefficients.

## Distributions and Fourier transform
[Lecture 8](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture08/lecture08.pdf)

We introduced the notion of distributions: a generalization of functions that make sense under integration. Specifically, we discussed Dirac delta distribution that can be used to extract a value of a function at a given point by using an inner product between the Delta distribution and the function. We covered some of the basic features of the delta distribution. 

We formally defined the Fourier transform and discussed its properties. We showed how one can see the Fourier transform as a representation of a function in the Fourier basis i.e. as superposition of plane waves. We derived some identities for the Fourier transform using the properties of the delta distribution. We also covered the connection between Fourier coefficients and the Fourier transform.

## Wave equation and resonances
[Lecture 9](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture09/lecture09.pdf)

We introduced Laplace transform and discussed some of its properties (i.e. the connection with the Fourier transform). We talked about the good and the bad of the Laplace transform and pointed out some problems with calculating the inverse Laplace transform. In the end we came up with a reasonable formula for the inverse Laplace transform. 

We considered the time evolution of the driven harmonic oscillator (time evolution of the eigenbasis coefficients of the wave equation with an oscillating input). We solved this equation with an ansatz explaining the phenomenon of _resonance_. 

 ## Wave equation and resonances vol. 2
[Lecture 10](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture10/lecture10.pdf)

We did a recap of complex analysis: the basic properties of complex numbers, what it means to be complex differentiable and how to do contour integrals on the complex plane. These contour integrals are useful for dealing with (inverse) Fourier and Laplace transforms. We used the Laplace transform to solve the driven harmonic oscillator ODE from last time and covered also the case of resonant driving frequencies. 

 ## Finite differences in higher dimensions
[Lecture 11](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture11/lecture11.pdf)

We discussed the Poisson equation with general Dirichlet boundary conditions. After briefly discussing the 1d case and how we deal with the boundary conditions we formulated the problem in 2d. We showed how to build _block matrices_ that let us write down the Poisson equation as a simple linear system in the discretized setting. We also talked about boundary conditions in the 2d case. 

 ## Spectral methods and the weak formulation
[Lecture 12](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture12/lecture12.pdf)

We talked about the weak (variational) formulation of a boundary value problem. The idea is that a certain inner product has to hold for _all_ test functions that go smoothly to zero at the boundary. For computational applications we can use a basis for the function to be solved and the test functions that gives a sparse linear operator for a truncated set of functions. Usually these functions are such that it is easy to calculate their derivatives. 

 ## Some important basis functions
[Lecture 13](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture13/lecture13.pdf)

We solved the Laplace equation in spherical coordinates introducing the spherical harmonics and Associated Legendre polynomials. We discussed the properties of these basis functions in the framework of weak solutions to a partial differential equation. 