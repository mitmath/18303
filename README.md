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
**Office Hours:** Wednesday 10-11am [https://mit.zoom.us/j/92763506665](https://mit.zoom.us/j/92763506665).

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
method not treated in class, and write a 5–10 page academic-style paper that
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

Vector spaces and linear operators
-------------------
During the first week we covered basic properties of vector spaces and linear operators including norms and inner products. We defined the notion of linear operators and introduced adjoints of linear operators. We went through some examples of smooth functions on the interval \[0,1\]. We used these tools and notions to solve the Poisson equation with Dirichlet boundary conditions on this interval. We also talked about bases for vector spaces and introduced the notion of an orthogonal and orthonormal bases. We introduced the basic properties of Fourier transform on a finite interval.
