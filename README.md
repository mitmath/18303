# 18.303: Linear Partial Differential Equations: Analysis and Numerics

Fall 2019, Dr. [Chris Rackauckas](http://chrisrackauckas.com/), Dept. of Mathematics.

Overview
--------

This is the home page for the 18.303 course at MIT in Spring 2019, where the syllabus, lecture materials, problem sets, and other miscellanea are posted.

> **Course description**
>
> _Provides students with the basic analytical and computational tools of linear partial differential equations (PDEs) for practical applications in science engineering, including heat/diffusion, wave, and Poisson equations. Analytics emphasize the viewpoint of linear algebra and the analogy with finite matrix problems. Studies operator adjoints and eigenproblems, series solutions, Green's functions, and separation of variables. Numerics focus on finite-difference and finite-element techniques to reduce PDEs to matrix problems, including stability and convergence analysis and implicit/explicit timestepping. [Julia](http://julialang.org/) (a Matlab-like environment) is introduced and used in homework for simple examples._
>
> Prerequisite: linear algebra ([18.06](http://web.mit.edu/18.06), 18.700, or equivalent).

Syllabus
--------

**Lectures**: MWF 1–2pm (2-146). **Office Hours:** Wednesday 2–3pm (2-347).

**Grading**: 45% homework, 25% mid-term (March 20th), 30% final project
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
summarizing what you intend to do, by Wednesday, March. 13th. Some suggestions of
possible projects will be given before then.

Tentative Schedule
--------------------

- Why PDEs are interesting.
- The Fourier series and eigenfunction expansions for the Poisson equation
- Optional: **Julia Tutorial** (5-7PM 32-141)
- Spectral methods for numerically solving PDEs (Problem Set 1)
- Finite difference discretizations
- Properties of Hermitian operators (Problem Set 2)
- (Semilinear) Heat : Equation
- Basic time stepping methods
- Method of Lines (MOL) Solutions (Problem Set 3)
- Lax equivalence, stability, Von Neumann Analysis
- Higher dimensional PDEs
- Generalized boundary conditions (Problem Set 4)
- Advection Equation
- Upwinding operators, Lax-Wendorff
- (Brief!) 1D conservation laws and Berger's Equation (Problem Set 5)
- Separation of Variables
- Wave Equation
- Traveling waves and D'Alembert's solution
- Numerical Dispersion
- Sturm-Liouville Operators (Problem Set 6)
- Distributions
- Min-Max Theorem
- Green's Functions (Problem Set 7)
- Stochastic (Partial) Differential Equations and Diffusion-Advection Equations
- Weak form and Galerkin expansions
- Finite Element Methods

Lecture Summary
-------------------

## Lecture 1

In lecture one we reviewed linear algebra to build a perspective for infinite
dimensional linear algebra. Topics like "What is a vector space?" were revisited
and the abstraction away from real numbers to an algebraic structure was
emphasized. Questions as to the difference between a matrix and a linear
operator were settled. With this in mind, we introduced linear algebra on
functions. The space of smooth functions forms a vector space. What are
some linear operators for smooth functions? What is a basis for smooth functions?
How do you write the infinite matrix for the derivative operator in the polynomial
basis? We then used these ideas to solve the Poisson equation u_xx = f by
diagonalizing the derivative operator and getting a solution in terms of
sine functions.

## Lecture 2

[Lecture Notes](https://github.com/mitmath/18303/blob/master/lecture_notes/2_fourier_series.pdf) ||
[Problem Set 1](https://github.com/mitmath/18303/blob/master/problem_sets/ps1.pdf)

In this lecture we revisit Fourier series and make the ideas more concrete.
We prove some properties of the Fourier series and transform and then write
down the logic for how we solve the Poisson equation using operator notation.
Using this formalism, we continue onto the Heat equation u_t = u_xx + f and
show how to solve the Heat equation by using the eigenfunction basis of the
Laplacian. Some simple facts about the Heat equation are then revealed by
this solution and plotted.

## Lecture 3

[Lecture Notes](https://github.com/mitmath/18303/blob/master/lecture_notes/3_finite_differences.pdf)

In this lecture some computational mathematics tools were demonstrated. Lyx,
Jupyter notebooks, and Weave.jl were shown as options for writing mathematical
documents and incorporating numerical results. A recap of the Heat equation was
shown, and an emphasis on problem conversion was given. The same tools were then
used to solve the Wave Equation. Then moved from global (spectral) bases to
forming a discrete local basis. For this basis, we took evenly-spaced points
at which to represent the function. The forward and central difference
approximations were derived and order of convergence was discussed.

## Lecture 4

[Lecture Notes](https://github.com/mitmath/18303/blob/master/lecture_notes/4_hermitian.pdf)

In this lecture we continued our discussion of differencing operators, proving
the second order convergence of the second derivative central differencing
operator.  Alternative derivations from polynomial interpolation and
Fornberg's algorithm were discussed. From there, we looked at the properties
of the discretization matrices and identified similarities with the Fourier
case, leading to the idea of self-adjoint operators as a generalization of
symmetric to infinite bases.

## Lecture 5

[Lecture Notes](https://github.com/mitmath/18303/blob/master/lecture_notes/4_hermitian.pdf)

In this lecture we continued the discussion of self-adjoint operators and found
it to be a bit more nuanced than we thought. The issue, was boundary conditions.
Thus, we worked to clarify the standard boundary conditions (Dirchlet, Neumann,
Robin), and showed how this effects the properties of the linear operator.

## Lecture 6

[Problem Set 2](https://github.com/mitmath/18303/blob/master/problem_sets/ps2.pdf)

In this lecture we began to explore the convergence of numerical methods on
the Heat Equation. First we derived the forward time centered space and
backwards time centered space approximations by using our previous work on
finite difference operators. We showed that we could alternatively think about
the discretization of these PDEs as a two part problem: discretizing just space
gives an ODE that can then be analyzed. From here, the Lax Equivalence Theorem
was introduced, and the consistency of these methods was demonstrated. But
are they stable?

## Lecture 7

[Problem Set 3](https://github.com/mitmath/18303/blob/master/problem_sets/ps3.pdf)
