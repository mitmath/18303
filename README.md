# 18.303: Linear Partial Differential Equations: Analysis and Numerics

This is the main repository of course materials for 18.303 at MIT, taught by Dr. Andrew Horning, in Fall 2023. The syllabus is attached in the README below. Lecture notes, problem sets, exams, and supplementary materials will be posted in the repo directory.

> **Course description**
>
> _Provides students with the basic analytical and computational tools of linear partial differential equations (PDEs) for practical applications in science engineering, including heat/diffusion, wave, and Poisson equations. Analytics emphasize the viewpoint of linear algebra and the analogy with finite matrix problems. Studies operator adjoints and eigenproblems, series solutions, Green's functions, and separation of variables. Numerics focus on finite-difference and finite-element techniques to reduce PDEs to matrix problems, including stability and convergence analysis and implicit/explicit timestepping. Some programming required for homework and final project._
>
> Prerequisite: linear algebra ([18.06](http://web.mit.edu/18.06), 18.700, or equivalent).

## Syllabus

**Lectures**: MW 11:00 am - 12:30 pm in Room 2-142. 

**Office Hours:** 4pm on Tuesday and Wednesday in Room 2-238C.

**Grading**: 50% homework, 15% mid-term, 35% final project
(due the last day of class). Problem sets are due in class on the due date. Missed
midterms require a letter from Student Support Services or Student Disabilities
Services to justify accommodations. Justifiable reasons for absence include sports,
professional obligations, or illness. In the event of a justified absence, a make-up assignment will be provided.

**Collaboration policy**: Set aside time to work on each problem independently before
discussing it with classmates. Always write up the solution on
your own and acknowledge your collaborators.

**Books**: [Introduction to Partial Differential Equations](https://www-users.cse.umn.edu/~olver/pde.html) by Olver.

**Final project**: Instead of a final exam, study a PDE and/or numerical
solver not covered in class, and write a 5-10 page academic-style paper that
includes:

> *Review*: why is this PDE/method important, what is its history, and what are
the important publications and references? (A comprehensive bibliography is
expected: not just the sources you happened to consult, but a complete set of
sources you would recommend that a reader consult to learn a fuller picture
>
> *Analysis*: what are the important general analytical properties? e.g.
conservation laws, algebraic structure, nature of solutions (oscillatory,
decaying, etc.). Analytical solution of a simple problem.
>
> *Numerics*: what numerical method do you use, and what are its convergence
properties (and stability, for timestepping)? Implement the method (e.g. in
Julia) and demonstrate results for some test problems. Validate your solution
(show that it converges in some known case).

You must submit a one-page proposal of your intended final-project topic,
summarizing what you intend to do. Potential projects will be suggested as the course progresses.

## Assignments

- [Homework 1](https://github.com/mitmath/18303/blob/master/problem_sets/hw1.pdf) is due at 10pm on Friday, September 22 ([Solutions](https://github.com/mitmath/18303/blob/master/problem_sets/hw1_soln.pdf)).
- [Homework 2](https://github.com/mitmath/18303/blob/master/problem_sets/hw2.pdf) is due at 10pm on Friday, October 6 ([Solutions](https://github.com/mitmath/18303/blob/master/problem_sets/hw2_soln.pdf)).
- [Homework 3](https://github.com/mitmath/18303/blob/master/problem_sets/hw3.pdf) is due at 10pm on Friday, October 27.

## Tentative Schedule

- Why study *linear* PDEs?
- The "linear algebra" of taking derivatives
- Optional: **Julia Tutorial**
- Poisson's equation: diagonalizing a linear operator
- Numerical approximation: global or local?
- The heat equation: operator exponentials
- Numerical approximations of operator exponentials
- Convergence, stability, and conditioning
- Wave equations: analysis and numerics
- Higher dimensional PDEs and systems of PDEs
- Generalized boundary conditions and well-posedness
- Green's function and singular integral equations
- Weak form and Galerkin's method
- Finite Element Methods

## Lecture Summaries

### Lecture 1

- What are linear partial differential equations (PDEs) and why should we solve them?
- Example: the transport equation and its characteristic curves.

[Notes](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture_1.pdf) | Olver, Chapter 1

### Lecture 2

- First-order linear PDEs.
- Reduction to an ODE system: the method of characteristics.
- Computation on a grid: finite differences.

[Notes](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture_2.pdf) | Olver, Chapter 2.1-2.2

### Lecture 3

- Finite differences: approximation errors.
- Time-stepping and spatial discretizations.
- A linear algebraic roadmap for linear PDEs.

[Notes](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture_3.pdf) | Olver, Chapter 5.1

### Lecture 4

- Function spaces and bases.
- The nullspace of a differential operator.
- Boundary conditions and invertibility.
- Inner products, norms, and integrability.

[Notes](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture_4.pdf) | [Linear Algebra with Functions](https://github.com/mitmath/18303/blob/master/supp_material/linear_algebra_with_functions.pdf)


### Lecture 5

- Diagonalizing differential operators.
- Fourier's basis and its properties.
- Poisson's equation: separation of variables.

[Notes](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture_5.pdf) | Olver, Chapter 3.2

### Lecture 6

- Poisson's equation with homogeneous boundary conditions.
- Enforcing homogeneous boundary conditions through a basis.
- Enforcing inhomogeneous boundary conditions via nullspace.

[Notes](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture_6.pdf) | Olver, Chapter 4.3

### Lecture 7

- Diagonalization in separable domains.
- Laplace's equation in a disk.
- Poisson's formula for the disk.

[Notes](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture_7.pdf) | Olver, Chapter 4.3

### Lecture 8

- Maximum and minimum principles.
- Finite differences in 2 dimensions.

[Notes](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture_8.pdf) | Olver, Chapter 4.3 | [poissonFD.ipynb](https://github.com/mitmath/18303/blob/master/supp_material/poissonFD.ipynb)

### Lecture 9

- Kronecker products and the "vec" identity.
- Finite difference discretizations on separable domains.
- Fourier spectral discretizations on periodic domains.

[Notes](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture_9.pdf)

### Lecture 10

- The heat equation and its derivation.
- Boundary conditions and their physics.
- Fourier series solution and the operator exponential.

[Notes](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture_10.pdf) | Olver, Chapter 4.1

### Lecture 11

- Operator exponentials of self-adjoint operators.
- Smoothness, asymptotics, and equilibrium solutions.
- Homogeneous vs inhomogeneous boundary conditions.

[Notes](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture_11.pdf) | Olver, Chapter 4.1

### Lecture 12

- Heat flow in a rod with temperature fixed at endpoints.
- Heat flow in a disk with temperature prescribed on boundary.
- Dirichlet eigenfunctions of the Laplacian in a disk.

[Notes](https://github.com/mitmath/18303/blob/master/lecture_notes/lecture_12.pdf) | Olver, Chapter 4.1
