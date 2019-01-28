# 18.303: Linear Partial Differential Equations: Analysis and Numerics

Fall 2013, Prof. [Steven G. Johnson](http://math.mit.edu/~stevenj/), Dept. of Mathematics.

Overview
--------

This is the home page for the 18.303 course at MIT in Fall 2013, where the syllabus, lecture materials, problem sets, and other miscellanea are posted.

You can also download the course announcement [flyer one](poster.pdf) or [two](poster2.pdf). Another version of 18.303 was previously offered as [18.303 by Prof. Freedman](http://math.mit.edu/classes/18.303/) (also [on OpenCourseWare](http://ocw.mit.edu/OcwWeb/Mathematics/18-303Fall-2006/CourseHome/index.htm)), but the content and focus of the previous 18.303 differs substantially from Prof. Johnson's syllabus.

> **Course description**
>
> _Provides students with the basic analytical and computational tools of linear partial differential equations (PDEs) for practical applications in science engineering, including heat/diffusion, wave, and Poisson equations. Analytics emphasize the viewpoint of linear algebra and the analogy with finite matrix problems. Studies operator adjoints and eigenproblems, series solutions, Green's functions, and separation of variables. Numerics focus on finite-difference and finite-element techniques to reduce PDEs to matrix problems, including stability and convergence analysis and implicit/explicit timestepping. [Julia](http://julialang.org/) (a Matlab-like environment) is introduced and used in homework for simple examples._
>
> Prerequisite: linear algebra ([18.06](http://web.mit.edu/18.06), 18.700, or equivalent).

Syllabus
--------

**Lectures**: MWF 1–2pm (4-159). **Office Hours:** Tues. 4–5pm (E17-416).

**Grading**: 45% homework, 25% mid-term (Nov. 1), 30% final project (due the last day of class). Problem sets are **due in class** on the due date, no excuses, but your lowest pset score will be dropped at the end of the term.

**Collaboration policy**: Talk to anyone you want to and read anything you want to, with two exceptions: First, make a solid effort to solve a problem on your own before discussing it with classmates or googling. Second, no matter whom you talk to or what you read, write up the solution on your own, without having their answer in front of you.

**Books**: [_Computational Science and Engineering_](http://math.mit.edu/cse/) by Strang is a recommended but not required text (emphasizing more the numerical part of the course). You may also find the free online book [Introduction to Partial Differential Equations](http://www.math.umn.edu/~olver/pdn.html) by Olver to be helpful.

**Final project:** There is a final project instead of a final exam. In your project, you should should consider a PDE (ideally in 2d; I would recommend against 3d for time reasons, and 1d problems may be too simple) or possibly a numerical method not treated in class, and write a 5–10 page academic-style paper that includes:

*   _Review:_ why is this PDE/method important, what is its history, and what are the important publications and references? (A comprehensive bibliography is expected: not just the sources you happened to consult, but a complete set of sources you would recommend that a reader consult to learn a fuller picture.)
*   _Analysis:_ what are the important general analytical properties? e.g. conservation laws, algebraic structure, nature of solutions (oscillatory, decaying, etcetera). Analytical solution of a simple problem.
*   _Numerics:_ what numerical method do you use, and what are its convergence properties (and stability, for timestepping)? Implement the method (e.g. in Matlab) and demonstrate results for some test problems. Validate your solution (show that it converges in some known case).

You must submit a **one-page proposal** of your intended final-project topic, summarizing what you intend to do, by **Monday, October 28**. Some suggestions of possible projects will be given before then.

* * *

Lecture Summaries and Handouts
------------------------------

### Lecture 1 (Sep 4)

**Handouts:** syllabus (this web page), [Lecture 1 notes](Lecture1.pdf), [pset 1](pset1-f13.pdf) (due Wed. Sep. 11)

General overview of what a PDE is and why they are important. Discussed examples of some typical and important PDEs (see handout, page 1). With non-constant coefficients (the most common case in real-world physics and engineering), even the simplest PDEs are rarely solvable by hand; even with constant coefficients, only a relative handful of cases are solvable, usually high-symmetry cases (spheres, cylinders, etc.) solvable. Therefore, although we will solve a few simple cases by hand in 18.303, the emphasis will instead be on two things: learning to _think_ about PDEs by recognizing how their _structure_ relates to concepts from finite-dimensional linear algebra (matrices), and learning to _approximate_ PDEs by actual matrices in order to solve them on computers.

Went through 2nd page of handout, comparing a number of concepts in finite-dimensional linear algebra (ala 18.06) with linear PDEs (18.303). The things in the "18.06" column of the handout should already be familiar to you (although you may need to review a bit if it's been a while since you took 18.06)—this is the kind of thing I care about from 18.06 for this course, not how good you are at Gaussian elimination or solving 2×2 eigenproblems by hand. The things in the "18.303" column are perhaps unfamiliar to you, and some of the relationships may not be clear at all: what is the dot product of two functions, or the transpose of a derivative, or the inverse of a derivative operator? Unraveling and elucidating these relationships will occupy a large part of this course.

Covered the concept of **nondimensionalization**: rescaling the units so that dimensionful constants and other odd numbers disappear, making as many things "1" as possible. Gave an example of a heat equation κ∇2T = ∂T/∂t in an L×L box in SI units, where we have a thermal conductivity κ in m2/s. By rescaling the spatial coordinates to x/L and y/L, and rescaling the time coordinate to κt/L2, we obtained a simplified equation of the form ∇2T = ∂T/∂t in a 1×1 box. Not only does this simplify the equations, but it can also improve our understanding: by rescaling with _characteristic times and distances_, we are left with distance and time units where 1 is the characteristic time and distance, and so in these units it is immediately obvious what we should consider "big" and "small". For example, in the rescaled time units, 0.01 is a small time in which probably not much happens, while 100 is a big time in which the solution has probably changed a lot. In the original SI units we would have had to explicitly compare to the characteristic time L2/κ.

### Lecture 2 (Sep 6)

**Handouts:** [Fourier sine series examples](http://web.mit.edu/18.06/www/Spring09/sines.pdf)

Started with a very simple vector space V of functions: functions u(x) on \[0,L\] with u(0)=u(L)=0 (Dirichlet boundary conditions), and with one of the simplest operators: the 1d Laplacian Â=d2/dx2. Explained how this describes some simple problems like a stretched string, 1d electrostatic problems, and heat flow between two reservoirs.

Inspired by 18.06, we begin by asking what the null space of Â is, and we quickly see that it is {0}. Thus, any solution to Âu=f must be unique. We then ask what the eigenfunctions are, and quickly see that they are sin(nπx/L) with eigenvalues -(nπ/L)2. If we can expand functions in this basis, then we can treat Â as a number, just like in 18.06, and solve lots of problems easily. Such an expansion is precisely a Fourier sine series (see handout).

In terms of sine series for f(x), solve Âu=f (Poisson's equation) and Âu=∂u/∂t with u(x,0)=f(x) (heat equation). In the latter case, we immediately see that the solutions are decaying, and that the high-frequency terms decay faster...eventually, no matter how complicated the initial condition, it will eventually be dominated by the smallest-n nonzero term in the series (usually n=1). Physically, diffusion processes like this smooth out oscillations, and nonuniformities eventually decay away. Sketched what the solution looks like in a typical case.

As a preview of things to come later, by a simple change to the time-dependence found a solution to the wave equation Âu=∂2u/∂t2 from the same sine series, which gives "wavelike" behavior. (This is an instance of what we will later call a "separation of variables" technique.)

**Further reading:** Section 4.1 of the Strang book (Fourier series and solutions to the heat equation).

### Julia Tutorial (Sep 9)

**Handout:** [Julia cheat-sheet](https://github.com/mitmath/julia-mit/blob/master/Julia-cheatsheet.pdf)

On Monday, 9 September, at 6pm in 32-123, I will give an (attendance-optional!) Julia tutorial, introducing the [Julia programming language and environment](http://julialang.org/) that we will use this term. I'll be posting some [tutorial notes online](https://github.com/mitmath/julia-mit/blob/master/README.md) as well (currently in draft form).

Feel free to bring your laptops; try to download Julia first, and I will try to help answer installation questions as well.

### Lecture 3 (Sep 9)

**Handouts:** [notes on difference approximations](difference-approx.pdf)

Now, we will go back to the happy land of finite-ness for a while, by learning to approximate a PDE by a matrix. This will not only give us a way to compute things we cannot solve by hand, but it will also give us a different perspective on certain properties of the solutions that may make certain abstract concepts of the PDE clearer. We begin with one of the simplest numerical methods: we replace the continuous space by a grid, the function by the values on a grid, and derivatives by differences on the grid. This is called a **finite-difference method**.

Went over the basic concepts and accuracy of approximating derivatives by differences; see handout.

Armed with center differences (see handout), went about approximating the 1d Laplacian operator d2/dx2 by a matrix, resulting in a famous tridiagonal matrix known as a _discrete Laplacian_. The properties of this matrix will mirror many properties of the underlying PDE, but in a more familiar context. We already see by inspection that it is real-symmetric, and hence we must have real eigenvalues, diagonalizability, and orthogonal eigenvectors—much as we observed for the d2/dx2 operator—and in the next lecture we will show that the eigenvalues are negative, i.e. that the matrix is negative-definite.

The negative eigenvalues mean that the discrete Laplacian is negative definite, and also suggest that it can be written in the form -DTD for some D. Reviewed the proof that this means the matrix is negative definite, which also relies on D being full column rank. Began to show that our A indeed has this form: we derived the discrete Laplacian by turning two derivatives into differences, one by one, and now by writing the first step as a matrix we get D, while writing the second step as a matrix shows that it is -DT. To get a negative _definite_ matrix (as opposed to just negative semidefinite), we additionally require that D be full column rank; showed that this is easy to see from DT since it is upper-triangular.

**Further reading:** [notes on finite-difference approximations from 18.330](http://homerreid.ath.cx/teaching/18.330/Notes/NumericalDifferentiation.pdf)

### Lecture 4 (Sep 11)

**Handouts:** [pset 2](pset2-f13.pdf) (due next Wed.), [pset1 solutions](pset1sol-f13.pdf) (see also IJulia [notebook](pset1sol-f13.ipynb); in [html](pset1sol-f13.html))

To do a similar analysis of the actual Laplacian, we first have to have a dot product, or **inner product**. Defined an abstract ⟨u,v⟩ notation (a map from _functions_ u and v to _scalars_ ⟨u,v⟩) for inner products, as well as three key properties. First, ⟨u,v⟩ = complex conjugate of ⟨v,u⟩. Second, |u|2\=⟨u,u⟩ must be nonnegative, and zero only if u=0. Third, it must be linear: ⟨u,αv+βw⟩=α⟨u,v⟩+β⟨u,w⟩. (Note: some textbooks, especially in functional analysis, put the conjugation on the second argument instead of the first.) For functions, the most common inner product (though _not the only choice_ and not always the _best_ choice, as we will see next time) is a simple integral ∫uv (conjugating u for complex functions); we will look at this more next time.

Reviewed inner products of functions. A vector space with an inner product (plus a technical criterion called "completeness" that is always satisfied in practice) is called a **Hilbert space**. Note that we include only functions with finite norm ⟨u,u⟩ in the Hilbert space (i.e. we consider only [square-integrable](Square-integrable
function) functions), which throws out a lot of divergent functions and means that everything has a convergent Fourier series. (Another omitted technicality: we have to ignore finite discrepancies at isolated points, or otherwise you can have ⟨u,u⟩=0 for u(x) nonzero; there is a rigorous way to do this, which we will come back to later.)

Defined the **adjoint** Â\* of a linear operator: whatever we have to do to move it from one side of the inner product to the other, i.e. whatever Â\* satisfies ⟨u,Âv⟩=⟨Â\*u,v⟩ for all u,v. (Omitted techicality: we must further restrict ourselves to functions that are sufficiently differentiable that ⟨u,Âu⟩ is finite, which is called a **Sobolev space** for this Â, a subset of the Hilbert space.) For matrices and ordinary vector dot products, this is equivalent to the "swap rows and columns" definition. For differential operators, it corresponds to integration by parts, and depends on the boundary conditions as well as on the operator and on the inner product.

Showed that with u(0)=u(L)=0 boundary conditions and this inner product. (d/dx)T\=-(d/dx)...very closely analogous to what we found with the finite-difference derivative matrix D. The "transpose of a derivative" is minus a derivative!

Furthermore, showed that (d2/dx2)T is real-symmetric (also called "Hermitian" or "self-adjoint").

Not only that, but showed that d2/dx2 is negative-definite on this space, since ⟨u,u''⟩=-∫|u'|2, and u'=0 only if u=constant=0 with these boundary conditions.

Showed that the proof of real eigenvalues from 18.06 carries over without modification for Hermitian operators; similarly for the proof of orthogonal eigenvectors, hence the orthogonality of the Fourier sine series. Similarly for the proof of negative eigenvalues.

So, many of the key properties of d2/dx2 follow "by inspection" once you learn how to transpose operators (integrate by parts). And this immediately tells us key properties of the solutions, if we assume the spectral theorem: Poisson's equation has a unique solution, the diffusion equation has decaying solutions, and the wave equation has oscillating solutions.

**Further reading:** [Notes on function spaces, Hermitian operators, and Fourier series](http://web.mit.edu/18.06/www/Fall07/operators.pdf) that I once wrote for 18.06 (slightly different notation). Textbook, section 3.1: transpose of a derivative. The general topic of linear algebra for functions leads to a subject called _functional analysis_; a rigorous introduction to functional analysis can be found in, for example, the book _Basic Classes of Linear Operators_ by Gohberg et al. There are some technicalities that I omit: a differential operator is only called "self-adjoint" if it is equal to its adjoint and is "densely defined", and showing that an operator equals its adjoint furthermore requires an extra step of showing that Â and Â\* act on the same domains.

### Lecture 5 (Sep. 13)

Discussed diagonalizability of infinite-dimensional Hermitian operators. Unlike the proof of real eigenvalues, etcetera, we cannot simply repeat the proof from the matrix case (where one can proceed by induction on the dimension). In practice, however, real-symmetric operators arising from physical systems are almost always diagonalizable; the precise conditions for this lead to the "spectral theorem" of functional analysis.) In 18.303, we will typically just assume that that all functions of interest lie in the span of the eigenfunctions, and focus on the consequences of this assumption. Diagonalizability and spectral theorems. Self-adjointness of more operators. Need for modified inner products. From discrete systems to PDEs.

Not only do we now understand d2/dx2 at a much deeper level, but you can obtain the same insights for many operators that _cannot_ be solved analytically. For example, showed that the operator d/dx \[c(x) d/dx\], which is the 1d Laplacian operator for a non-uniform "medium", is also real-symmetric positive definite if c(x)>0, given the same u(0)=u(L)=0 boundary conditions.

As another example, considered the operator c(x)d2/dx2 for real c(x)>0. This is _not_ self-adjoint under the usual inner product, but _is_ self-adjoint if we use the _modified_ inner product ⟨u,v⟩=∫uv/c with a "weight" 1/c(x). (This modified inner product satisfies all of our required inner-product properties for positive c(x).) Therefore, c(x)d2/dx2 indeed has real, negative eigenvalues, and has eigenfunctions that are orthogonal under this new inner product. Later on, we will see more examples of how sometimes you have to change the inner product in order to understand the self-adjointness of Â.

Fortunately, it's usually pretty obvious how to change the inner product, typically some simple weighting factor that falls out of the definition of Â. (In fact, for matrices, it turns out that _every_ diagonalizable matrix with real eigenvalues is Hermitian under some modified inner product. I didn't prove this, however.)

**Handout:** [handwritten notes](lecture-5.5.pdf)

Previously, we started with the continuous PDE equations and derived a discrete/matrix version as an approximation. Now we will do the reverse: we will _start_ with a truly discrete (finite-dimensional) system, and then derive the continuum PDE model as a limit or approximation. This is one of the common ways by which PDE models are derived in the first place, will give yet another perspective on similar mathematics, and will also shed some light on the origin of the variable c(x) coefficients from the last lecture.

In particular, we look at the 8.01 (or 8.03) system of a set of N masses m sliding without friction and connected by springs of constant k, anchored to walls on both sides.

If the displacements of each mass are given by a vector **u**, we show that d2**u**/dt2\=A**u** where A = (k/m) (-DTDT) Δx2, with D being _exactly_ our "first derivative" matrix from previous lectures. Thus, up to a scale factor, our discrete Laplacian A from previous lectures is not only an approximation for d2/dx2, but it is also _exactly_ the matrix describing this coupled-spring system!

(In the next lecture, we will cover the remainder of the notes.)

**Further reading:** Same as previous lecture. Sections 2.1 and 2.2 in the Strang book cover very similar material on discrete vibration problems.

### Lecture 6 (Sep. 16)

Continued with lecture 5.5 notes from last lecture...

By inspection, A is real-symmetric negative-definite as before, and hence we have N real eigenvalues λn<0 and orthonormal eigenvectors **u**n. By expanding the solution **u**(t) in terms of these eigenvectors, we see that we obtain oscillating solutions: a set of (orthogonal) _normal modes_ with real "eigenfrequencies" ωn\=√-λn. The negative-definiteness is critical to have oscillation, as otherwise we would get complex ωn and exponential growth! Showed a few examples for N=1,2,3 to get an intuition for these modes.

Took the N→∞ limit keeping the total length L fixed, which corresponds to Δx→0 (breaking the system into smaller and smaller pieces). Correspondingly, we decrease the masses proportional to Δx, so that m=ρΔx where ρ is a density. On the other hand, reminded ourselves that cutting springs in half _increases_ the spring constant, so we should let k=c/Δx for some constant c. With these definitions, our matrix A is precisely (c/ρ) (-DTD), where -DTD is our approximate (center-difference) d2/dx2 from before, and hence this limit gives precisely the scalar wave equation ∂2u/∂t2\=(c/ρ) ∂2u/∂x2, with the boundary conditions u(0,t)=u(L,t)=0 (fixed ends). As before, this operator is self-adjoint negative definite and so we get real λn<0 etcetera. Exactly as in the finite-dimensional case above, we therefore get oscillating "normal mode" solutions, just with an infinite series of eigenfunctions instead of a finite sum.

Finally, considered the "inhomogeneous medium" case where we allow all the masses and spring constants to be different. Showed that this corresponds to inserting some diagonal matrices into the above expressions, and in the continuum limit gives ∂2u/∂t2\=(1/ρ) ∂/∂x (c ∂u/∂x)=Âu where ρ(x) and c(x) are positive functions. As in the previous lecture, we can see that this Â is indeed self-adjoint negative-definite, and hence we get _oscillating normal-mode solutions_, if we define the modified inner product ⟨u,v⟩=∫ρuv. And we can see this even though we probably cannot solve this PDE by hand except for very special cases of ρ and c!

**Handout:** [handwritten notes](lecture-6.pdf)

Now that we have seen several specific examples, we are equipped to consider more general problems, ones that are even harder to solve analytically.

Started in 1d with the "Sturm-Liouville operator" w(x)\-1 \[ - d/dx (c(x) d/dx) + p(x) \], with Dirichlet (0) boundaryies on Ω=\[0,L\]. Showed that it is self-adjoint under the weighted inner product ⟨u,v⟩=∫wuv, assuming w is real and positive and c and p are real. If we further assume that c≥0 and p>0, the operator is positive-definite, a so-called **elliptic operator** (although the technical definition of elliptic operators is slightly more restrictive in order to exclude pathological cases of the coefficient functions).

**Further reading:** sections 2.1 and 2.2 in the Strang book cover very similar material on the ball/spring system. Much of the theory of these kinds of 1d operators is traditionally called "Sturm-Liouville theory", and can be found under that name in many books (e.g. _Methods of Applied Mathematics_ by Hildebrand, _Methods of Mathematical Physics_ by Courant & Hilbert, and many similar titles). Even Wikipedia has a [decent article](http://en.wikipedia.org/wiki/Sturm%E2%80%93Liouville_theory) on the topic under that name.

### Lecture 7 (Sep. 18)

**Handouts:** [notes on elliptic operators](lecture-6.pdf) (from last time), [notes on separation of variables](separation.pdf), [pset 2 solutions](pset2sol-f13.pdf) (see also [IJulia notebook](pset2-f13.ipynb); [html format](pset2-f13.html)), [pset 3](pset3-f13.pdf) (due Sep 25)

Generalized Sturm-Liouville operators to multiple dimensions: Â=w(**x**)\-1 \[ - ∇⋅(c(**x**)∇) + p(**x**) \], with Dirichlet (0) boundaryies on some finite domain Ω, and again we will show that this is self-adjoint for real coefficients and w>0, and positive-definite (elliptic) for c≥0 and p>0. The key to the proof is the divergence theorem, which generalizes "integration by parts" to multiple dimensions.

We can now analyze three important cases, and give them their conventional historical names:

*   Âu=f where Â is self-adjoint and positive (or negative) definite. This is sometimes called an _elliptic_ problem.
*   Âu=∂u/∂t where Â is self-adjoint and negative definite (hence exponentially decaying solutions), or possibly semidefinite. This is sometimes called a _parabolic_ problem.
*   Âu=∂u2/∂t2 where Â is self-adjoint and negative definite (hence oscillating solutions), or possibly semidefinite. This is sometimes called a _hyperbolic_ problem.

(I won't spend any times on the analogies between these equations and those of parabolas, ellipses, and hyperbolas. That only works well in simple cases with scalar functions u, and I find it much clearer and more general just to talk about the definiteness and self-adjointness of Â.)

For fun, spent a couple of minutes on the Schrödinger wave equation, which is a different scalar wave equation that we will come back to later.

### Lecture 8 (Sep. 23)

**New topic: Separation of variables:** (See [notes](separation.pdf).) This is a technique to _reduce the dimensionality_ of a PDE by representing the solution as a product of lower-dimensional functions. It _only works in a handful of cases_, usually as a consequence of _symmetry_, but those cases are common enough that it is important to know them. It also gives us our only analytically solvable PDE examples in more than 1d; otherwise we will have to use the computer.

**Separation of Time**: The most important case is the one we've already done, under another name. We solved Au=∂u/∂t by looking for eigenfunctions Au=λu, and then multiplying by exp(λt) to get the time dependence. Similarly for Au=∂u2/∂t2 except with sines and cosines. In both cases, we wrote the solution as a sum of products of purely spatial functions (the eigenfunctions) and purely temporal functions like exp(λt). The key point here is that we aren't assuming that the _solution_ is separable, only that it can be decomposed into _linear combination_ of separable functions.

**Separation of Space**: Here, we try to solve problems in more than one _spatial_ dimension by factoring out 1d problems in one or more dimension. In particular, we will try to find _eigenfunctions_ in separable form, and then write any solution as a linear combination of eigenfunctions as usual. In practice, this mainly works only in a few important cases, especially when one direction is _translationally invariant_ or when the problem is _rotationally invariant_. In the former case, translational invariance in one direction (say z) allows us to write the eigenfunctions in separable form as X(x,y)Z(z), where it turns out that Z(z)=exp(ikz) for some k (and X and λ will then depend on k). In the latter case, we get separable eigenfunctions R(r)exp(imθ) where m is an integer, in 2d, and R(r)Yl,m(θ,φ) in 3d, where Yl,m(θ,φ) is a [spherical harmonic](http://en.wikipedia.org/wiki/Spherical_harmonics). Also, we can _sometimes_ get separable solutions for finite "box-like" domains, i.e. translationally invariant problems that have been truncated to a finite length in z.

To start with, we looked at ∇2u=λu in a 2d Lx×Ly box with Dirichlet boundary conditions, and looked for separable solutions of the form X(x)Y(y). Plugging this in and dividing by XY (the standard techniques), we get 1d eigenproblems for X and Y, and these eigenproblems (X''=X×constant and Y''=Y×constant) just give us our familiar sine and cosine solutions. Adding in the boundary condition, we get sin(nxπx/Lx) sin(nyπx/Ly) eigenfunctions with eigenvalues λ=-(nxπ/Lx)2\-(nyπ/Ly)2. As expected, these are real and negative, and the eigenfunctions are orthogonal...giving us a 2d Fourier sine series. For example, this gives us the "normal modes" of a square drum surface.

Finished consideration of separability of ∇2u=λu in a 2d box, from notes: discussed orthogonality of these eigenfunctions. Also showed that separability breaks down, in general, for non-constant coefficients in the box.

### Lecture 9 (Sep 25)

**Handouts:** [notes on separation of variables](separation.pdf) (from last time), [notes on Bessel functions and cylindrical problems](lecture-8.pdf), [pset 3 solutions](pset3sol-f13.pdf), [pset 4](pset4-f13.pdf) (due Oct 2); also [IJulia Bessel-function notebook](lecture-8.ipynb) (also in [html](lecture-8.html))

More separation of variables: cylindrical case of a cylinder of radius R with Dirichlet boundary conditions. Show that the Laplace eigenequation here is indeed separable into a function of θ multiplied by a function of r, satisfying separate 1d ODEs. Show that the θ dependence is sin(mθ) or cos(mθ) (or any linear combination), where m is an integer (in order to be periodic in θ). The r dependence satisfies a more complicated 2nd-order ODE that we can't solve by hand (even if you have taken 18.03).

At this point, it's more a historical question than a mathematical one: has someone solved this equation before, and if so is there a standard name (and normalization, etc) for the solutions? In fact, that is the case here (not surprisingly, since the Laplacian is so important): our r equation is an instance of **Bessel's equation**, and the solutions are called **Bessel functions**. The canonical two Bessel functions are Jm and Ym: there is a standard convention defining the normalization, etcetera, of these, but the important thing for our purposes is that J is finite at r=0 and Y blows up at the origin. In Matlab, these are supplied as the built-in functions [besselj](http://www.mathworks.com/help/techdoc/ref/besselj.html) and [bessely](http://www.mathworks.com/help/techdoc/ref/bessely.html), and we use Matlab to plot a few of them to get a feel for what they look like: basically, sinusoidal functions that are slowly decaying in r.

To get eigenfunctions, we have to impose boundary conditions. Finite-ness of the solution at r=0 means that we can only have Jm(kr) solutions, and vanishing at r=R means that kR must be a root of Jm. We have to find these roots numerically, but this is easy to do, and we obtain a discrete set of eigenfunctions and eigenvalues.

From the general orthogonality of the Laplacian eigenfunctions, we can derive an orthogonality relation for Bessel functions, and by evaluating the integral numerically we can see that this orthogonality is indeed the case.

By looking at Bessel's equation asymptotically, we find that it reduces to sines and cosines for large r; more careful considerations show that it must actually reduce to sines and cosines multiplied by 1/√r, and we can verify this from the plot. Conversely, for small r we show that it goes as either rm (Jm) or 1/rm (Ym, except for m=1 where Y0 is proportional to log r); this is why we have one finite solution and one divergent one at r=0. (There are many, many more properties of Bessel functions that one can derive analytically, but that is not our major concern here.)

**Further reading:** The Wikipedia page on [Bessel functions](http://en.wikipedia.org/wiki/Bessel_function) has many plots, definitions, and properties, as well as links to standard reference works.

### Lecture 10 (Sep 27)

Discussed **boundary conditions** more generally than we have done in the past. Up to now, we have mostly considered u=0 (Dirichlet) or **n**⋅∇u=0 (Neumann) on the boundary, and mostly the former.

More generally, we can consider "general Dirichlet" and "general Neumann" boundary conditions, where either the values u(x) or the normal derivatives **n**⋅∇u are some given function g(**x**,t) on the boundary. (For example, general Dirichlet boundary conditions arise for a drum head where the edges are not held flat, and you may even be warping the edges as a function of time.) If g≠0, these functions u do _not_ form a vector space because they do not include u=0, so we must transform the problem somehow to get back to linear algebra.

For general Dirichlet, one simple approach is to write u=v+g, where our new unknown v has zero-Dirichlet boundaries (similar to pset 1). Showed how this transforms e.g. a wave equation Âu=∂2u/∂t2\-f(x,t) into wave equation in v but with an additional "force" term modifying f: Âv=∂2v/∂t2\-\[f+Âg-∂2g/∂t2\]. For example, considered the steady-state 1d version of this problem d2u/dx2\=0 (Laplace's equation) with u(0)=a, u(L)=b boundaries and showed that the solution is a straight line, which is physically obvious for a string stretched from (0,a) to (L,b).

Intuitively, this makes a certain amount of sense: warping the boundary corresponds to an external force. But intuitively, the "physical" boundary force is only applied at the boundary, not everywhere in Ω as it is for a general g(x) above. It turns out that we can do this, too. It is easier to see this in the discrete case, for the same 1d problem as above. In this case, showed that we obtained the same Dirichlet A matrix as we do for 0 boundary conditions, while the (a,b) boundary conditions just turned into terms added to the right hand side, but only in the first and last rows: an "external force" applied at the boundaries. The PDE version of this technique involves delta functions, which we aren't prepared to handle yet. (In fact, this generalizes to cases where we want to specify jumps and other discontinuities of u in the _interior_ of Ω as well, in which case one can again use new surface-localized terms on the right-hand-side and it is sometimes called an "immersed boundary" or "imbedded boundary" method, especially in fluid mechanics.)

To better understand how Neumann boundary conditions arise, we have to better understand the meaning of ∇u. Considered the case of a diffusion equation, where u=mass/volume of some solute, and ∇⋅c∇u=∂u/∂t. The total mass M within some volume V is just ∫Vu, and showed by applying the divergence theorem we obtain dM/dt equal to a surface integral of c∇u. Since dM/dt>0 when mass is flowing _in_ to the volume, this means that -c∇u is a mass "flux" vector (mass/time⋅area).

If we have diffusion in a closed container, so that no mass can flow in or out of Ω, we then immediately see that we should apply (0) Neumann boundary conditions. Furthmore, total mass = ∫Ωu is _conserved_ (constant in time) for any solution u.

More generally, for any equation Âu=∂u/∂t, showed that we obtain a **conservation law** ∂/∂t ⟨v,u⟩=0 for any v(**x**) in the _left null space_ N(Â\*).

For the case of diffusion with Neumann boundary conditions, reviewed the fact that Â=Â\* but Â is only negative _semidefinite_: N(Â)=N(Â\*) contains any _constant function_, and is spanned by v(**x**)=1. Hence ⟨1,u⟩ is conserved. i.e. total mass, or total heat, or average temperature, is conserved in a closed/insulated Ω.

Another example of a (0) Neumann boundary condition arises when we are considering u(**x**) that are mirror-symmetric (even) around some mirror plane, which is equivalent to imposing a Neumann boundary condition on the mirror plane. (Similar, antisymmetric/odd symmetry is equivalent to a zero Dirichlet boundary.) Another example is a stretched string where one end can slide freely up and down a rod with no friction: that end has a Neumann condition.

There are many other possible boundary conditions, of course. The most complicated ones can arise for PDEs with multiple unknowns (e.g. pressure, temperature, velocity, ...), in which case the boundary conditions may be equations relating several different unknowns or their derivatives.

One can also have _nonlocal_ boundary conditions, in which u at one point on ∂Ω is related to u at a _different_ point. The most common example of this are _periodic_ boundary conditions. e.g. considered Â=d2/dx2 on \[0,L\] for u(0)=u(L). Showed that Â is still self-adjoint, but not because the boundary terms are _individually_ zero, but rather because the x=0 and x=L boundary terms _cancel_. The eigenfunctions are now sines _and_ cosines of 2πnx/L, and give a general Fourier series (not just a sine or cosine series)! Also, Â is now negative _semidefinite_ because constant u are allowed. Hence, for example, diffusion on a periodic domain still conserves total mass, because any mass that exits one side comes back in through the other side.

**Further reading:** The u=v+g trick is closely related to the standard proof of the uniqueness of solutions to Laplace's/Poisson's equation with general Dirichlet boundaries (google "Laplace uniqueness" or "Poisson uniqueness", e.g. [this page](http://www-solar.mcs.st-and.ac.uk/~andy/LectureNotes/Fundamentals1/node52.html)). The trick of moving boundary conditions over to the right-hand side is so obvious for finite-difference methods that it hardly has a name, but it is often commented on explicitly for finite-element methods where things are less obvious (e.g. section 3.6 of the book). There is a [review of immersed boundary methods](www.stanford.edu/group/uq/pdfs/journals/annurev_05.pdf) by Mittal and Iaccarino that is fairly readable, but oriented mainly towards fluid mechanics. Periodic domains arise in many cases, the most obvious being equations on a torus (e.g. waves on a membrane that loops back to itself, diffusion in a circular tube, or masses and springs connected into a ring). They also arise for systems that repeat periodically, e.g. a periodic crystal in solid-state physics, in which case you can write the solutions as [Bloch waves](http://en.wikipedia.org/wiki/Bloch_wave) of the form u(**x**)=u**k**(**x**)ei**k**⋅**x** where u**k** is a periodic function that solves a PDE with periodic boundary conditions (and plotting the eigenvalues as a function of **k** gives a [band structure](http://en.wikipedia.org/wiki/Electronic_band_structure)).

### Lecture 11 (Sep 30)

**Handouts:** [handwritten notes](lecture-10.pdf)

2d finite-difference discretizations: discretized the 2d Laplacian operator ∇2 in an Lx×Ly box with Nx×Ny points, for Dirichlet (0) boundaries, so that u(mΔx,nΔy) is approximated by NxNy degrees of freedom um,n. Showed that simply applying the 1d center-difference rule along x and y results in a (famous) "[5-point stencil](Five-point stencil)" approximation for -∇2 in which the Laplacian at (nx,ny) depends on u at (nx,ny) and the 4 nearest-neighbor points.

In order to express this as a matrix A, however, we need to "flatten" the 2d grid of points unx,ny into a single column vector **u** with NxNy components. There are multiple ways to do this, but a standard and easy scheme is the "[column-major](http://en.wikipedia.org/wiki/Row-major_order)" approach in which **u** is formed from the contiguous columns (x) unx,: concatenated in sequence. (This is the approach used internally within Matlab to store matrices.)

Given this, we then see how to form ∂2/∂x2 by operating one Nx\-column at a time, using the the Nx×Nx discrete 1d Laplacian Ax (=-DxTDx). The ∂2/∂x2 matrix is simply a matrix with Ax along the diagonal Ny times, which differentiates each Nx\-column block by block. The ∂2/∂y2 matrix is a little tricker, but if we think of operating on whole columns then we see that it is just the Ay matrix with the entries "multiplied" by Nx×Nx identity matrices Ix.

In order to have a convenient way to express this, we use the [Kronecker product](http://en.wikipedia.org/wiki/Kronecker_product) notation A⊗B \[[kron](http://www.mathworks.com/help/techdoc/ref/kron.html)(A,B) in Matlab\], which multiplies the _entries_ of A by the _matrix_ B to create a _matrix of matrices_. In this notation, A = Iy⊗Ax + Ay⊗Ix.

Using this machinery, constructed A for Nx\=10 and Ny\=15 for Lx\=1 and Ly\=1.5 in Matlab. Visualized the pattern of nonzero entries with [spy](http://www.mathworks.com/help/techdoc/ref/spy.html). Solved for the eigenfunctions, and plotted a few; to convert a column vector **u** back into a 2d matrix, used [reshape](http://www.mathworks.com/help/techdoc/ref/reshape.html)(**u**,Nx,Ny), and plotted in 3d with the [surf](http://www.mathworks.com/help/techdoc/ref/surf.html) command. The first few eigenfunctions can be seen to roughly match the sin(nxπx/Lx) sin(nyπx/Ly) functions we expect from separation of variables. However, Nx\=10, Ny\=15 is rather coarse, too coarse a discretization to have a really nice (or accurate) picture of the solutions.

In order to increase Nx and Ny, however, we have a problem. If the problem has N=NxNy degrees of freedom, we need to store N2 numbers (8N2 bytes) just to store the matrix A, and even just solving Ax=b by Gaussian elimination takes about N3 arithmetic operations. Worked through a few numbers to see that even Nx\=Ny\=100 would have us waiting for 20 minutes and needing a GB of storage, while 3d grids (e.g. 100×100×100) seem completely out of reach. The saving grace, however, is [sparsity](http://en.wikipedia.org/wiki/Sparse_matrix): the matrix is mostly zero (and in fact the 5-point stencil A has < 5N nonzero entries). This means that, first, you can store only the nonzero entries, greatly reducing storage. Second, it turns out there are ways to exploit the sparsity to solve Ax=b much more quickly, and there are also quick ways to find a _few_ of the eigenvalues and eigenvectors.

In Julia, you exploit sparsity by using the `sparse` command and friends to create sparse matrice. Once you have a sparse matrix, Matlab automatically uses algorithms to exploit sparsity if you solve Ax=b by x=A\\b and use the `eigs` function to find a few eigenvalues (instead of `eig`).

Starting with the ∇2 operator on a square grid, showed how we can convert to any other Ω shape with Dirichlet boundaries just by taking a subset of the rows/cols (as in problem 2 of pset 5). Recovered the Bessel solutions for a circular domain. See the IJulia notebook:

*   [lecture-10.ipynb](lecture-10.ipynb) (also in [HTML](lecture-10.html))

**Further reading**: Section 3.5 of the Strang book on 2d finite differences, section 7.1 on sparsity. See, for example [min-max theorem in Wikipedia](http://en.wikipedia.org/wiki/Min-max_theorem), although this presentation is rather formal. Unfortunately, most of the discussion you will find of this principle online and in textbooks is either (a) full of formal functional analysis or (b) specific to quantum mechanics \[where the operator is A=-∇2+V for some "potential-energy" function V(x)\]. You can find another [Laplacian demo here](http://facstaff.unca.edu/mcmcclur/class/LinearII/presentations/html/2.08.02.TwoDVibes.html).

### Lecture 12: October 2

**Handouts:** [notes on the min–max theorem](minmax.pdf), [pset 5](pset5-f13.pdf), [pset 4 solutions](pset4sol-f13.pdf) and [notebook](pset4-f13.ipynb) (also [HTML](pset4-f13.html))

Starting with the ∇2 operator on a square grid (from last lecture), showed how we can convert to any other Ω shape with Dirichlet boundaries just by taking a subset of the rows/cols. Looked at a couple of triangular domains, and recovered the Bessel solutions for a circular domain.

In order to get an intuitive feel for what the eigenfunctions should look like, a powerful tool is the **min–max theorem**. See handout for notes.

As a final example corresponding to the -c∇2 operator in the notes, considered an "L"-shaped domain Ω with c=1/w(x). In particular, suppose that w(x)=1 everywhere except for a small region where w(x)=w0\>1. In order to concentrate in this small region, u(x) will have to have bigger slope (sacrificing the numerator). As w0 increases, we expect the denominator of the Rayleigh quotient to "win" and the concentration to increase, while for w0 close to 1 the eigenfunctions should be similar to the case of -∇2.

**Further reading:** See, for example [min-max theorem in Wikipedia](http://en.wikipedia.org/wiki/Min-max_theorem), although this presentation is rather formal. Unfortunately, most of the discussion you will find of this principle online and in textbooks is either (a) full of formal functional analysis or (b) specific to quantum mechanics \[where the operator is A=-∇2+V for some "potential-energy" function V(x)\]. You can find another [Laplacian demo here](http://facstaff.unca.edu/mcmcclur/class/LinearII/presentations/html/2.08.02.TwoDVibes.html).

### Lecture 13: October 4

**Handouts:** [introduction to Green's functions](green.pdf)

(See notes.) Introduced Green's functions by analogy with matrix inverses, and constructed Green's function of -d2/dx2 with Dirichlet boundaries as an example.

We had to jump through some hoops to avoid a problematic-looking "delta function" that keeps appearing, a limit of a function whose area is "infinitely concentrated" at a "single point". This is possible, but becomes more and more painful as we go on, motivating us to find an alternate definition of "function" in the future, a **distribution**.

**Further reading:** Strang book, section 1.4. Many PDE books introduce Green's functions and delta functions in various ways; see, e.g. section 9.3.4 of _Elementary Applied Partial Differential Equations_ by Haberman.

### Lecture 14: October 7

**Handouts:** [notes on reciprocity and positivity of Green's functions](reciprocity.pdf), [explicit check of 1d Green's function solution](Green-explicit.pdf), [notes on delta functio ns and distributions](delta-notes.pdf)

Reviewed definition and construction of Green's function via limits of narrow "box" sources from last lecture.

For the 1d example, we can explicitly check that u(x)=∫G(x,x')f(x')dx' solves -u''=f. (See 2nd handout.)

Went through [notes](reciprocity.pdf) on reciprocity and positivity of Green's functions.

This is all extremely cumbersome though—we have to go through lots of contortions to avoid differentiating discontinuities and avoid delta functions. More generally, there are a number of difficulties that continually arise when we deal with classical functions (mapping numbers to numbers): went through section 1 of the distribution handout.

### Lecture 15: October 9

Delta functions and distributions: finished [notes](delta-notes.pdf) from previous lecture.

**Further reading:** See the books _Generalized Functions_ by Gel'fand and Shilov and _A Guide to Distribution Theory and Fourier Transforms_ by Strichartz referenced at the end of the notes. Wikipedia has a decent article on [distributions](http://en.wikipedia.org/wiki/Distribution_%28mathematics%29). The idea that there are functions φ(x) which are infinitely differentiable but are zero outside of a finite region is a bit counterintuitive if you think about the interface between the zero and nonzero regions, but it is quite possible; see [bump function](http://en.wikipedia.org/wiki/Bump_function) on Wikipedia for an elaboration on the example I gave in class, and a proof that the derivatives are continuous [here](http://en.wikipedia.org/wiki/Non-analytic_smooth_function). In practice, however, we will almost never have to explicitly construct test functions to talk about distributions.

### Lecture 16: October 11

**Handouts:** [pset 5 solutions](pset5sol-f13.pdf) and [solutions notebook](pset5-f13.ipynb) (also in [HTML](pset5-f13.html)), [pset 6](pset6-f13.pdf) (due next Wed.)

Derived Green's function of ∇2 in 3d for infinite space (requiring solutions to → zero at infinity to get a unique solution), in three steps:

1.  Because the ∇2 operator is invariant under translations (changes of variables **x**→**x**+**y**), showed that G(**x**,**x**') can be written as G(**x**,**x**')=G(**x**\-**x**',0). Similarly, rotational invariance implies that G(**x**\-**x**',0)=g(|**x**\-**x**'|) for some function g(r) that only depends on the distance from **x**'.
2.  In spherical coordinates, solved -∇2g = 0 for r > 0 (away from the delta function), obtaining g(r)=c/r for some constant c to be determined.
3.  Took the distributional derivative (-∇2g){φ}=g{-∇2φ} ("integrating by parts" using the fact from Lecture 7 that ∇2 is self-adjoint) for an arbitrary test function φ(**x**), and showed by explicit integration that we get cφ(0). Therefore c=1/4π for us to solve -∇2g = δ(**x**\-**x**').

Hence G(**x**,**x**') = 1/4π|**x**\-**x**'| for this problem, and -∇2u=f is solved by u(**x**)=∫f(**x**')d3**x**'/4π|**x**\-**x**'|.

A physical example of this can be found in electrostatics, from 8.02: the potential V of a charge density ρ, satisfies -∇2V=ρ/ε0. A point charge q at **x**' is a charge density that is zero everywhere except for **x**', and has integral q, hence is ρ(**x**)=qδ(**x**\-**x**'). Solving for V is exactly our Green's function equation except that we multiply by q/ε0, and hence the solution is V(**x**) = q/4πε0|**x**\-**x**'|, which should be familiar from 8.02. Hence -∇2V=ρ/ε0 is solved by V(**x**)=∫ρ(**x**')d3**x**'/4πε0|**x**\-**x**'|, referred to in 8.02 as a "superposition" principle (writing any charge distribution as the sum of a bunch of point charges).

Perhaps the most important reason to solve for G(**x**,**x**') in empty space is that solutions for more complicated systems, with boundaries, are "built out of" this one.

An illustrative example is Ω given by the 3d half-space z>0, with Dirichlet boundaries (solutions=0 at z=0). For a point **x**' in Ω, showed that the Green's function G(**x**,**x**') of -∇2 is G(**x**,**x**')=(1/|**x**\-**x**'| - 1/|**x**\-**x**''|)/4π, where **x**'' is the same as **x**' but with the sign of the z component flipped. That is, the solution in the upper half-space z>0 looks like the solution from _two_ point sources δ(**x**\-**x**')-δ(**x**\-**x**''), where the second source is a "negative image" source in z<0. This is called the **method of images**.

Reviewed method-of-images solution for half-space. There are a couple of other special geometries where a method-of-images gives a simple analytical solution, but it is not a very general method ([complicated generalizations](http://www2.imperial.ac.uk/~dgcrowdy/PubFiles/Paper-20.pdf) for 2d problems notwithstanding). The reason we are covering it, instead, is that it gives an analytically solvable example of a principle that _is_ general: Green's functions (and other solutions) in complicated domains _look like solutions in the unbounded domain plus extra sources on the boundaries_.

**Further reading**: See e.g. sections 9.5.6–9.5.8 of _Elementary Applied Partial Differential Equations_ by Haberman for a traditional textbook treatment of Green's functions of ∇2 in empty space and the half-space. If you Google "method of images" you will find lots of links, mostly from the electrostatics viewpoint (see e.g. [these lecture notes](http://www.phys.ufl.edu/~dorsey/phy6346-00/lectures/lect04.pdf)); see also e.g. _Introduction to Electrodynamics_ by Griffiths for a standard textbook treatment; the only mathematical difference introduced by (vacuum) electrostatics is the multiplication by the [physical constant](http://en.wikipedia.org/wiki/Vacuum_permittivity) ε0 (and the identification of -∇V as the electric field).

### Lecture 17: October 16

**Handouts:** [pset 6 solutions](pset6sol-f13.pdf) and [notebook](pset6-f13.ipynb) (also in [HTML](pset6-f13.html)), [pset 7](pset7-f13.pdf) [notes on Green's functions in inhomogeneous media](inhomog-notes.pdf)

In the image method, the "extra source" is ostensibly not on the boundary, it is on the other side of the boundary. However, we can transform it to what we want by the following trick: consider the function u(**x**) in Ω=R3 that equals (1/|**x**\-**x**'| - 1/|**x**\-**x**''|)/4π \[the method-of-images solution\] for z>0 and u(**x**)=0 for z<0. What right-hand-side does -∇2u give? In z>0 -∇2u gives δ(**x**\-**x**') as before, and for z<0 -∇2u gives zero. _At_ z=0, however, there is a slope discontinuity in (1/|**x**\-**x**'| - 1/|**x**\-**x**''|)/4π, which means that -∇2u also gives a δ(z) term: δ(z) σ(x,y) for a σ(x,y) given by the amplitude of the slope discontinuity.

What does this mean? Our solution u(**x**) is due to the sum of a point source at **x**' and _sources at the interface_ (z=0). Worked out what these sources σ(x,y) are. Physically, in the electrostatics example they correspond to a _surface charge density_ on the surface of a conductor. Why are these sources there? They are there to _cancel_ the effect of the source at **x**' for z<0, enforcing the boundary condition u=0 at z=0.

More generally, we can do this for _any_ interface dΩ: we can write the solution from a point source δ(**x**\-**x**') in Ω as the sum of the solution from that point source plus an integral of _unknown point sources_ σ(**x**') for points **x**' the boundary dΩ. Formally, we determine σ(**x**') by requiring u(**x**) to satisy the boundary condition at dΩ, which gives a _surface integral equation_ (SIE) (of the "first kind") for σ(**x**'). Numerically, we discretize the surface in some way to get a finite number of unknowns approximating σ(**x**'), leading to an SIE numerical method.

SIE methods (most commonly the "boundary element method", BEM) are very powerful in several ways. Because they only have unknowns on the _boundaries_ (not everywhere in space like in a PDE method like finite differences), they can greatly reduce the number of unknowns: they handle the homogeneous regions analytically. They can handle infinite space (e.g. a cylinder surrounded by infinite space as in the pset) analytically, with no artificial truncation. Done properly, the matrices can have very nice properties. There are also some difficulties. SIE methods are not so efficient for problems that are not mostly homogeneous, especially continuously-varying media, and nonlinear or time-dependent problems are also difficult. Because the empty-space Green's function (1/4π|**x**\-**x**'| in 3d) blows up for nearby points, there are lots of tricky singularities that must be dealt with carefully in setting up SIE methods. Furthermore, because you have long-range interactions (every surface point interacts with every other point via the Green's function), the matrices are dense, not sparse. That means that developing fast solvers for large problems is tricky; remarkably, there are ways to do this (most famously the pioneering [fast multipole method](http://en.wikipedia.org/wiki/Fast_multipole_method) invented in 1985), but implementing them is not for the timid. Worse, the singularity-handling and fast-solver techniques depend heavily on the specific Green's function of empty space; for example, changing from 3d (1/|**x**\-**x**'|) to 2d (ln|**x**\-**x**'|) problems requires a completely different implementation, although the concepts are similar.

(Finished notes through 2.2.1.)

**Further reading:** There are many books on integral-equation methods, e.g. _Boundary Integral Equation Methods for Solids and Fluids_ by Bonnet is a reasonably general introduction.

### Lecture 18: October 18

Finished notes from previous lecture.

### Lecture 19: October 21

**Handouts:** this summary

Handout from previous lecture, section 3.1.

**New topic: Time-stepping and stability.** Before, we turned operator equations Âu=f into matrix equations A**u**\=**f** by discretizing in space. Now, we want to turn time-dependent operator equations Âu=f+∂u/∂t into discrete equations in both time and space. This will involve a new concern: **stability**.

Began with a trivial example of an operator Â=a, a single number a<0, which for f=0 gives the ODE du/dt=au, and has the exponentially decaying (for a<0) solution u(t)=u(0)eat. Now we will discretize u(t) in time by u(nΔt)≈un — we will always use _superscripts_ to denote the _timestep_ n. Approximating the time derivative by a _forward difference_ yields un+1≈(1+aΔt)un\=(1+aΔt)n+1u0. Even though the exact ODE has decaying solutions, this discretization may have _exponentially growing_ solutions unless Δt<2/|a|: the discretization is only **conditionally stable**. In contrast, a _backward difference_ yields un+1≈(1-aΔt)\-1un\=(1+aΔt)\-1-nu0, which is always exponentially decaying for a<0: the scheme is **unconditionally stable**.

For a more general operator Â, we proceed conceptually in two steps. First, we discretize in space only to yield a system of ODEs A**u**\=∂**u**/∂t for a matrix A (e.g. by finite differences in space). Then we discretize in time and ask what happens to the eigenvectors of A. Focused on the case where A (and Â) are self-adjoint and negative-definite (negative eigenvalues λ<0), as for the heat equation (Â=∇2) with Dirichlet boundaries. In this case, showed that forward differences give an **explicit timestep** un+1≈(1+AΔt)un and are conditionally stable: we must have Δt<2/|λ|. In contrast, backward differences give an **implicit timestep** un+1≈(1-AΔt)\-1un where we must solve a linear system at each step, but are unconditionally stable (decaying for any Δt).

Some definitions:

*   Âu=∂u/∂t is **well posed** if the solution u(**x**,t) is finite for any finite t and for any initial condition u(**x**,0). (Note that PDEs with diverging solutions can still be well-posed, as long as they are finite at finite times, even if they are exponentially large.)
*   A discretization is **consistent** if the discretization goes to Âu=∂u/∂t as Δx and Δt → 0.
    *   If the difference between the discrete equations and Âu=∂u/∂t, the **local truncation error**, goes to zero as Δxa and as Δtb, then we say the scheme is "a-th order in space and b-th order in time."
*   A discretization is **stable** if **u**t/Δt≈u(**x**,t) does _not_ blow up as Δt → 0 and Δx → 0. (Informally, we often say it is "stable" if the solution does not blow up for any Δt, but a more precise definition has to take into account that the original PDE may have solutions that blow up as t→∞.)
    *   it is **conditionally stable** if it is stable only when Δt has a certain relationship to the spatial discretization A, and in particular this usually means that Δt is constrained by some relationship with Δx.
    *   it is **unconditionally stable** if it is stable for all Δt independent of Δx or A (or at least as long as A has some property like negative-definiteness).
*   A discretization is **convergent** if **u**t/Δt→u(**x**,t) as Δx, Δt → 0.

A very important result (stated here without proof) is the **Lax equivalence theorem**: for any consistent discretization of a well-posed linear initial-value problem, **stability implies convergence and vice versa**. If it is unstable, then it is obvious that it cannot converge: the discretization blows up but the real solution doesn't. Less obvious is the fact that _if it does not blow up, it must converge_.

The Lax theorem is very reassuring, because it turns out that it is quite difficult to prove stability in general (we usually prove necessary but not sufficient conditions in conditionally stable schemes), but if you run it and it doesn't blow up, you know it must be converging to the correct result.

The tricky case to analyze is that of conditionally stable schemes. We need to relate the eigenvalues of A to Δx in some way to obtain a useful condition on Δt.

For explicit timestepping of the heat/diffusion equation with forward differences, Δt is proportional to Δx2, so even though the discretization is second-order in space (errors ~ Δx2) and first-order in time (errors ~ Δt), the time and space discretization errors are comparable (or at least proportional).

On the other hand, for implicit timestepping with backward differences, Δt is independent of Δx, so the first-order accuracy in time can really limit us. Instead, presented a second-order scheme in time by considering (**u**n+1\-**u**n)/Δt to be a _center_ difference around step n+0.5 \[t=(n+0.5)Δt\]. In this case, we evaluate the right-hand side A**u** at n+0.5 by averaging: A(**u**n+1+**u**n)/2. This gives a **Crank-Nicolson** scheme: **u**n+1\=(1-AΔt/2)\-1(1+AΔt/2)**u**n. This is an implicit scheme, but is second-order accurate in both space and time (assuming a 2nd-order A in space). Showed that it is unconditionally stable if A is negative-definite.

For conditionally stable schemes, we need the eigenvalues of A. Gave a crude argument that the biggest |λ| for ∇2 and similar operators is proportional to Δx2, based on the fact that the solution cannot oscillate faster than the grid. To do better than this, we need to consider simplified cases that we can analyze analytically.

### Lecture 20: October 23

**Von Neumann analysis**. The idea of Von Neumann analysis is to analyze the eigenvalues of the space discretization, A, in a simple case that can be solved analytically: ∞ space and constant coefficients. In this case the eigensolutions will be sinusoids (Fourier modes), which are most conveniently written as complex exponentials.

In particular, considered Â=d2/dx2 in one dimension, discretized by the usual center difference. We try a solution um\=eikm, and show that it is indeed an eigenvector (with infinitely many components in ∞ space!) of the discretized second derivative. (Briefly reviewed properties of complex exponentials and the equivalence to sines and cosines by [Euler's identity](http://en.wikipedia.org/wiki/Euler%27s_formula).) Showed that the corresponding eigenvalues are λ(k)=−4 sin2(k/2) / Δx2. Hence, the maximum |λ| is for k=π (a solution that oscillates with every grid point).

Applying the conditional stability result from last lecture for forward-difference timestepping, we find Δt<Δx2/2. This is necessary and sufficient for stability in the ∞-space case, because any (polynomially bounded) initial condition can be written as a sum of these eikm functions. In fact, this is a kind of reverse Fourier series, a [discrete-"time" Fourier transform](http://en.wikipedia.org/wiki/Discrete-time_Fourier_transform) (DTFT, although here it is space not time). Reviewed how Fourier series can be written in terms of complex exponentials rather than sines and cosines. Noted that k and k+2π give equivalent solutions: this is called [aliasing](http://en.wikipedia.org/wiki/Aliasing), and means that we only need to consider k in \[-π,π\].

When we have boundaries, inhomogeneities, etcetera, then it is usually too hard to compute the eigenvalues exactly; in this case Von Neumann analysis usually gives us at best a necessary condition for stability, but not a sufficient condition. In practice, though, it works very well, although usually we err on the conservative side and make Δt a little bit smaller than the Von Neumann bound might strictly request.

Similarly, analyzed the 2d heat equation with center differences in space and forward differences in time, and showed that the maximum Δt is decreased by a factor of 2. In general, it is decreased proportional to the number of dimensions.

The important consequence of this is: when you refine the discretization in space (decreasing Δx), you must _also refine the discretization in time_ (decreasing Δt) in an explicit scheme (like forward differences in time).

New subject: **Wave equations.** Although we originally wrote the wave equation as a second derivative in time, in order to think about time evolution (either numerically or analytically) it is nicer to write it as a first-order derivative into time. The most obvious way to do this is to introduce a new variable v=∂u/∂t, but this turns out to lead to a somewhat unsymmetrical system of equations that is hard to analyze.

Instead, we will look at the scalar wave equation ∇2u=∂2u/∂t2 in a new way. We will introduct a new vector-valued unknown **v**, defined by ∂**v**/∂t=∇u and ∂u/∂t=∇⋅**v**; showed that this is equivalent to ∇2u=∂2u/∂t2. This leads to a new equation of the form ∂**w**/∂t=Â**w**, where **w**\=(u;**v**) and Â is the 2×2operator Â=(0,∇⋅;∇,0)=(0,div;grad;0). Now the problem looks superficially like the heat/diffusion equation: it is first-order in time, with some new operator Â. But this Â is very different from the old Â=∇2! In fact, we will see that this Â gives Â\*\=-Â: it is **anti-Hermitian**, and from this stems many of the important properties of wave equations.

**Further reading** See for example chapter 6 of the Strang book for a typical textbook treatment of Von Neumann analysis. Something I dislike about this and many textbooks is that it does Von Neumann analysis right away. I prefer considering the dependence of Δt on the eigenvalues of A in general first (where things are both simpler and more general than diving into a specific discretization), and only then finding the eigenvalues of A in the Von Neumann case where we can solve for them exactly.

### Lecture 21: October 25

**Handout:** [notes on the algebraic structure of wave equations](http://math.mit.edu/~stevenj/18.369/wave-equations.pdf)

Âefined the inner product ⟨**w**,**w**'⟩ in the obvious way as ∫(uu'+**v**⋅**v**') and showed that Â\*\=−Â, i.e. that it is anti-Hermitian, for either Dirichlet or Neumann boundary conditions. From this, reprised the proof of real eigenvalues for the Hermitian case to show that now the eigenvalues are purely imaginary. Alternatively, showed that iÂ is a Hermitian operator, so all of the properties of Hermitian operators carry over to Â except that the eigenvalues are multiplied by i.

Showed that an anti-Hermitian operator Â is simply _i_ times a self-adjoint/Hermitian operator −_i_Â, and therefore it inherits the nice properties of self-adjoint operators with one difference: the eigenvalues are _imaginary_ instead of real. If we call the eigenvalues −_i_ω for a real ω, then it is clear that we obtain oscillating solutions with time dependence _e−iωt_. Furthermore, showed that ⟨w,w⟩ is conserved in time, which in the next lecture we will interpret as conservation of "energy".

With the wave equation in a new form ∂**w**/∂t=D**w**, we had derived important properties of D: anti-Hermitian, imaginary eigenvalues, unitary time evolution, conservation of energy, similar to the handout (but at a somewhat simpler level, as the handout is from a graduate class).

As in the notes from the previous lecture, considered the general case of the scalar wave equation with non-constant coefficients a,b>0: b∇⋅(a∇u)=∂2u/∂t2, splitting this up as ∂**v**/∂t=a∇u and ∂u/∂t=b∇⋅**v**. As in the notes, showed that the resulting D operator is still anti-Hermitian but under a modified inner product ⟨**w**,**w**'⟩ = ∫(uu'/b+**v**⋅**v**'/a).

Gave simple example of compression waves in a 1d system that is the limit of springs and masses, from lecture 5.5. If h is the displacement, it is convenient (to obtain conservation of energy in a familar form) to write u=∂h/∂t and v=∂h/∂x. We then get b=1/ρ, a=κ (a "spring constant"), and find that ∫(ρu2+κv2) is conserved. We interpreted this as kinetic+potential energy.

Example: pressure waves in a fluid or gas. In this case, u=P (pressure), **v** is a velocity, a=1/ρ (ρ is density), and b=K (bulk modulus: dP=-KdV/V, relating change in pressure dP to fractional change in volume dV/V). Again this gives a wave equation, with a conserved kinetic+potential energy ∫(ρ|**v**|2+P2/K).

### Lecture 22: October 28

Continued notes from last lecture: Considered the case of Maxwell's equation in vacuum (which you already proved is anti-Hermitian in homework) and gave the corresponding energy density in the EM fields. Mentioned the case of Schrodinger's equation in quantum mechanics, where we only have one time derivative, and the conserved norm is interpreted as conservation of probability. The other cases in the notes are more complicated.

**Traveling waves: D'Alembert's solution**. Considered the 1d scalar wave equation c2∂2u/∂x2\=∂2u/∂t2 on an infinite domain with a constant coefficient c. Showed that any f(x) gives possible solutions u(x,t)=f(x±ct). This is called D'Alembert's solution, and describes the function f(x) "moving" to the left or right with speed c. That is, wave equations have travelling solutions, and the constant c can be interpreted as the speed of these solutions. Adding a hard wall (Dirichlet boundary) is equivalent to looking for an odd solution f(x±ct)−f(−x±ct), which gives an _inverted reflection_ off the wall. (Neumann boundary conditions correspond to even solutions and give non-inverted reflections.) If we have two Dirichlet boundaries, as in a finite stretched string, then we obtain an infinite sequence of inverted reflections which we can write as an infinite series.

Given these solutions, it is attractive to try to write any solution u(x,t) as a superposition of D'Alembert solutions. We can do this if we pick a convenient basis of f(x) functions, and the most convenient basis will turn out to be f(x)=eikx for real k: this leads to [Fourier transforms](http://en.wikipedia.org/wiki/Fourier_transform), which we will return to later. In particlar, we then obtain **planewave** solutions ei(kx±ωt) where ω=±ck (the _dispersion relation_). 2π/k is a spatial wavelength λ, and ω/2π is a frequency f, and from this we find that λf=c, a relation you may have seen before.

There is something suspiciously unphysical about D'Alembert solutions: they travel _without changing shape_, even if f(x) is a very non-smooth shape like a triangle wave. Real waves on strings, etcetera, don't seem to do this. The problem is that real wave equations incorporate a complication that we have not yet considered: the speed c, in reality _depends on ω_, an effect called [dispersion](http://en.wikipedia.org/wiki/Dispersion_%28disambiguation%29), so that different frequency components travel at different speeds and the solution will distort as it travels. Physically, it turns out that this comes down to the fact that materials do not respond instaneously to stimuli, which is mathematically expressed by the fact that the Fourier transformation of the frequency-domain equation ∂2u/∂x2\=-ω2c(&omega)−2u Fourier transform to ∂2u/∂x2\=∂2u/∂t2\*(some function of time) where "\*" is a [convolution](http://en.wikipedia.org/wiki/Convolution) operation. We will come back to this later.

### Lecture 23: October 30

**Handouts:** [pset 7 solutions](pset7sol-f13.pdf). [Notes on Fourier transforms, wave velocity, and dispersion](fourier-dispersion.pdf)

Discretization of the (1d scalar) wave equation: staggered grids and leap-frog schemes. Von Neumann and CFL analysis. Dispersion relation.

Discretization of the (1d scalar) wave equation, simplifying for now to an infinite domain (no boundaries) and constant coefficients (c=1). This corresponds to the equations ∂u/∂t=∂v/∂x and ∂v/∂t=∂u/∂x.

The obvious strategy is to make everything a center difference. First concentrating on the spatial discretization, showed that this means that u and v should be discretized on different grids: for integers m, we should discretize u(mΔx)≈um and v(\[m+0.5\]Δx)≈vm+0.5. That is, the u and v spatial grids are offset, or **staggered**, by Δx/2.

For discretizing in time, one strategy is to discretize u and v at the same timesteps nΔt. Center-differencing then leads to a Crank-Nicolson scheme, which can easily show to be unconditionally stable (albeit implicit) for anti-Hermitian spatial discretizations.

Alternatively, we can ues an explicit **leap-frog** scheme in which u is discretized at times nΔt and v is discretized at times \[n-0.5\]Δt. Sketched out the corresponding staggered grids, difference equations, and leap-frog process.

Went through Von Neumann stability analysis of this leap-frog scheme, and derived the **dispersion relation** ω(k) for **planewave** solutions eikΔx m - iωΔt n. Compared to dispersion relation ω(k)=±c|k| of the analytical equation: matches for small k, but a large mismatch as k approaches π/Δx.

Went over handout, first two pages.

**Further reading:** Strang book, section 6.4 on the leapfrog scheme for the wave equation.

### Midterm: November 1

Midterm [exam](midterm-f13.pdf) and [solutions.](midtermsol-f13.pdf)

### Lecture 24: November 4

Handout from last lecture, through group-velocity derivation.