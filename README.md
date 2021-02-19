# THE FUNDAMENTAL THEOREM OF FOURIER SERIES
*Supplment and Problems*

Watch the video! https://youtu.be/7BJgmC0T5xQ

Using computation tools and seeing lots of examples is a very important part of groking fourier series.  This repository contains the source code used to create the youtube video "The Fundamental Theorem of Fourier Series"(https://www.youtube.com/watch?v=7BJgmC0T5xQ), and requires math inspector version 0.9.2+, you can download the latest version from the website https://mathinspector.com, or by running the command

```
$ pip install mathinspector
```

To open this project in math inspector, select File > Open from the main menu and choose the `.math` file.  Once the project is open, you should see a large number of items in the node editor.  These items have all been set up so that you can right click on them and choose "Plot" and it will work.  Choose "View Doc" from the right click menu to learn more about any of the items.

Feel free to open a PR if you find ways to improve the code, or to open an issue if you have any questions or find bugs.

## Why should we care about the reason Fourier Series converge?
Fourier series are kind of a big deal.  They are an extremely popular topic not only among mathematicians and physicists, but also among programmers.  The main reason for this is probably that the animations for their convergence are quite beautiful.  There is something very mysterious when you see this phenomenon for the first time.  It's natural to see this animation, go home and think about it, and then ask the question: is the fourier series just approximating the function, or is the fourier series what the function **is**?  To answer that question, you need to know if, how, when, and why they converge.

## The discontinous function h(x)
`hfunction.py`
The function `h(x)` serves as a great introduction to computational fourier series because its fourier expansion is one of the simplest possible.  

`h(x) = 0.5 * (pi - x)`
`h(0) = 0`
`h(x + 2pi) = h(x)`
It's expansion is given by

```
h(x) = sin(nx) + (1/2)sin(2x) + (1/3)sin(3x) + ...
```

`h(x)` plays a central role in the convergence proof.  The tl;dr is that you use `h(x)` to glue together an arbitrary periodic and piecewise smooth function `f(x)` and then you can use the glued together version to prove that `f(x)` equals it's own Fourier Series.

Each function can be viewed as a point in an infinite dimensional vector space which is called "Hilbert Space", and we can do geometry in this space.  Just like an arbitrary irrational number has a non-repeating decimal expansion in terms of the symbols 0, 1, ..., 9; functions also have their own "decimal expansion".  The complex number in the n-th decimal place of an arbitrary function `f(x)` is it's n-th fourier coefficient.  Instead of 0-9 we are using the continously infinite set of all complex numbers for each decimal place.

There is a kind of algorithmic beauty to the way the proof works on a technical level.  `h(x)` along with it's super simple fourier expansion is plays a fundamental role.  Its sort of like a really well written computer program, where all the parts of the code work together very cleanly but in the way only a really brilliant programmer could think of.

## The Weierstass approximation
`weierstrass.py`
The code in this file is capable of plotting the weierstrass approximation for a function of a single real variable, `f(x)`, and running animations.

Can anyone optimize the weierstrass implementation details and optimize everything for performance during animations?

## 3D visualization of the 2D Weierstrass approximation
`3d_weierstrass.py`
I thought it would be super cool to represent the convergence of the 2-dimensional Weierstrass approximation by representing the value of P_n(x,y), for a fixed number of terms n, as the height of a surface in 3-dimensional space.  This is sort of working(?)

Does anyone know what is going on with the regions on the boundary of the 2-dimensional appriximation and make the animation not have all these nasty artifacts?

## Fourier series of an arbitrary function
`fourier.py`
There are two functions in this file, one of them computes the fourier coefficients of an arbitrary function, and the other one computes the partial sums of the fourier series of an arbitrary function for a fixed number of terms.

It can be interesting to see how the complex coefficients change depending on small perterbations of the test function you are using.  It would be great to have some more widgets in the node editor which can transform functions continuously, and then be able to watch how the fourier coefficients change.  For now there is only one widget that does this, which is the "scale" function (it's in the file `test_functions.py`, see the docstring for more information about how to use it)

At this point, the performance is pretty slow for the partial sums when the number of terms is too large.  This is because the implementation is quick and dirty, and it's not taking advantage of the stuff numpy can do under the hood.  There are also some issues with the implementation where it no longer approximates the function when the number of terms becomes very large, most likely this is due to rounding errors from the coefficients.  It would be nice to figure out what's going on for large values of N here and fix this function up so it works properly.