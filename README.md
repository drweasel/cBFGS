# cBFGS

A small, self-contained implementation of the *Broyden–Fletcher–Goldfarb–Shanno
algorithm* (BFGS) with line-search for unconstraint nonlinear optimisation.

The implementation avoids any memory allocation and is well-suited for
embedded applications.

### Warning

> *This code is experimental. Use it on your own risk and don't blame me if
> something goes wrong. Don't use it in critical applications - there are more
> battle-tested implementations of the BFGS algorithm around.*

If this doesn't scare you, give it a try. I'd be happy to get feedback.

## Literature

See [Wikipedia](https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm)
for a brief introduction. An in-depth description can be found in the following
textbook:

> Jorge Nocedal and Stephen J. Wright: *Numerical Optimization*, Second
> Edition, 2006, Springer Verlag

## License

Copyright (c) 2022 Michael Weitzel <mich@elweitzel.de>

This file is licensed under the terms of the standard MIT License; see the
file LICENSE for details.

