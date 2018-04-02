#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# michael a.g. aïvázis <michael.aivazis@para-sim.com>
#
# (c) 2013-2018 parasim inc
# (c) 2010-2018 california institute of technology
# all rights reserved
#


def test():
    # get the class
    from altar.mcmc import annealer
    # make an instance
    catmip = annealer()
    # and publish it
    return catmip


# bootstrap
if __name__ == "__main__":
    # run the driver
    test()
    # report success
    raise SystemExit(0)


# end of file
