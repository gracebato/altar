# -*- python -*-
# -*- coding: utf-8 -*-
#
# michael a.g. aïvázis <michael.aivazis@para-sim.com>
#
# (c) 2013-2020 parasim inc
# (c) 2010-2020 california institute of technology
# all rights reserved
#

# expansion placeholders
date = "2020-04-04T21:10:47Z"
major = 2
minor = 0
micro = 1
revision = "168be39"


# assemble the version tuple
version = (major, minor, micro, revision)


# decorations for the help system
copyright = """
    (c) 2013- ParaSim Inc.
    (c) 2010- California Institute of Technology
    ALL RIGHTS RESERVED
"""

banner = f"""
    reverso {major}.{minor}.{micro} revision {revision}
""" + copyright

header = banner + """
authors:
    Michael Aïvázis     <michael.aivazis@para-sim.com>
    Grace Bato          <mary.grace.p.bato@jpl.nasa.gov>
    Eric Gurrola        <eric.m.gurrola@jpl.nasa.gov>
"""

license = header + """
license:
    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    * Neither the name altar nor the names of its contributors may be used
      to endorse or promote products derived from this software without
      specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

acknowledgments = header + """
acknowledgments:
    This work has been funded in part by NASA and the National Science
    Foundation.

    Sarah Minson wrote the original C version of CATMIP in 2010 as part of
    her Ph.D. thesis.

    In late 2011, Michael Aivazis rewrote CATMIP in python and C++, and
    refactored it into a collection of pyre components.

    altar was born in late 2012 as an attempt to reorganize the previous
    efforts into an extensible framework for exploring earthquake models
    using a broad spectrum of stochastic techniques.
"""


# end of file
