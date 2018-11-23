# -*- python -*-
# -*- coding: utf-8 -*-
#
# michael a.g. aïvázis <michael.aivazis@para-sim.com>
#
# (c) 2013-2018 parasim inc
# (c) 2010-2018 california institute of technology
# all rights reserved
#


# externals
import altar
# the pure python implementation of the Mogi source
from altar.models.mogi import libmogi


# declaration
class Fast:
    """
    A strategy for computing the data log likelihood that is written in pure python
    """

    # interface
    def initialize(self, model):
        """
        Initialize the strategy with {model} information
        """
        # build the calculator and attach it
        self.source = source = libmogi.newSource(model.nu)
        # get the number of observations
        observations = model.observations
        # the locations on the ground where the observations were made
        locations = model.points
        # the observed displacements
        displacements = model.d
        # the array with the lines of sight to the observation locations
        los = model.los
        # and the data set id for each observation
        oid = model.oid

        # attach the coordinates of the observation points
        libmogi.locations(source, locations)
        # attach the observed displacements
        libmogi.data(source, displacements.data)
        # attach the LOS vectors
        libmogi.los(source, los.data)
        # attach the map of observations to their set
        libmogi.oid(source, oid)
        # inform the source about the parameter layout; assumes contiguous parameter sets
        libmogi.layout(source, model.xIdx, model.dIdx, model.sIdx, model.offsetIdx)

        # nothing to do
        return self


    def dataLikelihood(self, model, step):
        """
        Fill {step.data} with the likelihoods of the samples in {step.theta} given the available
        data.
        """
        # grab my calculator
        source = self.source
        # compute the portion of the sample that belongs to this model
        θ = model.restrict(theta=step.theta)
        # get the norm
        norm = model.norm
        # the inverse of the data covariance matrix
        cd_inv = model.cd_inv
        # the normalization
        normalization = model.normalization
        # and the data likelihood vector
        dataLLK = step.data

        # allocate a matrix to hold the predicted displacements
        predicted = altar.matrix(shape=(step.samples, model.observations))

        # compute the predicted displacements
        libmogi.displacements(source, θ.capsule, predicted.data)
        # compute the residuals (in place)
        libmogi.residuals(source, predicted.data)

        # find out how many samples in the set
        samples = θ.rows
        # go through the samples
        for sample in range(samples):
            # get the residuals
            residuals = predicted.getRow(sample)
            # compute the norm, normalize, and store it
            dataLLK[sample] = normalization - norm.eval(v=residuals, sigma_inv=cd_inv) / 2

        # all done
        return self


    # meta-methods
    def __init__(self, **kwds):
        # chain up
        super().__init__(**kwds)
        # all done
        return


    # private data
    source = None


# end of file
