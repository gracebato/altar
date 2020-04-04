// -*- C++ -*-
//
// michael a.g. aïvázis <michael.aivazis@para-sim.com>
// (c) 2013-2020 parasim inc
// all rights reserved
//

// configuration
#include <portinfo>
// external
#include <cmath>
#include <pyre/journal.h>
#include <gsl/gsl_matrix.h>

// support
#include "reverso.h"

// my declarations
#include "Source.h"

// meta-methods
// destructor
altar::models::reverso::Source::
~Source() {
    // make a channel
    pyre::journal::debug_t channel("reverso.source");

    // sign on
    channel
        << pyre::journal::at(__HERE__)
        << "deleting source " << this
        << pyre::journal::endl;

    // if we were handed a matrix of the locations of the observation points
    if (_locations) {
        // release it
        gsl_matrix_free(_locations);
        // tell me
        channel
            << pyre::journal::at(__HERE__)
            << "  released {_locations} matrix at" << _locations
            << pyre::journal::endl;
    }

    // all done
    channel
        << pyre::journal::at(__HERE__)
        << "done deleting source " << this
        << pyre::journal::endl;
}


// interface
void
altar::models::reverso::Source::
displacements(gsl_matrix_view * samples, gsl_matrix * predicted) const {
    // pull the number of samples and parameters from the shape of the {sample} matrix
    auto nSamples = samples->matrix.size1;
    auto nParameters = samples->matrix.size2;

    // clean up the resulting matrix
    gsl_matrix_set_zero(predicted);

    // go through all the samples
    for (auto sample=0; sample<nSamples; ++sample) {
        // unpack the parameters
        // the flow rate
        auto Qin = gsl_matrix_get(&samples->matrix, sample, _QinIdx);
        // the depth of the deep reservoir
        auto H_d = gsl_matrix_get(&samples->matrix, sample, _HdIdx);
        // the ratio of the radii
        auto k = gsl_matrix_get(&samples->matrix, sample, _kIdx);
        // the hydraulic pipe radius
        auto a_c = gsl_matrix_get(&samples->matrix, sample, _acIdx);
        // the magma density difference
        auto drho = gsl_matrix_get(&samples->matrix, sample, _drhoIdx);

        // compute the displacements
        reverso(sample, _locations,
                drho, H_d, k,
                Qin, a_c,
                _G, _v, _mu, _Hs, _as, _g,
                predicted);     
  
    }

    // all done
    return;
}


void
altar::models::reverso::Source::
residuals(gsl_matrix * predicted) const {
    // make a channel
    pyre::journal::debug_t channel("reverso.source");

    // sign on
    channel
        << pyre::journal::at(__HERE__)
        << "computing residuals in place"
        << pyre::journal::endl;

    // unpack the number of samples and number of observations
    auto nSamples = predicted->size1;
    auto nObservations = predicted->size2;

    // go through all observations
    for (auto obs=0; obs < nObservations; ++obs) {
        // get the corresponding measurement
        auto actual = gsl_vector_get(_data, obs);
        // go though the samples
        for (auto sample=0; sample < nSamples; ++sample) {
            // get the predicted displacement
            auto pred = gsl_matrix_get(predicted, sample, obs);
            // compute the difference
            auto residual = pred - actual;
            // and store back into the matrix we were handed
            gsl_matrix_set(predicted, sample, obs, residual);
        }
    }

    // all done
    return;
}


// end of file
