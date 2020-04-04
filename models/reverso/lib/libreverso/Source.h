// -*- C++ -*-
//
// michael a.g. aïvázis <michael.aivazis@para-sim.com>
// (c) 2013-2020 parasim inc
// all rights reserved
//

// code guard
#if !defined(altar_models_reverso_Source_h)
#define altar_models_reverso_Source_h

// external
#include <vector>

// forward declarations
namespace altar::models::reverso {
    // forwards declarations
    class Source;

    // type aliases
    using source_t = Source;
}


// the reverso source
class altar::models::reverso::Source {
    // types
public:
    using size_type = std::size_t;

    // meta-methods
public:
    virtual ~Source();
    inline Source(double G, double v, double mu, double H_s, double a_s, double g);

    // interface
public:
    inline void data(gsl_vector * data);
    inline void locations(gsl_matrix * locations);
    inline void layout(size_type QinIdx,
                       size_type drhoIdx, size_type HdIdx,
                       size_type kIdx, size_type acIdx);

    void displacements(gsl_matrix_view * samples, gsl_matrix * predicted) const;
    void residuals(gsl_matrix * predicted) const;

    // implementation details
private:
    gsl_vector * _data;        // borrowed reference
    gsl_matrix * _locations;   // owned

    // Fixed parameters
    double _G;
    double _v;
    double _mu;
    double _Hs;
    double _as;
    double _g;

    // the layout of the various parameters within a sample
    size_type _QinIdx;
    size_type _drhoIdx;
    size_type _HdIdx;
    size_type _kIdx;
    size_type _acIdx;
};

// the implementations of the inline methods
#define altar_models_reverso_Source_icc
#include "Source.icc"
#undef altar_models_reverso_Source_icc

// code guard
#endif

// end of file
