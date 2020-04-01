// -*- C++ -*-
//
// michael a.g. aïvázis <michael.aivazis@para-sim.com>
//
// (c) 2013-2020 parasim inc
// all rights reserved
//

#if !defined(altar_extensions_models_reverso_source_h)
#define altar_extensions_models_reverso_source_h


// place everything in my private namespace
namespace altar::extensions::models::reverso {
    // make a new source
    extern const char * const newSource__name__;
    extern const char * const newSource__doc__;
    PyObject * newSource(PyObject *, PyObject *);

    // attach the observations
    extern const char * const data__name__;
    extern const char * const data__doc__;
    PyObject * data(PyObject *, PyObject *);

    // attach the locations of the observation points
    extern const char * const locations__name__;
    extern const char * const locations__doc__;
    PyObject * locations(PyObject *, PyObject *);

    // the structure of the parameter sets
    extern const char * const layout__name__;
    extern const char * const layout__doc__;
    PyObject * layout(PyObject *, PyObject *);

    // compute the predicted displacements that correspond to a set of samples
    extern const char * const displacements__name__;
    extern const char * const displacements__doc__;
    PyObject * displacements(PyObject *, PyObject *);

    // compute the residuals
    extern const char * const residuals__name__;
    extern const char * const residuals__doc__;
    PyObject * residuals(PyObject *, PyObject *);
}

#endif

// end of file
