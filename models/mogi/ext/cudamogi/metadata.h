// -*- C++ -*-
//
// michael a.g. aïvázis <michael.aivazis@para-sim.com>
//
// (c) 2013-2018 parasim inc
// all rights reserved
//

#if !defined(altar_extensions_models_cudamogi_metadata_h)
#define altar_extensions_models_cudamogi_metadata_h


// place everything in my private namespace
namespace altar {
    namespace extensions {
        namespace models {
            namespace cudamogi {
                // version
                extern const char * const version__name__;
                extern const char * const version__doc__;
                PyObject * version(PyObject *, PyObject *);
            } // of namespace cudamogi
        } // of namespace models
    } // of namespace extensions
} // of namespace altar

#endif

// end of file
