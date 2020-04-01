// -*- C++ -*-
//
// michael a.g. aïvázis <michael.aivazis@para-sim.com>
//
// (c) 2013-2020 parasim inc
// all rights reserved
//

// for the build system
#include <portinfo>
// external dependencies
#include <string>
#include <Python.h>

// the module method declarations
#include "exceptions.h"
#include "source.h"


// put everything in my private namespace
namespace altar::extensions::models::reverso {
    // the module method table
    extern PyMethodDef module_methods[];
    extern PyModuleDef module_definition;
}

PyMethodDef
altar::extensions::models::reverso::
module_methods[] = {
    // source methods
    // constructor
    { newSource__name__, newSource, METH_VARARGS, newSource__doc__ },
    // user supplied information
    { data__name__, data, METH_VARARGS, data__doc__ },
    { locations__name__, locations, METH_VARARGS, locations__doc__ },
    { layout__name__, layout, METH_VARARGS, layout__doc__ },
    // the calculation of the displacements
    { displacements__name__, displacements, METH_VARARGS, displacements__doc__ },
    // and the residuals
    { residuals__name__, residuals, METH_VARARGS, residuals__doc__ },

    // sentinel
    {0, 0, 0, 0}
};

// the module definition structure
PyModuleDef
altar::extensions::models::reverso::
module_definition = {
    // header
    PyModuleDef_HEAD_INIT,
    // the name of the module
    "reverso",
    // the module documentation string
    "the reverso extension module",
    // size of the per-interpreter state of the module; -1 if this state is global
    -1,
    // the methods defined in this module
    module_methods
};

// initialization function for the module
// *must* be called PyInit_altar
PyMODINIT_FUNC
PyInit_reverso()
{
    // create the module
    PyObject * module = PyModule_Create(&altar::extensions::models::reverso::module_definition);
    // check whether module creation succeeded
    if (!module) {
        // and raise an exception if not
        return 0;
    }
    // otherwise, we have an initialized module
    // return the newly created module
    return module;
}

// end of file
