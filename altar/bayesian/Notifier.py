# -*- python -*-
# -*- coding: utf-8 -*-
#
# michael a.g. aïvázis <michael.aivazis@para-sim.com>
#
# (c) 2013-2018 parasim inc
# (c) 2010-2018 california institute of technology
# all rights reserved
#


# the package
import altar


#
class Notifier(altar.component,
               family="altar.simulations.dispatchers.notifier",
               implements=altar.simulations.dispatcher):

    """
    A dispatcher of events generated during the annealing process
    """


    # constants
    start = "start"
    finish = "finish"


    # protocol obligations
    @altar.export
    def initialize(self, application):
        """
        Initialize me given an {application} context
        """
        # nothing to do
        return self


    # interface
    def register(self, monitor):
        """
        Enable {monitor} as an observer of simulation events
        """
        # go through all known event types
        for event in self.events.keys():
            # check
            try:
                # whether the {monitor} implements a handler for this particular event
                handler = getattr(monitor, event)
            # if not
            except AttributeError:
                # no worries
                continue
            # otherwise, add the handler to the correct event table entry
            self.events[event].addObserver(handler)
        # all done
        return


    def notify(self, event, controller):
        """
        Notify all handlers that are waiting for {event}
        """
        # find the observable associated with this event
        observable = self.events[event]
        # and ask it to notify its observers
        observable.notifyObservers(controller=controller)
        # all done
        return


    # meta-methods
    def __init__(self, **kwds):
        # chain up
        super().__init__(**kwds)

        # establish the table of handled events
        self.events = {
            self.start: altar.patterns.observable(),
            self.finish: altar.patterns.observable(),
        }

        # all done
        return


# end of file
