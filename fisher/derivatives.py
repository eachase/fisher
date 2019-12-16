from copy import copy
import numpy as np

def perturb_param(event, param, epsilon):
    """
    Derivative w.r.t. some parameter
    """
    
    # Original value
    param0 = getattr(event, param)

    assert abs(epsilon) < param0
    # FIXME: add some parameter-specific assertions
    #assert event.mass1 - epsilon > 0 

    # Create new event
    copied_event = copy(event)
    setattr(copied_event, param, param0 + epsilon)

    # Compute waveform for new event
    copied_event.waveform(flow=event.flow, deltaf=event.deltaf, 
        fhigh=event.fhigh)


    return copied_event



def derivative(event, parameter, epsilon):
    """
    For a given event, compute the derivative of the waveform 
    with respect to some GW parameter.

    Parameters:
    -----------
    event: gw_event_gen.eventgen.event.Event
        Reference GW event for derivatives

    parameter: string
        Parameter to compute derivatives for

    epsilon: float
        Difference uses in finite differencing

    """

    assert epsilon > 0

    # Check that waveform exists
    assert np.all(np.asarray([hasattr(event, attr) for attr \
        in event._waveattrs+('hx','hp')]))

    # Compute waveform for event one step forward
    event_forward = perturb_param(event, parameter, epsilon)

    # Compute waveform for event one step backward
    event_backward = perturb_param(event, parameter, -epsilon)

    # Report the derivative
    assert (event_backward.freqs == event_forward.freqs).all()
    deriv = (event_forward.hp.data.data - event_backward.hp.data.data) / (2 * epsilon)

    return deriv

