from copy import copy
import numpy as np

def perturb_mass1(event, epsilon):
    """
    Derivative w.r.t. primary component mass
    """

    assert abs(epsilon) < abs(event.mass1)
    assert event.mass1 - epsilon > 0 

    # Create new event
    copied_event = copy(event)
    copied_event.mass1 = event.mass1 + epsilon

    return copied_event

def perturb_mass2(event, epsilon):
    """
    Derivative w.r.t. primary component mass
    """

    assert abs(epsilon) < abs(event.mass1)
    assert event.mass2 - epsilon > 0 

    # Create new event
    copied_event = copy(event)
    copied_event.mass2 = event.mass2 + epsilon

    return copied_event


# Functions for numerical derivative computation
deriv_functions = {
    'mass1': perturb_mass1,
    'mass2': perturb_mass2,
    }



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

    # Select function for differencing step
    deriv_func = deriv_functions[parameter]

    # Compute waveform for event one step forward
    event_forward = deriv_func(event, epsilon)

    # Compute waveform for event one step backward
    event_backward = deriv_func(event, -epsilon)

    # Report the derivative
    assert (event_backward.freqs == event_forward.freqs).all()
    deriv = (event_forward.hp.data.data - event_backward.hp.data.data) / (2 * epsilon)

    return deriv

