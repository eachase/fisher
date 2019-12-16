
def perturb_mass1(event, epsilon):
    """
    Derivative w.r.t. primary component mass
    """

    return 1



# Functions for numerical derivative computation
deriv_functions = {
    'mass1': perturb_mass1,
    'mass2': None,
    'distance': None,
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

    # Select function for differencing step
    deriv_func = deriv_functions[parameter]

    # Compute waveform for event one step forward
    event_forward = deriv_func(event, epsilon)

    # Compute waveform for event one step backward
    event_backward = deriv_func(event, -epsilon)

    print(event_forward)

    return 0



