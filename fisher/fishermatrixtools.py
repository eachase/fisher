import numpy as np
import pycbc.psd

# Epsilons for finite differencing
epsilons = {
    'mass1': 1e-5,
    'mass2': 1e-5,
    }

class FisherMatrix(object):

    def __init__(self, event, params,
        ifo='H1', flow=10, deltaf=0.125, fhigh=2048):

        # Compute waveform, if not already done
        if not np.all(np.asarray([hasattr(event, attr) for attr \
            in event._waveattrs+('hx','hp')])):
            event.waveform(flow=flow, deltaf=deltaf, fhigh=fhigh)

        # Define PSD
        flen = int(fhigh / deltaf) + 1
        self.psd = pycbc.psd.from_string('aLIGODesignSensitivityP1200087',
            flen, deltaf, flow)
        # FIXME: allow for different PSDs

        # Compute necessary derivatives

        # Form Fisher Matrix



class FisherParameter(object):

    def __init__(self, parameter):
        self.epsilon = epsilons[parameter]
        self.derivative = 0




