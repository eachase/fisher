import numpy as np
import pycbc.psd

from fisher import derivatives

# Epsilons for finite differencing
epsilons = {
    'mass1': 1e-5,
    'mass2': 1e-5,
    'distance': 1e-2,
    }

class FisherMatrix(object):

    def __init__(self, event, params,
        ifo='H1', flow=10, deltaf=0.125, fhigh=2048):


        self.params = np.asarray(params)
        self.numparams = len(params)


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
        fisherparams = {}
        for param in params:
            fisherparams[param] = FisherParameter(event, param)

        # Form Fisher Matrix
        FIM = np.zeros((self.numparams, self.numparams))
        for i, param_i in enumerate(params):
            for j, param_j in enumerate(params):
                deriv_i = fisherparams[param_i].derivative
                deriv_j = fisherparams[param_j].derivative
                FIM[i, j] = 4 * np.real(np.nansum((
                    deriv_i*np.conjugate(deriv_j) / self.psd * deltaf).numpy()))
        self.fishermatrix = FIM    

        # Pre-compute basic Fisher Matrix statistics
        self._covariance()
        self._inverse()


    def __getitem__(self, param):
        """
        Return FIM for this parameter
        """

        param_idx = np.where(self.params == param)[0][0]
        return self.fishermatrix[param_idx][param_idx]



    def _covariance(self):
        self.covariance_matrix = np.cov(self.fishermatrix)
        return self.covariance_matrix


    def _inverse(self):
        if not np.all(self.covariance_matrix == 0):
            self.inverse_matrix = np.linalg.inv(self.fishermatrix)
        else:
            # FIXME: make this a warning
            print('Fisher Matrix not invertible')
            self.inverse_matrix = None

        return self.inverse_matrix


    def stdev(self, param):

        stdev = None

        param_idx = np.where(self.params == param)[0][0]
        try:    
            stdev = np.sqrt(self.inverse_matrix[param_idx][param_idx])
        except:
            raise ValueError('Fisher Matrix not invertible')


        return stdev



# FIXME: make a general matrix class and initiate the covariance, inverse, and FIM as Matrix objects
# (Might be overkill given what numpy already has available)



class FisherParameter(object):

    def __init__(self, event, parameter):
        self.epsilon = epsilons[parameter]
        self.derivative = derivatives.derivative(event, parameter, self.epsilon)




