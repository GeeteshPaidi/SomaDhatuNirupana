import numpy as np
from scipy import optimize
from element_handler import ElementHandler


class XRFFitter:
    def __init__(self, counts, std_dev=0.1):
        self.counts = np.array(counts)
        self.num_channels = len(self.counts)
        self.energy_factor = 0.0135 if self.num_channels == 2048 else 0.0277
        self.energies = np.arange(self.num_channels) * self.energy_factor

        self.handler = ElementHandler(std_dev=std_dev, num_channels=self.num_channels)
        self.elements = self.handler.elements
        self._setup_params()

    def _setup_params(self):
        self.bounds = [(0, 200) for _ in self.elements] + [(1e-8, 1), (0.01, 0.5)]
        self.initial_guess = [self.handler.conc[el] for el in self.elements] + [self.handler.scale, 0.1]

    def calculate_model(self, params):
        concentrations, scale, std_dev = params[:len(self.elements)], params[-2], params[-1]

        for i, el in enumerate(self.elements):
            self.handler.conc[el] = concentrations[i]
            self.handler.std_dev[el] = std_dev
            self.handler.element_models[el].conc = concentrations[i]
            self.handler.element_models[el].std_dev = std_dev
            for lt in self.handler.element_models[el].std_dev_dict:
                self.handler.element_models[el].std_dev_dict[lt] = std_dev

        self.handler.scale = scale
        return self.handler.calculate_folded_intensity(self.energies) * scale

    def residuals(self, params):
        return self.counts - self.calculate_model(params)

    def fit(self, method='leastsq'):
        if method == 'leastsq':
            result = optimize.least_squares(self.residuals, self.initial_guess, bounds=tuple(zip(*self.bounds)), method='trf')
        else:
            result = optimize.minimize(lambda p: np.sum(self.residuals(p)**2), self.initial_guess, method=method, bounds=self.bounds)
        self.result = result
        return result.x

