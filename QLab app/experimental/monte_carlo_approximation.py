# -----------------------------------------------------------------------------
# Copyright (c) 2020, Miguel Angel Avila Torres. All Rights Reserved.
# -----------------------------------------------------------------------------
import model.hydrogen as hy
import model.numerical_methods as methods

import pandas as pd
import numpy as np

data = np.empty(shape=(methods.MONTE_CARLO_SAMPLING, 5), dtype='float32', order='C')

idx: int = 0
r, theta, phi = methods.uniform_spherical_interval(0, 6.7E-11, 0, np.pi * 2, 0, np.pi)
summation: float = 0
for r_i, theta_i, phi_i in zip(r, theta, phi):
    data[idx] = int(idx + 1), r_i, theta_i, phi_i, (np.abs(hy.hydrogen_1s(r_i))**2)
    idx += 1

monte_carlo_approximation = pd.DataFrame(data=data, columns=['sample', 'r', 'theta', 'phi', '|\u03A8(r)|^2'])

monte_carlo_approximation.to_excel(excel_writer='monte-carlo.xlsx', sheet_name='hydrogen 1s',
                                   encoding='utf-8', index=False)
