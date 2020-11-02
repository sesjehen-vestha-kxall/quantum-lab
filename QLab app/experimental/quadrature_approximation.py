# -----------------------------------------------------------------------------
# Copyright (c) 2020, Miguel Angel Avila Torres. All Rights Reserved.
# -----------------------------------------------------------------------------
from model import hydrogen as hy
import model.numerical_methods as methods

import pandas as pd
import numpy as np

# Error Bound
# print(np.sqrt((2.629029177212380E-07 ** 3)/(20*np.pi*(hy.a_0**5))))

quadrature_approximation = pd.DataFrame(columns=['iteration',
                                                 'inferior limit',
                                                 'superior limit',
                                                 'quadrature of \u03A8^2',
                                                 'difference error',
                                                 'relative error'])

sup: float = hy.a_0
inf: int = 0

ans_quad = quad = np.inf
iteration: int = int(1)  # Type Conversion
while quad > 1:
    quad = methods.quadrature_rule(f=lambda r: np.abs(hy.hydrogen_1s(r)) ** 2,
                                   f_second=hy.hydrogen_1s_second,
                                   a=inf,
                                   b=sup,
                                   n=100)

    quadrature_approximation.loc[iteration - 1] = [int(iteration),
                                                   inf,
                                                   sup,
                                                   quad,
                                                   abs(ans_quad - quad),
                                                   abs((1 - quad))]

    ans_quad = quad
    sup += 1E-11
    iteration += 1

quadrature_approximation.to_excel(excel_writer='quadrature.xlsx', sheet_name='hydrogen 1s',
                                  encoding='utf-8', index=False)
