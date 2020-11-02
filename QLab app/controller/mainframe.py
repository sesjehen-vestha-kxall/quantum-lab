# -----------------------------------------------------------------------------
# Copyright (c) 2020, Miguel Angel Avila Torres. All Rights Reserved.
# -----------------------------------------------------------------------------
from tkinter import *

import view.MonteCarloProbabilityGraph as MonteCarlo
import view.QuadratureProbabilityGraph as Quadrature
import model.numerical_methods as methods

import numpy as np


app = Tk()
main_frame = Frame()


def fill_data_captors(arguments: dict, default_values: list, plotter: Frame):
    """
    Fills the data captors (Entries) that will be passed to
    the plotter so that this can request input data to update the graphs
    :param arguments: the arguments for the plotter
    :param default_values: the arguments default values
    :param plotter: the Frame which will plot a required numerical method
    """
    y = 0
    for arg in arguments.keys():
        arguments[arg] = StringVar(app)
        Label(master=plotter, text=arg, padx=10, pady=5).grid(row=y, column=0)
        Entry(master=plotter, width=20, textvariable=arguments[arg]).grid(row=y, column=1)
        arguments[arg].set(default_values[y])
        y += 1


"""
Assemble section for Hydrogen's 1s wave function probability
"""
monte_carlo_plotter = Frame(master=app, relief='raised', borderwidth=5)
monte_carlo_plotter.grid(row=0, column=0)

spatial_integral_arguments = {'r_a': None, 'r_b': None, '\u03B8_a': None,
                              '\u03B8_b': None, '\u03D5_a': None, '\u03D5_b': None, 'samples': None}
spatial_integral_values = [0, 6.73E-11, 0, np.pi * 2, 0, np.pi, methods.MONTE_CARLO_SAMPLING]

fill_data_captors(spatial_integral_arguments, spatial_integral_values, monte_carlo_plotter)

monte_carlo_graph = Frame(master=app, width=750, height=470, relief='raised', borderwidth=5)
monte_carlo_graph.grid(row=0, column=1)
MonteCarlo.MonteCarloProbabilityGraph(monte_carlo_graph, spatial_integral_arguments)


"""
Assemble section for quadrature integral of Hydrogen's 1s squared wave function
"""
quadrature_plotter = Frame(master=app, relief='raised', borderwidth=5)
quadrature_plotter.grid(row=0, column=2)

one_dimensional_arguments = {'r_0': None, 'r_N': None, 'N': None}
one_dimensional_values = [0, 2.6289491772072073653583319623183722768544612335972487926483154296875E-7, 100]

fill_data_captors(one_dimensional_arguments, one_dimensional_values, quadrature_plotter)

quadrature_graph = Frame(master=app, width=750, height=470, relief='raised', borderwidth=5)
quadrature_graph.grid(row=0, column=3)
Quadrature.QuadratureProbabilityGraph(quadrature_graph, one_dimensional_arguments)

app.mainloop()
