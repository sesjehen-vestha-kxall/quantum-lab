# -----------------------------------------------------------------------------
# Copyright (c) 2020, Miguel Angel Avila Torres. All Rights Reserved.
# -----------------------------------------------------------------------------
import model.hydrogen as hy
import model.numerical_methods as methods

from tkinter import *

import numpy as np

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure


class QuadratureProbabilityGraph(object):

    def __init__(self, window, arguments: dict):
        """
        :param window: The frame which will draw the plot
        :param arguments: a dictionary with {str: StringVar} entries
        """
        self.window = window
        self.arguments = arguments
        self.button = Button(master=window, text='Simulate', command=self.plot_quadrature)
        self.button.pack()
        self.probability = Label(master=window)
        self.probability.pack()
        self.fig = Figure(figsize=(5, 5))
        self.points = self.fig.add_subplot(projection='polar')
        self.canvas = FigureCanvasTkAgg(self.fig, master=window)
        self.plot_quadrature()

    def plot_quadrature(self):
        r_a, r_n, n = [float(e.get()) for e in self.arguments.values()]
        integral_prob = methods.quadrature_rule(f=lambda r: np.abs(hy.hydrogen_1s(r)) ** 2,
                                                f_second=hy.hydrogen_1s_second, a=r_a, b=r_n, n=int(n))

        self.points.clear()
        self.points.set_title('\nQuadrature Rule of |\u03A8|^2 over r', fontsize=18)

        rad = np.linspace(start=0, stop=r_n, num=int(n))
        azm = np.linspace(start=0, stop=2 * np.pi, num=int(n))
        r, th = np.meshgrid(rad, azm)
        z = np.abs(hy.hydrogen_1s(r)) ** 2

        self.points.pcolormesh(th, r, z, cmap='inferno', shading='auto')
        # self.fig.savefig(fname='quadrature-prob.png')  # it saves locally the plotted image
        self.probability['text'] = '\nProbability: ' + str(integral_prob) + '\n'
        self.canvas.get_tk_widget().pack()
        self.canvas.draw()
