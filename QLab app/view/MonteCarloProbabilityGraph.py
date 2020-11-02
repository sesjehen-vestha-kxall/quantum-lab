# -----------------------------------------------------------------------------
# Copyright (c) 2020, Miguel Angel Avila Torres. All Rights Reserved.
# -----------------------------------------------------------------------------
import model.hydrogen as hy
import model.numerical_methods as methods

import numpy as np

from tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib import pyplot as plt


class MonteCarloProbabilityGraph(object):

    def __init__(self, window: Frame, arguments: dict, meridians=50):
        """
        :param window: The frame which will draw the plot
        :param arguments: a dictionary with {str: StringVar} entries
        :param meridians: the number of meridians enclosing the radius r_b in the integral
        """
        self.window = window
        self.arguments = arguments
        self.meridians = meridians
        self.button = Button(master=window, text='Simulate', command=self.__plot_monte_carlo_probability)
        self.button.pack()
        self.probability = Label(master=window)
        self.probability.pack()
        self.error = Label(master=window)
        self.error.pack()
        self.fig = Figure(figsize=(5, 5))
        self.points = self.fig.add_subplot(111, projection='3d')
        self.canvas = FigureCanvasTkAgg(self.fig, master=window)
        self.__plot_monte_carlo_probability()

    @staticmethod
    def wire_frame_sphere(radius=1E-10, n_meridians=20, n_circles_latitude=None):
        """
        :param radius: the radius for the sphere limiting with the meridians
        :param n_meridians: the number of meridians
        :param n_circles_latitude: the latitude of the meridians
        :return: the coordinates for the meridians
        """
        if n_circles_latitude is None:
            n_circles_latitude = max(n_meridians / 2, 4)
        u, v = np.mgrid[0:2 * np.pi:n_meridians * 1j, 0:np.pi:n_circles_latitude * 1j]
        sphere_x = radius * np.cos(u) * np.sin(v)
        sphere_y = radius * np.sin(u) * np.sin(v)
        sphere_z = radius * np.cos(v)
        return sphere_x, sphere_y, sphere_z

    def __plot_monte_carlo_probability(self):
        r_a, r_b, theta_a, theta_b, phi_a, phi_b, samples = [float(e.get()) for e in self.arguments.values()]

        def psi_100_(r_i, theta_i, phi_i):
            return np.abs(hy.hydrogen_1s(r_i, theta_i, phi_i)) ** 2

        p_1s, r, theta, phi, volume = methods.spherical_monte_carlo_integral(
            psi_100_, r_a, r_b, theta_a, theta_b, phi_a, phi_b, int(samples)
        )

        self.points.clear()
        self.points.set_title('\nMonte Carlo Integral of |\u03A8|^2\nin three dimensions', fontsize=18)

        x_ = r * np.sin(theta) * np.cos(phi)
        y_ = r * np.sin(theta) * np.sin(phi)
        z_ = r * np.cos(theta)

        self.points.scatter3D(x_, y_, z_, color='blue', linewidths=0.1, alpha=0.1)
        self.points.plot_wireframe(*self.wire_frame_sphere(r_b, self.meridians), color='orange', alpha=0.2)
        # self.fig.savefig(fname='monte-carlo.png')  # it saves locally the plotted image
        self.probability['text'] = '\nProbability: ' + str(p_1s) + '\n'
        self.error['text'] = 'Error: ' + str(methods.spherical_monte_carlo_estimation_error(psi_100_, r, theta, phi,
                                                                                            volume, samples))
        self.canvas.get_tk_widget().pack()
        self.canvas.draw()
