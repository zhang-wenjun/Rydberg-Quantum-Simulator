# -*- coding: utf-8 -*-
"""
Created on Sun Jan 10 16:45:52 2021

@author: Wenjun Zhang
"""

# %% Import modules

import matplotlib.pyplot as plt
import numpy as np
from tenpy.models import lattice
from tenpy.models.model import CouplingMPOModel
from tenpy.tools.params import asConfig

from tenpy.networks.site import SpinSite
from tenpy.networks.mps import MPS

from tenpy.algorithms import dmrg
from tenpy.algorithms.exact_diag import ExactDiag

# %% Plot lattices

plt.figure(figsize=(5, 4))
ax = plt.gca()
lat = lattice.Honeycomb(4, 4, None, bc='periodic')
lat.plot_coupling(ax, linewidth=3.)
lat.plot_order(ax, linestyle=':')
lat.plot_sites(ax)
lat.plot_basis(ax, origin=-0.25*(lat.basis[0] + lat.basis[1]))
ax.set_aspect('equal')
ax.set_xlim(-1)
ax.set_ylim(-1)
plt.show()

# %% Define a spin-s Heisenberg chain

class SpinSHeisenberg(CouplingMPOModel):
    def init_sites(self, model_params):
        conserve = model_params.get('conserve', None)
        spin = model_params.get('spin', 0.5)
        site = SpinSite(S=spin, conserve=conserve)
        return site
        
    def init_lattice(self, model_params):
        site = self.init_sites(model_params)
        L = model_params.get('L', 2)
        bc_MPS = model_params.get('bc_MPS', 'infinite')
        bc = 'periodic' if bc_MPS == 'infinite' else 'open'
        bc = model_params.get('bc', bc)
        assert bc in ['open', 'periodic']
        lat = lattice.Chain(L, site)
        return lat
        
    def init_terms(self, model_params):
        J = np.asarray(model_params.get('J', 1.))
        for u1, u2, dx in self.lat.pairs['nearest_neighbors']:
            self.add_coupling(-0.5*J, u1, 'Sp', u2, 'Sp', dx)
            self.add_coupling(-0.5*J, u1, 'Sm', u2, 'Sm', dx)
            self.add_coupling(-J, u1, 'Sz', u2, 'Sz', dx)
            
    def __init__(self, model_params):
        model_params = asConfig(model_params, self.__class__.__name__)
        CouplingMPOModel.__init__(self, model_params)
        
# %% Initialize the model

model_params = {
    'J': 1. ,
    'L': 10 ,
    'bc_MPS': 'finite',
    'spin': 0.5,
}

M = SpinSHeisenberg(model_params)

# %% Exact diagonalization

ED = ExactDiag(M)
ED.build_full_H_from_mpo()
print("Hamiltonian = \n", ED.full_H.to_ndarray())
ED.full_diagonalization()  # the expensive part for large L
E0_ED, psi_ED = ED.groundstate()  # return the ground state
print("psi_ED = \n", psi_ED.to_ndarray())
print("ground state energy = ", E0_ED)

# %% DMRG

dmrg_params = {
    'mixer': False, 
    'max_E_err': 1.e-10,
    'trunc_params': {
        'chi_max': 100,
        'svd_min': 1.e-10,
    },
    'verbose': True,
    'combine': True
}

psi = MPS.from_lat_product_state(M.lat, [['up']])
eng = dmrg.TwoSiteDMRGEngine(psi, M, dmrg_params)
E, psi = eng.run()

# %% Compare and print

print("ground state energy = ", E)
Z = psi.expectation_value("Sz")
X = psi.expectation_value("Sx")
Y = psi.expectation_value("Sy")
# psi_ED_MPS = MPS.from_lat_product_state(M.lat, psi_ED)
# Z_ED = psi_ED_MPS.expectation_value("Sz")
x = np.arange(psi.L)
plt.figure()
plt.plot(x, Z, label="Z")
plt.plot(x, Y, label="Y")
plt.plot(x, X, label="X")
# plt.plot(x, Z_ED, label="Z_ED")
plt.xlabel("site")
plt.ylabel("onsite expectation value")
plt.legend()
plt.show()