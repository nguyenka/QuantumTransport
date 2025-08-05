Introduction
============

QuantumTransport is a package for calculating electron spin polarization and 
transmission of a molecular junction system. Calculations of electron
transmission through a molecule (central region) connected to (semi-infinite)
contacts (left and right electrodes) are based on the non-self-consistent
Green’s function method, an approximation to the non-equilibrium Green’s
function (NEGF) method.

Spin polarization in the electron transmission of closed- and open shell
systems (such as radicals or diradicals) can be performed by interfacing with
electronic structure packages. For closed-shell systems, it requires
electronic structure calculations taking spin−orbit coupling (SOC) into
account.  Greeen's functions of the left and right electrodes are approximated
within the wide-band limit, using a constant local density of states, cLDOS.
Coupling terms of electrodes are neglected.

The Hamiltonian for transmission function for non-interacting coherent tunneling of electrons through a left (L) electrode to central (C) molecule and right (R) electrode system is partitioned as:

.. math::

    H_{\text{system}} = \begin{pmatrix}
                    H_L         & H_{LC}      & 0           \\
                    H_{CL}      & H_C         & H_{CR}      \\
                     0          & H_{RC}      & H_{RC}      \\
    \end{pmatrix},

where `H_{L/R}` describes the left/right electrode layer, and `H_C`,
the central region. The `H_{CL}` and `H_{LR}` are the coupling matrices 
between left and right electrodes with the central region.  It is assumed that the electrodes (L and R) are not coupled each other in the Hamiltonian.  The overlap matrix has an analogous partition. The transmission function is given as:

.. math::

    T_{RL} = Tr[\Gamma_LG\Gamma_R G^\dagger].

To evaluate transmission function `T`, the Green’s function `G` and the bulk broadening
function Γ defined as

.. math::

    G = \frac{1} {ES_C-H_C-\Sigma_L-\Sigma_R},

and

.. math::

    \Gamma_X = -2\text{img}(\Sigma_X).

The self-energy of the electrodes X (left or right) is calculated using the overlap (`S_{XC}`) and Hamiltonian (`H_{XC}`) matrix elements between the central region and the electrode as 

.. math::

    \Sigma_X = (ES_{XC} - H_{XC})^\dagger \mathbf{g}_X (ES_{XC} - H_{XC}),

where :math:`\mathbf{g}_X` is the Green’s function of electrode X, which can be approximated within wide-band limit using cLDOS. The overlap matrix `S` obtained from one-component scalar-relativistic calculation, where there is no SO interaction.  The effective single-particle Hamiltonian matrices are approximate by KS matrices


References
----------
  | Jacob, D.; Palacios, J. J.
  | `Critical comparison of electrode models in density functional theory based quantum transport calculations.`__
  | *J. Chem. Phys.* **2011**, 134 (4), 044118.

__ http://dx.doi.org/10.1063/1.3526044

  | Zöllner, M. S.; Varela, S.; Medina, E.; Mujica, V.;Herrmann, C.
  | `Insight into the Origin of Chiral-Induced Spin Selectivity from a Symmetry Analysis of Electronic Transmission`__
  | *J. Chem. Theory Comput.* **2020**, 16, 2914-2929. 

__ https://dx.doi.org/10.1021/acs.jctc.9b01078?ref=pdf


