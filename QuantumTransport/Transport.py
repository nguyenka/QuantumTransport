#!/usr/bin/env python3
# pylint: disable=invalid-name
"""Module for electron spin polarization and transmission."""
import json
from os.path import isfile
import sys
from math import pi
import matplotlib.pyplot as plt
import numpy as np
#from scipy.linalg import sqrtm
#test 2 by Paul
from .cclib_interface import TurbomoleParse, GaussianParse
from .time_util import *

def banner():
    # pylint: multiple-statements
    """
    function for a banner.
    """
    print("----------------------------------------------------------------")
    print("Program to calculate electron spin polarization and transmission")
    print("----------------------------------------------------------------")
    #print("By using this software you accept the terms of the GNU public ")
    #print("license in COPYING, and agree to attribute the use of this software")
    #print("in publications as: \n")
    #print("K. A. Nguyen, QuantumTransport 0.1 (2022)")
    #print("----------------------------------------------------------------")

def plot_spinpol(fname,yscale):
    """
    plot spin polarization and four transmissions functions vs energy 
    """

    mat0 = np.loadtxt(fname)
    E = mat0[:,0]
    Polar = mat0[:,1]
    plt.plot(E,Polar)
    plt.xlabel(r'$E - E_f$ (eV)')
    plt.ylabel('% Polarization')
    plt.show() # show figure
    Tuu = mat0[:,2]
    Tdd = mat0[:,3]
    Tdu = mat0[:,4]
    Tud = mat0[:,5]
    # transmission
    if yscale is not None:
        plt.yscale(yscale) # linear,logit or symlog
    else:
        plt.yscale('log') # default
    legend= []
    legend.append(r"$T^{\uparrow\uparrow}$")
    plt.plot(E,Tuu)
    legend.append(r"$T^{\downarrow\downarrow}$")
    plt.plot(E,Tdd)
    legend.append(r"$T^{\downarrow\uparrow}$")
    plt.plot(E,Tdu)
    legend.append(r"$T^{\uparrow\downarrow}$")
    plt.plot(E,Tud)
    plt.ylim([0.001,10.0])
    plt.xlabel(r'$E - E_f$ (eV)')
    plt.ylabel('Transmission')
    plt.title('Electron Transmission')
    plt.legend(legend)
    plt.show() # show figure

def plot_spinpol2(fname,yscale):
    """
    plot spin polarization and Tuu, Tdd transmission functions vs energy 
    """

    mat0 = np.loadtxt(fname)
    E = mat0[:,0]
    Polar = mat0[:,1]
    plt.plot(E,Polar)
    plt.xlabel(r'$E - E_f$ (eV)')
    plt.ylabel('% Polarization')
    plt.show() # show figure
    Tuu = mat0[:,2]
    Tdd = mat0[:,3]
    # transmission
    legend= []
    legend.append(r"$T^\uparrow$")
    if yscale is not None:
        plt.yscale(yscale) # linear,logit or symlog
    else:
        plt.yscale('log') # default
    plt.plot(E,Tuu)
    legend.append(r"$T^\downarrow$")
    plt.plot(E,Tdd)
    plt.ylim([0.0001,10])
    plt.xlabel(r'$E - E_f$ (eV)')
    plt.ylabel('Transmission')
    plt.title('Electron Transmission')
    plt.legend(legend)
    plt.show() # show figure

def sqrtm_np(A):
    """
    function matrix square root.
    """
    #diagonalization
    evalues, evectors = np.linalg.eig(A)
    # Ensuring square root matrix exists
    assert (evalues >= 0).all()
    result = evectors @ np.diag(np.sqrt(evalues)) @ np.linalg.inv(evectors)
    return result

def dprint_mat(A,name='matrix'):
    """
    print a few elements of a matrix
    """
    print(f'{name}: {A.shape}')
    print(f'{name}: {A[0,0]:12.5e} {A[0,1]:12.5e} {A[1,0]:12.5e} {A[0,2]:12.5e}')

class TransportCal:
    """
    Class for electron spin polarization and transmission calculations.
    """
    # pylint: disable=too-many-instance-attributes

    def __init__(self,qcfile='ridft.out',ciss_calc=False,load_data=True,\
                 fock_calc = False,debug = False):

        """
        Parameters:

        ciss_calc: Chiral induced spin polarization transmission with SOC data.
        load_data: Load files produced by qt_dat.py.
        """
        # pylint: disable=too-many-arguments

        banner()
        self.debug=debug
        self.load_data=load_data
        self.fock_calc=fock_calc
        if load_data:
            with open('basis_data.json') as f:
                bparams = json.load(f)
                self.nL=bparams['left_basis']
                self.nC=bparams['central_basis']
                self.nR=bparams['right_basis']
                self.nmo=bparams['total_basis']
        else:
            if ciss_calc:
                tparse = TurbomoleParse(tout=qcfile,savefile = False)
                self.nL,self.nC,self.nR,self.nmo,self.Sf,self.Fs,self.Fsoc=tparse.parse_output()
            else:
                gparse = GaussianParse(gout=qcfile,savefile = False)
                self.nL,self.nC,self.nR,self.nmo,self.Sf,self.Fa,self.Fb = gparse.parse_output()

        print('Basis functions for the partitioned Hamiltonian/overlap matrix:')
        print('---------------------------------------------------------------')
        print(f'Left electrode: {self.nL}')
        print(f'Central: {self.nC}')
        print(f'Right electrode: {self.nR}')
        print(f'Total: {self.nmo}\n')
        if isfile('transport_params.json') :
            f2 = open('transport_params.json')
            tparams = json.load(f2)
            print('Input parameters for electron transport:')
            print('----------------------------------------')
            self.cLDOS = tparams['cLDOS'] # a constant local density of state
            self.Efermi= tparams['Efermi']
            self.eini  = tparams['eini']
            self.eend  = tparams['eend']
            self.estep = tparams['estep']
        else:
            print('Default input parameters for electron transport:')
            print('------------------------------------------------')
            self.cLDOS = 0.036 # a constant local density of state for Au
            self.Efermi= -5.0 
            self.eini  = -7.0
            self.eend  = -2.0
            self.estep = 0.01

        self.NE = int((self.eend-self.eini)/self.estep)
        print(f'Local density of states constant (cLDOS): {self.cLDOS:.4f}')
        print(f'Fermi enegy: {self.Efermi:.2f} eV')
        print(f'Scanning energy range from {self.eini:.2f} to {self.eend:.2f} eV.')
        print(f'Stepsize of {self.estep:.2f} eV.\n')
        # define some constant terms
        uL = np.eye(self.nL)
        uR = np.eye(self.nR)
        # Greeen's functions of the left and right electrodes approximated
        # within the wide-band limit.
        self.gL= -1j*pi*self.cLDOS*uL
        self.gR= -1j*pi*self.cLDOS*uR

    @PRTiming("transmission")
    def transmission(self,Scc,Fcc,SLC,HLC,SRC,HRC):
        """Evaluating transmission from Green's functions.

        Parameters:

        Scc    : Overlap matrix of the central region.
        Fcc    : Fock matrix of the central region.
        SLC    : Coupling block of the overlap matrix for the left electrode
                 with the central region.
        SRC    : Coupling block of the overlap matrix for the right electrode
                 with the central region.
        HLC    : Coupling block of the Fock for the left electrode with the
                 central region.
        HRC    : Coupling block of the Fock for the right electrode with the
                 central region.
        gL     : Green's function matrix of the isolated left electrode.
        gR     : Green's function matrix of the isolated right electrode.
        nC     : Number of basis functions in the central region.

        Return:

        spT    : Percentage of spin polarization transmission as a function
                 of E-Efermi.
        """
        # pylint: disable=too-many-arguments
        # pylint: disable=too-many-locals

        NE = self.NE
        nC = self.nC
        Efermi = self.Efermi
        spT=[]
        Iu=[[1,0],[0,0]]
        Id=[[0,0],[0,1]]
        Guu=np.zeros((nC,nC),dtype=complex)
        Gud=np.zeros((nC,nC),dtype=complex)
        Gdu=np.zeros((nC,nC),dtype=complex)
        Gdd=np.zeros((nC,nC),dtype=complex)
        f = open('transmission.txt','w')
        f.write('# E-Efermi, percentage polarization,Tuu,Tdd,Tdu,Tud \n')
        print('# E-Efermi, percentage polarization,Tuu,Tdd,Tdu,Tud')
        for sc in range(NE):
            E = self.eini + sc*self.estep
            sL1=E*SLC-HLC
            sigL1 = sL1.T @ self.gL @ sL1
            sR1=E*SRC-HRC
            sigR1 = sR1.T @ self.gR @ sR1
            GammaLu=-2*np.imag(sigL1)
            GammaRu=-2*np.imag(sigR1)
            GammaLd = GammaLu
            GammaRd = GammaRu
            sigL=np.kron(Iu,sigL1)+np.kron(Id,sigL1)
            sigR=np.kron(Iu,sigR1)+np.kron(Id,sigR1)

            Ginv=E*Scc-Fcc-sigL-sigR
            G=np.linalg.inv(Ginv) # symmetric
            for i in range(nC):
                for j in range(nC):
                    Guu[i,j]=G[i,j]
                    Gdu[i,j]=G[i+nC,j]
                    Gud[i,j]=G[i,j+nC]
                    Gdd[i,j]=G[i+nC,j+nC]

            Tuu=np.trace(GammaLu @ Guu @ GammaRu @ np.conj(Guu.T))
            Tdd=np.trace(GammaLd @ Gdd @ GammaRd @ np.conj(Gdd.T))

            Tdu=np.trace(GammaLu @ Gdu @ GammaRd @ np.conj(Gdu.T))
            Tud=np.trace(GammaLd @ Gud @ GammaRu @ np.conj(Gud.T))

            spin = (Tuu+Tdu-Tud-Tdd)/(Tuu+Tdu+Tud+Tdd)
            spT.append(spin)
            str1 = (f'{E-Efermi:6.2f} {spT[sc].real*100:11.8f} {Tuu.real:12.5e} ')
            str2 = (f'{Tdd.real:12.5e} {Tdu.real:12.5e} {Tud.real:12.5e}\n')
            str2p = (f'{Tdd.real:12.5e} {Tdu.real:12.5e} {Tud.real:12.5e}')
            f.write(str1+str2)
            print(str1+str2p)
            if self.debug:
                sys.exit()
    
    @PRTiming("radical_transmission")
    def radical_transmission(self):
        """Electron spin polarization transmission calculation for a radical system
        """
        # pylint: disable=too-many-locals
        # pylint: disable=too-many-branches
        # pylint: disable=too-many-statements

        print('Electron transport w/o spin-orbit coupling for a radical system')
        print('---------------------------------------------------------------')
        nC = self.nC
        nL = self.nL
        nR = self.nR

        if self.load_data is True:
            if isfile('overlap.txt'):
                Sf=np.loadtxt('./overlap.txt') # overlap maxtrix
            else:
                print('overlap.txt not found')
                sys.exit()

            if self.fock_calc is True:
                # load  alpha MO energies, and MO coeffs
                if(isfile('moeigv.txt') and isfile('mocoff.txt')):
                    Eig=np.loadtxt('./moeigva.txt') # MO energies
                    Cmo=np.loadtxt('./mocoffa.txt') # MO coeffs
                    Eigm = np.diagflat(Eig) # converted to 2d diagonal
                    # scalar Fock matrix
                    Ff=Sf@Cmo@Eigm@Cmo.T@Sf
                else:
                    print('moeigva.txt and mocoffa.txt not found')
                    sys.exit()
                # load  beta MO energies, and MO coeffs
                if(isfile('moeigvb.txt') and isfile('mocoffb.txt')):
                    Eigb=np.loadtxt('./moeigvb.txt') # MO energies
                    Cmob=np.loadtxt('./mocoffb.txt') # MO coeffs
                    Eigmb = np.diagflat(Eigb) # converted to 2d diagonal
                    # scalar Fock matrix
                    Ffb=Sf@Cmob@Eigmb@Cmob.T@Sf
                else:
                    print('moeigva.txt and mocoffa.txt not found')
                    sys.exit()

            else:
                if(isfile('fock.txt') and isfile('fockb.txt')):
                    Ff=np.loadtxt('./fock.txt')
                    Ffb=np.loadtxt('./fockb.txt')
                    print('alpha and beta Fock matrices loaded')
                else:
                    print('fock.txt and fockb.txt not found')
                    sys.exit()

        else:
            Ff = self.Fa
            Ffb = self.Fb
            Sf  = self.Sf

        Fcc=np.zeros((nC*2,nC*2),dtype=float)
        Scc=np.zeros((nC*2,nC*2),dtype=float)
        for i in range(nC):
            for j in range(nC):
                Fcc[i,j] = Ff[i+nL,j+nL]
                Scc[i,j] = Sf[i+nL,j+nL]
                Fcc[i+nC,j+nC] = Ffb[i+nL,j+nL]
                Scc[i+nC,j+nC] = Sf[i+nL,j+nL]
        Ff1=Ff*27.21165
        Sf1=Sf
        Fcc=Fcc*27.21165
        HLC1=np.zeros((nL,nC),dtype=float)
        SLC1=np.zeros((nL,nC),dtype=float)
        HRC1=np.zeros((nR,nC),dtype=float)
        SRC1=np.zeros((nR,nC),dtype=float)
        for i in range(nL):
            for j in range(nC):
                SLC1[i,j]= Sf1[i,j+nL]
                HLC1[i,j]= Ff1[i,j+nL]

        for i in range(nR):
            for j in range(nC):
                SRC1[i,j]= Sf1[i+nL+nC,j+nL]
                HRC1[i,j]= Ff1[i+nL+nC,j+nL]
        self.transmission(Scc,Fcc,SLC1,HLC1,SRC1,HRC1)

    @PRTiming("ciss_transmission")
    def ciss_transmission(self):
        """Electron chiral-induced spin selectivity (ciss) transmission
        """
        # pylint: disable=too-many-locals
        # pylint: disable=too-many-branches
        # pylint: disable=too-many-statements

        debug=self.debug
        nC = self.nC
        nL = self.nL
        nR = self.nR
        nmo = self.nmo
        if self.load_data is True:
            # loading MO energies, real MO and imginary coeffs
            if(isfile('mocoff_so_real.txt') and isfile('mocoff_so_imag.txt')
                and isfile('moeigv_so.txt')):
                print('Loading SOC MO energies and coeffs ...')
                VscaR=np.loadtxt('./mocoff_so_real.txt') # real MO coeffs
                VscaI=np.loadtxt('./mocoff_so_imag.txt') # img MO coeffs
                VscaE=np.loadtxt('./moeigv_so.txt') # MO energies
            else:
                print('Files for SOC MO energies and coeffs not found')
                sys.exit()
            Emat = np.diagflat(VscaE) # converted to 2d diagonal
            Vlsoc=VscaR+1j*VscaI
            Fsca= Vlsoc @ Emat @ (np.conj(Vlsoc)).T
            #Fsca_matmul=np.matmul(np.matmul(Vlsoc,Emat),np.conj(Vlsoc))

            # load overlap, MO energies, and MO coeffs
            if(isfile('overlap.txt') and isfile('moeigv.txt')
                and isfile('mocoff.txt')):
                print('Loading scalar overlap matrix, MO energies, and coeffs ...')
                Sf=np.loadtxt('./overlap.txt') # overlap maxtrix
                Eig=np.loadtxt('./moeigv.txt') # MO energies
                Cmo=np.loadtxt('./mocoff.txt') # MO coeffs
            else:
                print('Files scalar overlap matrix, MO energies, and coeffs not found')
                sys.exit()
            Eigm = np.diagflat(Eig) # converted to 2d diagonal
            # scalar Fock matrix
            Ff=Sf@Cmo@Eigm@Cmo.T@Sf
            I2=np.eye(2)
            Sf2=np.kron(I2,Sf)
            # SOC Fock matrix
            Fsfo=Sf2 @ Fsca @ Sf2
            if debug:
                dprint_mat(Emat,name='Emat')
                print(f'type(Vlsoc) {type(Vlsoc)}')
                print(f'Vlsoc.dtype) {Vlsoc.dtype}')
                dprint_mat(Vlsoc,name='Vlsoc')
                dprint_mat(VscaR,name='VscaR')
                dprint_mat(VscaI,name='VscaI')
                dprint_mat(Fsca,name='Fsca')
                dprint_mat(Sf2,name='Sf2')
                dprint_mat(Fsfo,name='Fsfo')
        else:
            Sf = self.Sf
            Ff = self.Fs
            Fsfo=self.Fsoc
            I2=np.eye(2)
            Sf2=np.kron(I2,Sf)

        Fcc=np.zeros((nC*2,nC*2),dtype=complex)
        Scc=np.zeros((nC*2,nC*2),dtype=complex)
        for i in range(nC):
            for j in range(nC):
                Fcc[i,j] = Fsfo[i+nL,j+nL]
                Scc[i,j] = Sf2[i+nL,j+nL]
                Fcc[i+nC,j] = Fsfo[i+nL+nmo,j+nL]
                Scc[i+nC,j] = Sf2[i+nL+nmo,j+nL]
                Fcc[i,j+nC] = Fsfo[i+nL,j+nL+nmo]
                Scc[i,j+nC] = Sf2[i+nL,j+nL+nmo]
                Fcc[i+nC,j+nC] = Fsfo[i+nL+nmo,j+nL+nmo]
                Scc[i+nC,j+nC] = Sf2[i+nL+nmo,j+nL+nmo]
        Ff1=Ff*27.21165
        Sf1=Sf
        Fcc=Fcc*27.21165
        HLC1=np.zeros((nL,nC),dtype=float)
        SLC1=np.zeros((nL,nC),dtype=float)
        HRC1=np.zeros((nR,nC),dtype=float)
        SRC1=np.zeros((nR,nC),dtype=float)
        for i in range(nL):
            for j in range(nC):
                SLC1[i,j]= Sf1[i,j+nL]
                HLC1[i,j]= Ff1[i,j+nL]

        for i in range(nR):
            for j in range(nC):
                SRC1[i,j]= Sf1[i+nL+nC,j+nL]
                HRC1[i,j]= Ff1[i+nL+nC,j+nL]
        self.transmission(Scc,Fcc,SLC1,HLC1,SRC1,HRC1)

