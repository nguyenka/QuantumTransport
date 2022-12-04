# pylint: disable=invalid-name
"""
   Module for parsing data from the output files generated electronic
   structure programs(ADF, Gaussian, Turbomole) to calculate electron
   quantum transport (QT) properties.
"""
import argparse
import json
from os.path import isfile
import re
import sys
import numpy as np
from munch import Munch
from cclib.parser import utils
from cclib.parser.utils import PeriodicTable
from cclib.parser import ccopen
from cclib.parser.turbomoleparser import Turbomole

def write_basis_data(left_basis,central_basis,right_basis,total_basis):
    """
    function for writting basis functions of the left, central, right regions.
    """

    basis_data = Munch(left_basis = left_basis, central_basis= central_basis,\
                 right_basis = right_basis, total_basis = total_basis)
    with open('basis_data.json','w') as f:
        json.dump(basis_data,f,indent=4)
        print(f"basis data written in {'basis_data.json'}")

def partition_mol(atomcoords,atomnos):
    """
    function for partitioning a system into the left, central, right regions.
    """

    t = PeriodicTable()
    elem_list = []
    for index, _ in enumerate(atomcoords[-1]):
        elem = t.element[atomnos[index]]
        elem_list.append(elem)
    latoms = 0
    catoms = 0
    for i in range(len(elem_list)):
        elem = elem_list[i].lower()
        if i == 0:
            latoms += 1
        if(i > 0 and elem == 'au' and elem_list[i-1].lower() == 'au' \
           and catoms == 0):
            latoms += 1
        if(i > 0 and elem != 'au'):
            catoms += 1
    ratoms = len(elem_list) - latoms - catoms
    print('Partitionng the system into three regions based on Au atoms:')
    print(f'left electrode: {latoms} atoms')
    print(f'central region: {catoms} atoms')
    print(f'right electrode: {ratoms} atoms')
    print('Note:')
    print('The first and last sets of heavy atoms are assummed be on left and right')
    print('electrodes, respectively.\n')
    return latoms, catoms

class ADFparse:
    """
    Class for extract data generate by ADF
    """
    def __init__(self,lout = '../left/left.out',rout = '../right/right.out',\
                cout= '../junc/junc.out'):
        """
        Parameters:

        inputfile: ADF outputs (for left, right, and central parts), ljr1c_overlap.dat,
                   ljr1c_fock.dat,ljr2c_moeigv.dat,ljr2c_mocoff.dat
        """

        self.lout = lout
        self.rout = rout
        self.cout = cout

    def parse_basis(self):
        """
        get basis functions for left,right,central part of the system from ADF outputs
        write data in basis_data.json
        """

        def get_basis(out):
            nF = 0
            fD = open(out,"r")
            for k in fD:
                if  "Nr. of SFOs " in k:
                    nF = int(k.split(':')[1])
                    break
            fD.close()
            return nF
        nL = get_basis(self.lout)
        nC = get_basis(self.cout)
        nR = get_basis(self.lout)
        nmo = nL + nC + nR
        write_basis_data(nL,nC,nR,nmo)
        return nL, nC, nR, nmo

    def overlap_scalar(self,nmo):
        """ Get  Overlap-matrix (Fragment orbital as basis set) """

        SA=open('overlap.dat','r')
        SAw=open('overlap.txt','w')
        lines = [line for line in SA]
        ln = len(lines)
        k=3
        sa={}
        io=0
        while k < ln:
            line=lines[k]
            nele=len(line.split())
            for i in range(nele):
    #           sa[io]=float(line.split()[i])
                sa[io]=' '+line.split()[i]
                io = io +1
            k=k+1

        if nmo*(nmo+1)/2 != io:
            print(' Frag Orb numbers donot match')
            sys.exit()
        sfo=[[' ' for i in range(nmo)] for j in range(nmo)]
        k=0
        for m in range(nmo):
            for n in range(m+1):
                sfo[m][n]=sa[k]
                if n < m:
                    sfo[n][m]=sa[k]
                k = k + 1

        for m in range(nmo):
            SAw.writelines(sfo[m])
            SAw.writelines("\n")
        print(f"{'overlap.txt'} is printed out!")

    def fock_scalar(self,nmo):
        """ get one-component scalar relativistic Fock-matrix """

        FA=open('fock.dat','r')
        FAw=open('fock.txt','w')
        lines = [line for line in FA]
        ln = len(lines)
        k=3
        fa={}
        io=0
        while k < ln:
            line=lines[k]
            nele=len(line.split())
            for i in range(nele):
                fa[io]=' '+line.split()[i]
                io = io +1
            k=k+1
        if nmo*(nmo+1)/2 != io:
            print(' Frag Orb number donot match')
            sys.exit()
        ffo=[[' ' for i in range(nmo)] for j in range(nmo)]
        k=0
        for m in range(nmo):
            for n in range(m+1):
                ffo[m][n]=fa[k]
                if n < m:
                    ffo[n][m]=fa[k]
                k = k + 1
        for m in range(nmo):
            FAw.writelines(ffo[m])
            FAw.writelines("\n" )
        print(f"{'fock.txt'} is printed out!")

    def MOvec_spinorb(self,nmo):
        """ Get MO eigenvalue matrix."""
        # pylint: disable=too-many-locals

        MER=open("moeigv_so.dat","r")
        MEW=open("moeigv_so.txt","w")

        lines = [line for line in MER]
        ln = len(lines)
        k=3
        moe=[' ' for i in range(nmo*2)]
        io=0
        while  k < ln:
            line=lines[k]
            nele=len(line.split())
            for i in range(nele):
                moe[io]=' '+line.split()[i]
                io = io +1
            k=k+1

        MEW.writelines(moe)
        MEW.writelines("\n" )
        print(f"{'_moeigv_so.txt'} is printed out!")

        ### get MO coefficient matrix
        MLR=open("mocoff_so.dat","r")
        MLWr=open("mocoff_so_real.txt","w")
        MLWi=open("mocoff_so_imag.txt","w")
        im=nmo*nmo*4

        lines = [line for line in MLR]
        ln = len(lines)
        k=3
        vlow = [' ' for i in range(nmo*nmo*8)]
        io=0
        while k < ln:
            line=lines[k]
            nele=len(line.split())
            for i in range(nele):
                vlow[io]=' ' + line.split()[i]
                io = io +1
            k=k+1

        for m in range(nmo*2):
            MLWr.writelines(vlow[m*nmo*2:(m+1)*nmo*2])
            MLWr.writelines("\n" )
        print(f"{'mocoff_so_real.txt'} is printed out!")

        for m in range(nmo*2):
            MLWi.writelines(vlow[im+m*nmo*2:(m+1)*nmo*2+im])
            MLWi.write(" \n" )
        print(f"{'mocoff_so_imag.txt'} is printed out!")

class GaussianParse:
    """
    Class for extract data from Gaussian output
    """
    def __init__(self,latoms = 0,catoms =0 ,gout = 'dft.out',fock_parse = True, \
                save_mos = False,savefile = True):
        """
        Parameters:

        latoms: number of atoms in the left electrode.
        catoms: number of atoms in the in the central region.
        gout: gaussian output.
        fock_parse: Fock matrix parsing.
        save_mos: save MO data.
        savefile: save all data for QT calculations.
        """
        # pylint: disable=too-many-arguments

        self.unrestricted = False
        self.latoms = latoms
        self.catoms = catoms
        self.gout = gout
        self.fock_parse = fock_parse
        self.save_mos = save_mos
        self.savefile = savefile

    def parse_output(self):

        """
        parse: gaussian output
        exact: basis functions left,central,and right electrodes,overlap matrix,
               mo coeffs and energies
        """
        # pylint: disable=too-many-locals
        # pylint: disable=too-many-statements

        gout = self.gout
        latoms = self.latoms
        catoms = self.catoms
        ccfile = ccopen(gout)
        data = ccfile.parse()
        program = data.metadata['package']
        if program != 'Gaussian':
            print(f"Expecting Gaussian output, got {program} ouput.")
            sys.exit()
        atombasis = data.atombasis
        # partition a molecular junction based on Au atoms
        if(latoms == 0 and catoms ==0):
            latoms, catoms = partition_mol(data.atomcoords,data.atomnos)
        left_basis = 0
        central_basis = 0
        right_basis = 0
        for i in range(data.natom):
            if i + 1 <= latoms:
                left_basis += len(atombasis[i])
            elif  i +1  >= latoms and i +1 <= catoms+latoms:
                central_basis += len(atombasis[i])
            else:
                right_basis += len(atombasis[i])
        total_basis = data.nbasis
        if total_basis != left_basis+central_basis+right_basis:
            print('inconsistency in number of basis functions!')
            sys.exit()
        overlap = data.aooverlaps
        if self.savefile:
            write_basis_data(left_basis,central_basis,right_basis,total_basis)
            np.savetxt('overlap.txt',overlap)
            print(f"overlap matrix written in {'overlap.txt'}")
        if self.fock_parse:
            focka, fockb = self.parse_fock()
        else:
            if len(data.moenergies) == 2:
                amoenergies = utils.convertor(data.moenergies[0], 'eV', 'hartree')
                bmoenergies = utils.convertor(data.moenergies[1], 'eV', 'hartree')
                if self.save_mos:
                    np.savetxt('moeigva.txt',amoenergies)
                    np.savetxt('moeigvb.txt',bmoenergies)
                    np.savetxt('mocoffa.txt',data.mocoeffs[0].T)
                    np.savetxt('mocoffb.txt',data.mocoeffs[1].T)
                    print(f"alpha MO energies written in {'moeigva.txt'}")
                    print(f"beta MO energies written in {'moeigvb.txt'}")
                    print(f"alpha MO coeffs written in {'mocoffa.txt'}")
                    print(f"beta MO coeffs written in {'mocoffb.txt'}")
                amoenergy2d= np.diagflat(amoenergies) # converted to 2d diagonal
                bmoenergy2d= np.diagflat(bmoenergies)
                focka = overlap@data.mocoeffs[0].T@amoenergy2d@data.mocoeffs[0]@overlap
                fockb = overlap@data.mocoeffs[1].T@bmoenergy2d@data.mocoeffs[1]@overlap
                np.savetxt('fock.txt',focka)
                np.savetxt('fockb.txt',fockb)
                print(f"alpha Fock matrix written in {'fock.txt'}")
                print(f"beta Fock matrix written in {'fock.txt'}")
            else:
                moenergies = utils.convertor(data.moenergies, 'eV', 'hartree')
                moenergy2d= np.diagflat(moenergies) # converted to 2d diagonal
                focka = overlap @ data.mocoeffs.T @ moenergy2d @data.mocoeffs @ overlap
                if self.save_mos:
                    np.savetxt('moeigv.txt',moenergies)
                    np.savetxt('mocoff.txt',data.mocoeffs.T)
                    print(f"MO energies written in {'moeigv.txt'}")
                    print(f"MO coeffs written in {'mocoff.txt'}")
                np.savetxt('fock.txt',focka)
                print(f"Fock matrix written in {'fock.txt'}")
        return left_basis,central_basis,right_basis,total_basis,overlap,focka,fockb

    def parse_fock(self):
        """ Parse Fock matrix from Gaussian output."""
        # pylint: disable=too-many-locals
        # pylint: disable=too-many-nested-blocks

        inputfile = open(self.gout,'r')
        lines = inputfile.readlines()
        for i in range(len(lines)):
            if "NBasis" in lines[i]:
                nbasis = int(lines[i].split('=')[1].split()[0])
                #print(f'number of bais functions {nbasis}')
            if "Fock matrix (alpha):" in lines[i] or "Fock matrix (beta):" in lines[i]:
                #if "Fock matrix (beta):" in lines[i]:
                beta = "Fock matrix (beta):" in lines[i]
                fls = i + 1
                #print(f'fock elements start at: {fls}')
                fock = np.zeros((nbasis, nbasis),dtype=float)
                ll = 0
                base = 0
                while base < nbasis:
                    for col in range(int(nbasis/5)):
                        for l in range(nbasis-base):
                            ll += 1
                            fl = ll+fls + col
                            line = lines[fl]
                            #print(f'{line} line# {fl}')
                            parts = line.split()
                            for j in range(len(parts)-1):
                                val = float(parts[j+1].replace("D", "E"))
                                fock[base+j, l+base] = val
                                fock[l+base, base+j] = val
                        base += 5
                        #print(f'base: {base}')
                if beta:
                    fockb = fock
                    if self.savefile:
                        print(f"beta Fock matrix written in {'fockb.txt'}")
                        np.savetxt('fockb.txt',fockb)
                else:
                    focka = fock
                    if self.savefile:
                        print(f"alpha Fock matrix written in {'fock.txt'}")
                        np.savetxt('fock.txt',focka)
        return focka, fockb

class TurbomoleParse:
    """
    Class for extract data from turbmole output and  mos files.
    """
    # pylint: disable=too-many-instance-attributes

    def __init__(self,latoms = 0,catoms =0 ,tout = 'ridft.out',savefile = True):
        """
        Parameters:

        latoms: number of atoms in the left electrode.
        catoms: number of atoms in the in the central region.
        tout: Turbomole output.
        savefile: save all data for QT calculations.
        """

        self.unrestricted = False
        self.latoms = latoms
        self.catoms = catoms
        self.tout = tout
        self.savefile = savefile

    def normalisesym(self, label):
        """Normalise the symmetries used by Turbomole.

        The labels are standardized except for the first character being lowercase.
        """
        # TODO more work could be required, but we don't have any logfiles
        # with non-C1 symmetry.
        return label.capitalize()

    def parse_output(self):
        """ Parse Turbomle output."""
        # pylint: disable=too-many-locals
        # pylint: disable=too-many-branches
        # pylint: disable=too-many-statements

        t19out = self.tout
        print(f"Reading Turbomole output file: {t19out}")
        latoms = self.latoms
        catoms = self.catoms
        ccfile = ccopen(t19out)
        data = ccfile.parse()
        program = data.metadata['package']
        if program != 'Turbomole':
            print(f"Expecting Turbomole output, got {program} ouput.")
            sys.exit()
        t = PeriodicTable()
        atombasis = data.atombasis # required local version of turbomoleparser.py
        # partition a molecular junction based on Au atoms
        if(latoms == 0 and catoms ==0):
            latoms, catoms = partition_mol(data.atomcoords,data.atomnos)
        left_basis = 0
        central_basis = 0
        right_basis = 0
        for index,_ in enumerate(data.atomcoords[-1]):
            elem = t.element[data.atomnos[index]]
            if index + 1 <= latoms:
                left_basis += atombasis[elem.lower()]
            elif index +1  >= latoms and index +1 <= catoms+latoms:
                central_basis += atombasis[elem.lower()]
            else:
                right_basis += atombasis[elem.lower()]
        total_basis = data.nbasis
        if total_basis != left_basis+central_basis+right_basis:
            print('inconsistency in number of basis functions!')
            sys.exit()
        overlap = data.aooverlaps
        if self.savefile:
            write_basis_data(left_basis,central_basis,right_basis,total_basis)
            np.savetxt('overlap.txt',overlap)
            print(f"overlap matrix written in {'overlap.txt'}")
            if isfile('mos'):
                mofile = 'mos'
                moenergy_file = 'moeigv.txt'
                mocoeffs_file = 'mocoff.txt'
                print("Reading scalar mos file ...")
                self.parse_mos(mofile,moenergy_file,mocoeffs_file)
            if isfile('alpha'):
                mofile = 'alpha'
                moenergy_file = 'moeigva.txt'
                mocoeffs_file = 'mocoffa.txt'
                self.parse_mos(mofile,moenergy_file,mocoeffs_file)
                print("Reading scalar alpha file ...")
                self.parse_mos()
            if isfile('beta'):
                mofile = 'beta'
                moenergy_file = 'moeigvb.txt'
                mocoeffs_file = 'mocoffb.txt'
                print("Reading scalar beta file ...")
                self.parse_mos(mofile,moenergy_file,mocoeffs_file)
            if isfile('spinor.r'):
                mofile = 'spinor.r'
                moenergy_file = 'moeigv_so.txt'
                mocoeffs_file = 'mocoff_so_real.txt'
                print("Reading SOC spinor.r file ...")
                self.parse_mos(mofile,moenergy_file,mocoeffs_file)
            if isfile('spinor.i'):
                mofile = 'spinor.i'
                moenergy_file = 'moeigvi_so.txt'
                mocoeffs_file = 'mocoff_so_imag.txt'
                print("Reading SOC spinor.i file ...")
                self.parse_mos(mofile,moenergy_file,mocoeffs_file)
        else:
            if isfile('mos'):
                mofile='mos'
                print("Reading scalar mos file ...")
                print("Calculating scalar Fock matrix ...")
                Eig,Cmo = self.parse_mos(mofile)
                Eigm = np.diagflat(Eig) # converted to 2d diagonal
                # scalar Fock matrix
                Fs=overlap@Cmo@Eigm@Cmo.T@overlap
            if isfile('spinor.r') and isfile('spinor.i'):
                print("Reading SOC spinor.r and spinor.i files ...")
                mofile = 'spinor.i'
                Eigsoc,CmosocI = self.parse_mos(mofile)
                mofile = 'spinor.r'
                Eigsoc,CmosocR = self.parse_mos(mofile)
                print("Calculating SOC Hamiltonian matrix ...\n")
                Eigsocm = np.diagflat(Eigsoc) # converted to 2d diagonal
                Vlsoc=CmosocR+1j*CmosocI
                Fsca= Vlsoc @ Eigsocm @ (np.conj(Vlsoc)).T
                I2=np.eye(2)
                Sf2=np.kron(I2,overlap)
                # SOC Fock matrix
                Fsoc = Sf2 @ Fsca @ Sf2
            return left_basis,central_basis,right_basis,total_basis,overlap,Fs,Fsoc

    def parse_mos(self,mofile='mos',moenergy_file='moeigv.txt',mocoeffs_file= \
                  'mocoff.txt'):
        """
        Modified cclib function to parse mocoeffs and moenergies from files:
        mos, spinor.r,spinor.i, alpha, and beta.
        # $scfmo    scfconv=6   format(4d20.14)
        # # SCF total energy is     -382.3457535740 a.u.
        # #
        #      1  a      eigenvalue=-.97461484059799D+01   nsaos=60
        # 0.69876828353937D+000.32405121159405D-010.87670894913921D-03-.85232349313288D-07
        # 0.19361534257922D-04-.23841194890166D-01-.81711001390807D-020.13626356942047D-02
        # ...
        # ...
        # $end
        """
        # pylint: disable=too-many-locals

        inputfile = open(mofile,'r')
        line = inputfile.readline()
        if (line.startswith('$scfmo') or line.startswith('$uhfmo') or \
            line.startswith('$spinor_real') or line.startswith('$spinor_imag')) and \
            line.find('scfconv') > 0:
            #sys.exit()
            if line.strip().startswith('$uhfmo_alpha'):
                self.unrestricted = True

            # Need to skip the first line to start with lines starting with '#'.
            line = next(inputfile)
            while line.strip().startswith('#') and not line.find('eigenvalue') > 0:
                line = next(inputfile)

            moirreps = []
            moenergies = []
            mocoeffs = []
            mosyms = []

            while not line.strip().startswith('$'):
                number, sym = line.split()[:2]
                number = int(number)
                sym = self.normalisesym(sym)
                moenergy = re.split('=|nsaos',line)[1].strip()
                eigenvalue = utils.float(moenergy)
                #orbital_energy = utils.convertor(eigenvalue, 'hartree', 'eV')
                orbital_energy = eigenvalue

                moenergies.append(orbital_energy)
                mosyms.append(sym)
                moirreps.append((number, sym))

                single_coeffs = []
                nsaos = int(re.split('=',line)[-1].strip())
                #nsaos = int(info.group('count'))
                while len(single_coeffs) < nsaos:
                    line = next(inputfile)
                    single_coeffs.extend(Turbomole.split_molines(line))

                mocoeffs.append(single_coeffs)
                line = next(inputfile)

            max_nsaos = max([len(i) for i in mocoeffs])
            for i in mocoeffs:
                while len(i) < max_nsaos:
                    i.append(np.nan)

            # We now need to sort our orbitals (because Turbomole groups them by symm).
            mos = list(zip(moenergies, mocoeffs, mosyms, moirreps))
            mos.sort(key = lambda mo: mo[0])
            moenergies, mocoeffs, mosyms, moirreps = zip(*mos)
            np_mocoeffs = np.zeros((nsaos,nsaos),dtype=float)
            np_mocoeffs = np.array(mocoeffs)
            if self.savefile:
                #np.savetxt(self.moenergy_file,moenergies,newline=' ')
                np.savetxt(moenergy_file,moenergies,newline=' ')
                print(f"MO energies written in {moenergy_file}")
                #np.savetxt(self.mocoeffs_file,np_mocoeffs.T)
                np.savetxt(mocoeffs_file,np_mocoeffs.T)
                print(f"MO coeffs written in {mocoeffs_file}")
        return moenergies, np_mocoeffs.T

def main():
    """ Main function."""

    description = """
    Prepare input data for electron transport calculations from ADF, Gaussian, or Turbomle
    outputs.
    """
    epilog = """output: transport_params.json, overlap, fock matrices, and/or data files
                        for constructing fock matrix
    """
    parser = argparse.ArgumentParser(
        usage='%(prog)s [options] ',
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=epilog)
    parser.add_argument('filename', metavar='filename', type=str, help='for g16 or t19 \
            output, "none" for ADF')
    parser.add_argument('-la', '--latoms', type=int, default= 0,help='left electrode atoms')
    parser.add_argument('-ca', '--catoms', type=int, default= 0,help='central atoms')
    #parser.add_argument('-ot', '--otype', type=str, default= 'gaussian',
    parser.add_argument('-ot', '--otype', type=str, default= 'turbomole',
                   help='gausian or turbomole')
    args = parser.parse_args()
    if len(sys.argv) == 0:
        parser.print_help()
        sys.exit(1)
    filename = args.filename
    latoms = args.latoms
    catoms = args.catoms
    output_type = args.otype
    if output_type == 'turbomole':
        tparse = TurbomoleParse(latoms=latoms,catoms=catoms,tout=filename)
        tparse.parse_output()
    if output_type == 'gaussian':
        gparse = GaussianParse(latoms=latoms,catoms=catoms,gout=filename)
        gparse.parse_output()
    if output_type == 'adf':
        aparse = ADFparse()
        _,_,_,nmo = aparse.parse_basis()
        aparse.overlap_scalar(nmo)
        aparse.fock_scalar(nmo)
        aparse.MOvec_spinorb(nmo)
if __name__ == "__main__":
    main()
