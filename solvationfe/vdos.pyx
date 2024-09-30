
from scipy.fft import fft
import numpy as np
cimport numpy as np
cimport vdos
import logging

logger = logging.getLogger(__name__)

cdef class VDoS:

    cdef int nRes
    cdef int nCorr
    cdef vdos.t_residueList residueList
    cdef vdos.t_MDinfo MDinfo
    cdef float [:, :] atomPositions
    cdef float [:, :] atomVelocities

    cdef double [:] totCorr 
    cdef double [:] trCorr
    cdef double [:] rotCorr
    cdef double [:] rotBondCorr

    cdef double [:, :] resTotCorrs
    cdef double [:, :] resTrCorrs
    cdef double [:, :] resRotCorrs
    cdef double [:, :] resRotBondCorrs

    cdef object sel

    cdef float [:] tau
    cdef float [:] wavenumber
    cdef double [:, :] totVACF
    cdef double [:, :] totVDoS
    cdef double [:, :, :] trVACF
    cdef double [:, :, :] trVDoS
    cdef double [:, :, :] rotVACF
    cdef double [:, :, :] rotVDoS
    cdef double [:, :] rotBondVACF
    cdef double [:, :] rotBondVDoS

    # def __cinit__(self):
        

    def __init__(self,sel,nCorr):
        self.sel = sel
        self.nCorr = nCorr
        self.nRes = sel.residues.n_residues
        # time and frequency axes for correlation functions and VDoS
        self.tau = np.zeros(nCorr, dtype = np.float32)
        self.wavenumber = np.zeros(nCorr, dtype = np.float32)
        # numpy arrays for average (index 0) and per residue VACF & VDoS (in 3D-2PT for each voxel)
        # -> total VACF / VDoS (1D)
        # -> translation VACF / VDoS (3D)
        # -> rotation VACF / VDoS (3D)
        # -> rotatable bonds VACF / VDoS (1D)
        self.totVACF     = np.zeros((self.nRes+1,    nCorr), dtype = np.float64)
        self.totVDoS     = np.zeros((self.nRes+1,    nCorr), dtype = np.float64)
        self.trVACF      = np.zeros((self.nRes+1, 3, nCorr), dtype = np.float64)
        self.trVDoS      = np.zeros((self.nRes+1, 3, nCorr), dtype = np.float64)
        self.rotVACF     = np.zeros((self.nRes+1, 3, nCorr), dtype = np.float64)
        self.rotVDoS     = np.zeros((self.nRes+1, 3, nCorr), dtype = np.float64)
        self.rotBondVACF = np.zeros((self.nRes+1,    nCorr), dtype = np.float64)
        self.rotBondVDoS = np.zeros((self.nRes+1,    nCorr), dtype = np.float64)

        # Alloc correlation arrays (passed to C codes)
        self.totCorr = np.zeros(self.nCorr)
        self.trCorr = np.zeros(self.nCorr * 3)
        self.rotCorr = np.zeros(self.nCorr * 3)
        self.rotBondCorr = np.zeros(self.nCorr)

        # Alloc correlation arrays for each residue
        self.resTotCorrs = np.zeros((self.nRes, self.nCorr))
        self.resTrCorrs = np.zeros((self.nRes, self.nCorr * 3))
        self.resRotCorrs = np.zeros((self.nRes, self.nCorr * 3))
        self.resRotBondCorrs = np.zeros((self.nRes, self.nCorr))

        # Alloc pos, vel arrays to allow copying later on
        # self.atomPositions = np.zeros((self.sel.n_atoms, 3, 3), dtype = np.float32)
        # self.atomVelocities = np.zeros((self.sel.n_atoms, 3, 3), dtype = np.float32)
        
        # initialize residueList and MDinfo
        self.prep()


    def prep(self):
        '''construct list of residues for selection'''
        n_atoms = self.sel.n_atoms

        error = vdos.allocResidueList(
            &self.residueList,
            self.nRes,
            &self.totCorr[0],
            &self.trCorr[0],
            &self.rotCorr[0],
            &self.rotBondCorr[0],
        )

        error = vdos.allocMDinfo(
            &self.MDinfo,
            self.nCorr,
            n_atoms
        )

        
        cdef float [:] atomMasses

        r = 0
        for res in self.sel.residues:
            # Don't copy- assign, since each res has diff num atoms
            atomMasses = res.atoms.masses.astype(np.float32)
            error = vdos.allocResidue(
                &self.residueList.residues[r],
                len(res.atoms),
                # ptr to first element of atomMasses
                &atomMasses[0],
                res.mass.astype('float32'),
                n_atoms,
                self.nCorr,
                &self.resTotCorrs[r][0],
                &self.resTrCorrs[r][0],
                &self.resRotCorrs[r][0],
                &self.resRotBondCorrs[r][0]
            )
            if error != 0:
                print(f'ERROR reported by allocResidue\n')
            r += 1
        vdos.setArrayIndexOffsets(&self.residueList, &self.MDinfo)
        cdef int [:, :] dihed = self.sel.intra_dihedrals.indices.astype(np.int32)
        cdef int [:, :] resAtomRangeList = np.zeros((len(self.sel.residues),2), dtype = np.int32)
        i = 0
        for res in self.sel.residues:
            resAtomRangeList[i][0] = np.amin(res.atoms.indices)
            resAtomRangeList[i][1] = np.amax(res.atoms.indices)
            i += 1

        error = vdos.getRotBonds(
            &self.residueList,
            self.nRes,
            &dihed[0][0],
            len(dihed),
            &resAtomRangeList[0][0],
            self.nCorr
        )

        if error != 0:
            print(f'ERROR reported by \'clib.getRotBonds\'\n')
    
    cpdef processStep(self, int tStep, float time):
        if tStep < self.nCorr:
            self.tau[tStep] = time
        self.atomPositions = np.asarray(self.sel.atoms.positions).view(np.float32)
        self.atomVelocities = np.asarray(self.sel.atoms.velocities).view(np.float32)
        error = vdos.processStep(
            tStep,
            &self.MDinfo,
            &self.atomPositions[0][0],
            &self.atomVelocities[0][0],
            &self.residueList,
            self.nCorr
        )

    def postProcess(self):
        '''post-process correlation functions'''
        vdos.postProcess(&self.residueList, self.nCorr)

        period = (self.tau[1] - self.tau[0]) * (2 * self.nCorr - 1)
        wn0 = (1.0 / period) * 33.35641
        self.wavenumber = np.arange(0,self.nCorr, dtype=np.float32) * wn0
        self.totVACF[0] = self.totCorr[0:self.nCorr]
        self.symFT(self.totVACF[0], self.totVDoS[0])

        for i in range(3):
            self.trVACF[0][i] = self.trCorr[i:3*self.nCorr:3]
            self.symFT(self.trVACF[0][i], self.trVDoS[0][i])
            self.rotVACF[0][i] = self.rotCorr[i:3*self.nCorr:3]
            self.symFT(self.rotVACF[0][i], self.rotVDoS[0][i])
        self.rotBondVACF[0] = self.rotBondCorr[0:self.nCorr]
        self.symFT(self.rotBondVACF[0], self.rotBondVDoS[0])
        for i in range(1,self.nRes+1):
            self.totVACF[i] = self.resTotCorrs[i-1,0:self.nCorr]
            self.symFT(self.totVACF[i], self.totVDoS[i])
            for j in range(3):
                self.trVACF[i][j] = self.resTrCorrs[i-1,j:3*self.nCorr:3]
                self.symFT(self.trVACF[i][j], self.trVDoS[i][j])
                self.rotVACF[i][j] = self.resRotCorrs[i-1, j:3*self.nCorr:3]
                self.symFT(self.rotVACF[i][j], self.rotVDoS[i][j])
            self.rotBondVACF[i] = self.resRotBondCorrs[i-1, 0:self.nCorr]
            self.symFT(self.rotBondVACF[i], self.rotBondVDoS[i])

    def symFT(self,input,output):
        '''symmetrize, Fourier transform and return first half of spectrum'''
        data = np.array(input, dtype = np.float32)
        nData=len(data)
        tmp = np.zeros(2 * nData - 1, dtype = np.float32)
        for i in range(nData):
            tmp[i] = data[i]
        for i in range(1,nData):
            j = 2 * nData - i - 1
            tmp[j] = tmp[i]
        tmp = fft(tmp)
        for i in range(nData):
            output[i] = tmp[i].real
    
    def outputGeometry(self,outputFileName):
        totMass = 0.0
        for i in range(self.nRes):
            totMass += self.residueList.residues[i].resMass
        averMass = totMass / self.nRes
        averInertia = np.zeros(3)
        averLogInertia = np.zeros(3)
        for i in range(3):
            averInertia[i] = self.residueList.inertia[i]
            averLogInertia[i] = self.residueList.logInertia[i]
        averRotBondCount = 0
        for i in range(self.nRes):
            averRotBondCount += self.residueList.residues[i].nRotBonds
        averRotBondCount /= self.nRes
        if averRotBondCount > 0:
            averRotBondInertia = self.residueList.rotBondInertia / averRotBondCount
            averLogRotBondInertia = self.residueList.logRotBondInertia / averRotBondCount
        else:
            averRotBondInertia = 0.0
            averLogRotBondInertia = 0.0
        outFile = open(outputFileName,"w")
        outFile.write("#Average Mass (g/mol):\n%20.6e\n" % averMass)
        outFile.write("#Average Inertia (g/mol*A^2):\n%20.6e %20.6e %20.6e\n" % (averInertia[0],averInertia[1],averInertia[2]))
        outFile.write("#Average Log(Inertia):\n%20.6e %20.6e %20.6e\n" % (averLogInertia[0],averLogInertia[1],averLogInertia[2]))
        outFile.write("#Average Rotatable Bond Inertia (g/mol*A^2): (%d)\n%20.6e\n" % (averRotBondCount,averRotBondInertia))
        outFile.write("#Average Rotatable Bond Log(Inertia): (%d)\n%20.6e\n" % (averRotBondCount,averLogRotBondInertia))
        outFile.close()
    
    def outputVACF(self,outputFileName):
        outFile = open(outputFileName,"w")
        outFile.write("%-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s\n" % ("#Time (ps)","Translation x","Translation y","Translation z","Rotation x","Rotation y","Rotation z","Rotatable Bonds","Total"))
        for i in range(self.nCorr):
            outFile.write("%-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e\n" % (self.tau[i],self.trVACF[0][0][i],self.trVACF[0][1][i],self.trVACF[0][2][i],self.rotVACF[0][0][i],self.rotVACF[0][1][i],self.rotVACF[0][2][i],self.rotBondVACF[0][i],self.totVACF[0][i]))
        outFile.close()
        outFile.close()

    def outputVDoS(self,outputFileName):
        outFile = open(outputFileName,"w")
        outFile.write("%-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s\n" % ("#Wavenumber (cm^-1)","Translation x","Translation y","Translation z","Rotation x","Rotation y","Rotation z","Rotatable Bonds","Total"))
        for i in range(self.nCorr):
            outFile.write("%-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e\n" % (self.wavenumber[i],self.trVDoS[0][0][i],self.trVDoS[0][1][i],self.trVDoS[0][2][i],self.rotVDoS[0][0][i],self.rotVDoS[0][1][i],self.rotVDoS[0][2][i],self.rotBondVDoS[0][i],self.totVDoS[0][i]))
        outFile.close()
