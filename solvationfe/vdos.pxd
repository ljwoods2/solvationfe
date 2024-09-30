cdef extern from "vdos_c.h":

    # Data structures    
    ctypedef struct t_rotBond:
        int[2] bondAtomIndices
        int[4] sat0AtomIndices
        int nSat0
        int[4] sat3AtomIndices
        int nSat3
        double sumInertia
        double sumLogInertia
        float* wOmegaBuffer  

    ctypedef struct t_residue:
        # constant data
        int nAtoms
        float resMass
        int nRotBonds
        float* atomMasses
        float* atomMassesX3
        t_rotBond* rotBonds  
        # transient data
        float* atomCrdMol
        float* atomVelMol
        float[3] inertia
        float[9] prevAxes
        float[3] angMomLab
        float[3] angMomMol
        float[3] omegaMol
        float[3] omegaLab
        float* COMposBuffer
        float* COMposCurrent
        float* COMposCorrRef
        float* COMvelBuffer
        float* COMvelCurrent
        float* COMvelCorrRef
        float* wOmegaBuffer
        float* wOmegaCurrent
        float* wOmegaCorrRef

        # Accumulating data
        double[3] sumInertia
        double[3] sumLogInertia
        double* totCorr
        double* trCorr
        double* rotCorr
        double* rotBondCorr
        int propCnt
        int corrCnt

        # Offset in arrays above
        int offset

    ctypedef struct t_residueList:
        t_residue* residues  # Pointer to t_residue
        int nResidues
        double* totCorr
        double* trCorr
        double* rotCorr
        double* rotBondCorr
        double[3] inertia
        double[3] logInertia
        double rotBondInertia
        double logRotBondInertia

    ctypedef struct t_MDinfo:
        int frameSize  # Number of floats per frame in array
        float* atomCrd  # Single set of atomic coordinates for selection
        int bufferSize  # Number of floats in atomVelBuffer
        float* atomVelBuffer  # nCorr sets of atomic velocities for selection
        float* atomVelCurrent  # Pointer to current set of atomic velocities in buffer
        float* atomVelCorrRef  # Pointer to first set of atomic velocities in current buffer
        int nCorr  # nCorr sets of atomic velocities for selection
        int nAtomsSel  # Number of atoms in selection

    # Methods
    int allocResidueList(t_residueList* residueList, int nRes, double* totCorr, double* trCorr, double* rotCorr, double* rotBondCorr)
    int allocMDinfo(t_MDinfo* MDinfo, int nCorr, int nAtomsSel)
    int allocResidue(t_residue* res, int nAtoms, float* atomMasses, float resMass, int nAtomsSel, int nCorr, double* totCorr, double* trCorr, double* rotCorr, double* rotBondCorr)
    int setArrayIndexOffsets(t_residueList* residueList, t_MDinfo* MDinfo)
    int computeTotalVACF(t_residue* res, t_MDinfo* MDinfo)
    int computeCOM(int bufferShift, t_residue* res, t_MDinfo* MDinfo)
    int computeTransVACF(int bufferShift, t_residue* res, t_MDinfo* MDinfo)
    int subtractCOM(int bufferShift, t_residue* res, t_MDinfo* MDinfo)
    int computeRotation(int bufferShift, t_residue* res, t_MDinfo* MDinfo)
    int computeRotVACF(int bufferShift, t_residue* res, t_MDinfo* MDinfo)
    int subtractRot(t_residue* res)
    int computeRotBonds(int bufferShift, t_residue* res)
    int computeRotBondCorr(int bufferShift, t_residue* res, t_MDinfo* MDinfo)
    int processStep(int tStep, t_MDinfo* MDinfo, float* crds, float* vels, t_residueList* residueList, int nCorr)
    int postProcess(t_residueList* residueList, int nCorr)
    int getRotBonds(t_residueList* sets, int nSets, int* dihedAtomIndices, int nDih, int* resAtomIdxRange, int nCorr)


