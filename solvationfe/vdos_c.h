#ifndef VDOS_C_H
#define VDOS_C_H

#pragma message("Including vdos_c.h with t_MDinfo and t_residueList")

typedef struct
{
    int bondAtomIndices[2];
    int sat0AtomIndices[4];
    int nSat0;
    int sat3AtomIndices[4];
    int nSat3;
    double sumInertia;
    double sumLogInertia;
    float *wOmegaBuffer;
} t_rotBond;

typedef struct
{
    /* constant data */
    int nAtoms;
    float resMass;
    int nRotBonds;
    float *atomMasses;
    float *atomMassesX3;
    t_rotBond *rotBonds;
    /* transient data */
    float *atomCrdMol;
    float *atomVelMol;
    float inertia[3];
    float prevAxes[9];
    float angMomLab[3];
    float angMomMol[3];
    float omegaMol[3];
    float omegaLab[3];
    float *COMposBuffer;
    float *COMposCurrent;
    float *COMposCorrRef;
    float *COMvelBuffer;
    float *COMvelCurrent;
    float *COMvelCorrRef;
    float *wOmegaBuffer;
    float *wOmegaCurrent;
    float *wOmegaCorrRef;
    /* accumulating data */
    double sumInertia[3];
    double sumLogInertia[3];
    double *totCorr;
    double *trCorr;
    double *rotCorr;
    double *rotBondCorr;
    int propCnt;
    int corrCnt;
    /* offset in MD info arrays */
    int offset;
} t_residue;

typedef struct
{
    t_residue *residues;
    int nResidues;
    double *totCorr;
    double *trCorr;
    double *rotCorr;
    double *rotBondCorr;
    double inertia[3];
    double logInertia[3];
    double rotBondInertia;
    double logRotBondInertia;
} t_residueList;

typedef struct
{
    int frameSize;
    float *atomCrd;
    int bufferSize;
    float *atomVelBuffer;
    float *atomVelCurrent;
    float *atomVelCorrRef;
    int nCorr;
    int nAtomsSel;
} t_MDinfo;

typedef struct
{
    int *indices;
    long long rank;
} t_dih;

int allocResidueList(t_residueList *residueList, int nRes, double *totCorr, double *trCorr, double *rotCorr, double *rotBondCorr);
int allocMDinfo(t_MDinfo *MDinfo, int nCorr, int nAtomsSel);
int allocResidue(t_residue *res, int nAtoms, float *atomMasses, float resMass, int nAtomsSel, int nCorr, double *totCorr, double *trCorr, double *rotCorr, double *rotBondCorr);
int setArrayIndexOffsets(t_residueList *residueList, t_MDinfo *MDinfo);
int computeTotalVACF(t_residue *res, t_MDinfo *MDinfo);
int computeCOM(int bufferShift, t_residue *res, t_MDinfo *MDinfo);
int computeTransVACF(int bufferShift, t_residue *res, t_MDinfo *MDinfo);
int subtractCOM(int bufferShift, t_residue *res, t_MDinfo *MDinfo);
int computeRotation(int bufferShift, t_residue *res, t_MDinfo *MDinfo);
int computeRotVACF(int bufferShift, t_residue *res, t_MDinfo *MDinfo);
int subtractRot(t_residue *res);
int computeRotBonds(int bufferShift, t_residue *res);
int computeRotBondCorr(int bufferShift, t_residue *res, t_MDinfo *MDinfo);
int processStep(int tStep, t_MDinfo *MDinfo, float *crds, float *vels, t_residueList *residueList, int nCorr);
int postProcess(t_residueList *residueList, int nCorr);
int getRotBonds(t_residueList *sets, int nSets, int *dihedAtomIndices, int nDih, int *resAtomIdxRange, int nCorr);

#endif // VDOS_C_H