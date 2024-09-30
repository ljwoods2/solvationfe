from solvationfe import VDoS
import MDAnalysis as mda
import tqdm

TOPOL = "/scratch/mheyden1/POPC/run-NVE.tpr"
TRAJ = "/scratch/mheyden1/POPC/run-NVE.trr"
u = mda.Universe(TOPOL,TRAJ)
sel = u.select_atoms("resname POPC")
vdos = VDoS(sel,200)

tStep = 0
for ts in tqdm.tqdm(u.trajectory[:201]):
    vdos.processStep(tStep,ts.time)
    tStep += 1

vdos.postProcess()
vdos.outputGeometry("residueProperties.dat")
vdos.outputVACF("VACF.dat")
vdos.outputVDoS("VDoS.dat")

