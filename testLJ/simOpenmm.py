#! /usr/bin/env python

from simtk.openmm.app import GmxTopMartini, GromacsGroFile
import simtk.openmm as mm
from simtk import unit
import numpy as np

def _setSystemParams():
    systemParams = {
                    'nonbondedCutoff': 1.2*unit.nanometers,
                    'nonbondedMethod': mm.NonbondedForce.CutoffPeriodic,
                   }

    return systemParams


def calcEnergy(gmxtopfile, gmxgrofile, trajfile):
    ''' '''
    temperature = 320 * unit.kelvin
    fricCoeff = 1.0 / unit.picoseconds
    stepsize = 0.01 * unit.picoseconds
    platform = mm.Platform.getPlatformByName('Reference')
    properties = {}
    traj = mm.app.PDBFile(trajfile)
    gro = GromacsGroFile(gmxgrofile)
    top = GmxTopMartini(gmxtopfile,
                        unitCellDimensions=gro.getUnitCellDimensions())
    sysParams = _setSystemParams()
    system = top.createSystem(**sysParams)
    integrator = mm.LangevinIntegrator(temperature, fricCoeff, stepsize)
    simulation = mm.app.Simulation(top.topology, system, integrator, platform,
                                   properties)
    observ=[]
    for frame in range(traj.getNumFrames()):
        simulation.context.setPositions(traj.getPositions(frame=frame))
        state = simulation.context.getState(getEnergy=True)
        observ.append( (frame *0.01,
                        state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)))

    observ =np.array(observ)
    np.savetxt('potential.txt',observ)



def main():
    datadir = "data/"
    simdir = "sim/"
    gmxtopfile = datadir + "xqCG.top"
    gmxgrofile = datadir + "xqCG.gro"
    trajfile = simdir + "bondedPout.pdb"
    n=mm.Platform.getNumPlatforms()
    for i in range(n):
        platform = mm.Platform.getPlatform(i)
        print i,":",platform.getName()

    calcEnergy(gmxtopfile, gmxgrofile, trajfile)
    print "DONE"
main()
