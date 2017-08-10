from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
f = open('CD.csv', 'w')#file to write chain data in
print('Loading...')
pdb = PDBFile('CD.pdb')
forcefield = ForceField('amber99sb.xml', 'tip3p.xml')#create forcefield
modeller = Modeller(pdb.topology, pdb.positions)#initialize modeller
print('Adding hydrogens...')
modeller.addHydrogens(forcefield)
print('Adding solvent...')
modeller.addSolvent(forcefield, model='tip3p', padding=1*nanometer)#neutralize the system with ions
print('Minimizing...')
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.001*picoseconds)
platform = Platform.getPlatformByName('OpenCL')
simulation = Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('CDtest.pdb', 1000))
print('Saving...')

simulation.reporters.append(StateDataReporter(f, 100, step=True,kineticEnergy=True,
                                              potentialEnergy=True,totalEnergy = True,temperature=True,density = True))
#for every 100 steps, write into f 
simulation.step(30000)#run for 30000 steps
print('done')
