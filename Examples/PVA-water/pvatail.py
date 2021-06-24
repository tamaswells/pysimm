from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk

def tailer():
    s = system.read_pubchem_smiles('CCO')
    f = forcefield.Gaff2()
    s.apply_forcefield(f)
    
    c1 = s.particles[2]
    c1.linker = 'tail'
    c2 = s.particles[3]
    c2.linker = 'head'
         
    for b in c1.bonds:
        if b.a.elem == 'H' or b.b.elem == 'H':
            pb = b.a if b.b is c1 else b.b
            s.particles.remove(pb.tag, update=False)
            break
            
    s.remove_spare_bonding()
    
    lmps.quick_min(s, min_style='fire')
    
    s.add_particle_bonding()
    
    return s
    
def polymer_chain(length):
    mon = tailer()
    polym = random_walk(mon, length, forcefield=forcefield.Gaff2())
    return polym