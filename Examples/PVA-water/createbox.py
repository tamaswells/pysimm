from pysimm import system, lmps, forcefield
from pysimm.apps.random_walk import random_walk
from pva import monomer
from pvahead import header
from pvatail import tailer

def run(test=False):
    # we'll create a pe monomer from the pysimm.models database
    pe = monomer()
    # add head and tail
    head = header()
    tail = tailer()  
    # we'll instantiate a GAFF2 forcefield object for use later
    f = forcefield.Gaff2()
    
    # the monomers do not have any charges, so we will derive partial charges using the gasteiger algorithm
    pe.apply_charges(f, charges='gasteiger')
    head.apply_charges(f, charges='gasteiger')
    tail.apply_charges(f, charges='gasteiger')
    
    # then we use the pattern argument to define how many of each "monomers" to add. Let's add 5 more monomers to our chain
    
    # first import the copolymer function
    
    from pysimm.apps.random_walk import copolymer
    
    polymer = copolymer([head, pe, tail], 16, pattern=[1, 14, 1], forcefield=f)
    polymer.write_xyz('polymer.xyz')
    polymer.write_lammps('polymer.lmps')
        
    tip3p = forcefield.tip3p.Tip3p()
    water = system.read_pubchem_smiles('O')
    water.apply_forcefield(tip3p)
    molecule_list=[polymer, water]

    n_molecules=[2, 192]

    s=system.replicate(molecule_list, n_molecules , density=0.5)

    min_settings = {
        'name': 'fire_min',
        'min_style': 'fire',
        'print_to_screen': True
    }

    nvt_settings = {
        'name': 'nvt_md',
        'print_to_screen': True,
        'ensemble': 'nvt',
        'temperature': {
            'start': 300,
            'stop': 300
        },
        'new_v': True,
        'length': 10000
    }

    lmps.quick_min(s, **min_settings)
    lmps.quick_md(s, **nvt_settings)

    s.write_xyz('mixture.xyz')
    s.write_lammps('mixture.lmps')

if __name__ == '__main__':
    run()
