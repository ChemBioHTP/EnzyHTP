
__doc__='''
This module make configuration file and deploy input file for different next-stage software. (Amber, Gaussian(QM/ONIOM))
Take all job set up that not relate to the PDB itself
-------------------------------------------------------------------------------------
Config --- Amber
        |- Gaussian
-------------------------------------------------------------------------------------
'''
import re

# >>>>>> Resource <<<<<<
# -----------------------------
# Cores available (used in MD(Amber) and QM(Gaussian) calculations)
#
n_cores = 24
# -----------------------------
# Per core memory in MB
#
max_core = 2000

class Config:
    # >>>>>> Resource <<<<<<
    # -----------------------------
    # Cores available (used in MD(Amber) and QM(Gaussian) calculations)
    #
    n_cores = n_cores
    # -----------------------------
    # Per core memory in MB
    #
    max_core = max_core
    # -----------------------------
    # debug info level
    # 0: No info. 1: normal debug (warning) 2: verbose debug (running log)
    debug = 1
    # -----------------------------
    # command line name for parallel computing
    # 
    PC_cmd = 'mpirun -np '+str(n_cores) 

    
    # >>>>>> Software <<<<<<
    class Amber:
        # -----------------------------
        # Amber Home dir
        #
        AmberHome = '$AMBERHOME'
        # -----------------------------
        # Default configuration for MD
        #
        # -----------------------------
        # Water box type and size
        # availible type: box/oct
        box_type = 'box'
        box_size = '10'
        # -----------------------------
        #            min
        # Minimize
        #  &cntrl
        #   imin  = 1,  ntx   = 1,  irest = 0,
        #   ntc   = 2,    ntf = 2,
        #   cut   = 10.0,
        #   maxcyc= 20000, ncyc  = {0.5maxcyc},
        #   ntpr  = {0.01maxcyc},   ntwx  = 0,
        #   ntr   = 1,	restraint_wt = 2.0, restraintmask = '@C,CA,N',
        #  /        
        conf_min={'ntc': '2', 'ntf': '2',
                  'cut': '10.0',
                  'maxcyc': 20000, 'ncyc': '0.5maxcyc',
                  'ntpr': '0.01maxcyc',
                  'ntr': '1', 'restraintmask': "'@C,CA,N'", 'restraint_wt': '2.0'
                  }
        # -----------------------------
        #            heat
        # Heat
        #  &cntrl
        #   imin  = 0,  ntx = 1, irest = 0,
        #   ntc   = 2, ntf = 2,
        #   cut   = 10.0,
        #   nstlim= 20000, dt= 0.002,
        #   tempi = 0.0,  temp0=300.0,  
        #   ntpr  = {0.01nstlim},  ntwx={nstlim},
        #   ntt   = 3, gamma_ln = 5.0,
        #   ntb   = 1,  ntp = 0,
        #   iwrap = 1,
        #   nmropt= 1,
        #   ig    = -1,
        #   ntr   = 1,	restraint_wt = 2.0, restraintmask = '@C,CA,N',
        #  /
        #  &wt
        #   type  = 'TEMP0',
        #   istep1= 0, istep2={0.9nstlim},
        #   value1= 0.0, value2=300.0,
        #  /
        #  &wt
        #   type  = 'TEMP0',
        #   istep1= {A_istep2+1}, istep2={nstlim},
        #   value1= 300.0, value2=300.0,
        #  /
        #  &wt
        #   type  = 'END',
        #  /
        conf_heat={'ntc': '2', 'ntf': '2',
                   'cut': '10.0',
                   'nstlim': 20000, 'dt': '0.002',
                   'tempi': '0.0', 'temp0': '300.0',
                   'ntpr': '0.01nstlim', 'ntwx': 'nstlim',
                   'ntt': '3', 'gamma_ln': '5.0',
                   'iwarp': '1',
                   'ntr': '1', 'restraintmask': "'@C,CA,N'", 'restraint_wt': '2.0',
                   'A_istep2':'0.9nstlim', 'B_istep1':'A_istep2+1'
                  }
        # -----------------------------
        #            equi
        # Equilibration: constant pressure
        #  &cntrl
        #   imin  = 0,  ntx = 5,  irest = 1,
        #   ntf   = 2,  ntc = 2,
        #   nstlim= 500000, dt= 0.002,
        #   cut   = 10.0,
        #   temp0 = 300.0,
        #   ntpr  = {0.002nstlim}, ntwx = 5000,
        #   ntt   = 3, gamma_ln = 5.0,
        #   ntb   = 2,  ntp = 1,
        #   iwrap = 1,
        #   ig    = -1,
        #   ntr   = 1,	restraint_wt = 2.0,
        #   restraintmask = '@C,CA,N',
        #  /
        conf_equi={'ntx': '5', 'irest': '1',
                   'ntc': '2', 'ntf': '2',
                   'cut': '10.0',
                   'nstlim': 500000, 'dt': '0.002',
                   'temp0': '300.0',
                   'ntpr': '0.002nstlim', 'ntwx': '5000', # default 10ps
                   'ntt': '3', 'gamma_ln': '5.0',
                   'iwarp': '1',
                   'ntr': '1', 'restraintmask': "'@C,CA,N'", 'restraint_wt': '2.0', # the later two are only used when ntr = 1
                  }
        # -----------------------------
        #            prod
        # Production: constant pressure
        #  &cntrl
        #   imin  = 0, ntx = 1, irest = 0,
        #   ntf   = 2,  ntc = 2,
        #   nstlim= 50000000, dt= 0.002,
        #   cut   = 10.0,
        #   temp0 = 300.0,
        #   ntpr  = 50000, ntwx = 5000,
        #   ntt   = 3, gamma_ln = 5.0,
        #   ntb   = 2,  ntp = 1,
        #   iwrap = 1,
        #   ig    = -1,
        #  /
        conf_prod={'ntx': '5', 'irest': '1',
                   'ntc': '2', 'ntf': '2',
                   'cut': '10.0',
                   'nstlim': 50000000, 'dt': '0.002',
                   'temp0': '300.0',
                   'ntpr': '0.001nstlim', 'ntwx': '5000', # default 10ps
                   'ntt': '3', 'gamma_ln': '5.0',
                   'iwarp': '1',
                   'ntr': '0', 'restraintmask': None, 'restraint_wt': '2.0', # the later two are only used when ntr = 1
                  }


    class Gaussian:
        # -----------------------------
        #   >>>>>>>>ONIOM<<<<<<<<
        # -----------------------------

        # -----------------------------
        # Cores for gaussian job (higher pirority)
        # 
        n_cores = n_cores
        # -----------------------------
        # Per core memory in MB for gaussian job (higher pirority)
        # 
        max_core = max_core
        # -----------------------------
        # key words for ONIOM job (strategy: use work type as a key to indicate different keyword set.)
        # Use mechanic embedding as default / hardfirst means find user assigned parameters first / iops are for reduce the size of output. 
        #
        # for oniom method level
        # -- open for edit --
        om_lvl = ['wb97xd/6-31g(d)', 'amber=hardfirst'] # add in the order of theory level  e.g.: ['wb97xd/def2svp', 'PM7', 'amber'] 
        # -- open for edit --
        oniom_kw = 'oniom(' + ':'.join(om_lvl) + ')'

        # complete keywords
        keywords = {'spe'  :[ oniom_kw, 'nosymm', 'geom=connectivity', 
                             'iop(2/9=1000,3/33=0,3/36=-1,3/160=2,4/33=0,5/33=0,6/7=2,6/9=2,6/10=2,6/12=2,7/8=10,2/15=3)'
                            ],
                    'opt'  :[ oniom_kw, 'opt=(modredundant,calcfc)', 'freq=(selectnormalmodes,hpmodes)', 'nosymm', 'geom=connectivity',
                             'iop(2/9=1000,3/33=0,3/36=-1,3/160=2,4/33=0,5/33=0,6/7=2,6/9=2,6/10=2,6/12=2,7/8=10,2/15=3)'
                            ],
                    'tsopt':[ oniom_kw, 'opt=(calcfc,ts,noeigen,gdiis)', 'freq=(selectnormalmodes,hpmodes)', 'nosymm', 'geom=connectivity', 
                             'iop(2/9=1000,3/33=0,3/36=-1,3/160=2,4/33=0,5/33=0,6/7=2,6/9=2,6/10=2,6/12=2,7/8=10,2/15=3)'
                            ]
                   }
        # -----------------------------
        # layer settings for oniom
        # preset: 0: no preset (fill layer atom manually) 1: preset_1 -> xxx
        layer_preset = 0
        layer_atoms = []

class Layer:
	'''
	set / use preset of oniom layer
	------------
	PDB: related PDB object
	layer_atoms: list of layer's list of atom indexes 
	'''
	def __init__(self, PDB_obj, atom_lists, if_set=0):
		'''
		general way to assign layer: list of atom indexes
		'''
		self.PDB = PDB_obj
		if not if_set:
			self.layer_atoms = []
			for layer in atom_lists:
				atoms = []
				for atom_str in layer.split(','):
					if '-' in atom_str:
						if re.match(r'[0-9]+-[0-9]+',atom_str) == None:
							raise Exception('Wrong layer syntax of atom indexes. e.g.: 1,2,3-5')
						ID = atom_str.split('-')
						atoms.extend(range(int(ID[0]),int(ID[1])+1))
					else:
						atoms.append(int(atom_str))
				self.layer_atoms.append(atoms)
		else:
			#for preset
			self.layer_atoms = atom_lists

		if Config.debug >= 2:
			print('current layer atoms: ' + repr(self.layer_atoms))
									
	
	@classmethod
	def preset(cls, PDB_obj, set_id):
		'''
		preset layer settings for ONIOM
		set_id = 
		1: 只包含底物
		2: 包含底物和“关键残基”
		3: 包含一定半径内的所有残基、配体
		4: 根据某种参数选择QM区域
		'''
		layer_atoms = []
		if set_id == 1:
			想想怎么设置预设
		
		return cls(PDB_obj, layer_atoms)
