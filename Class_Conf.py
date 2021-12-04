
__doc__='''
This module make configuration file and deploy input file for different next-stage software. (Amber, Gaussian(QM/ONIOM))
Take all job set up that not relate to the PDB itself
-------------------------------------------------------------------------------------
Config --- Amber
        |- Gaussian
        |- Multiwfn
-------------------------------------------------------------------------------------
'''
from math import exp
import re
from helper import mkdir, line_feed


class Config:
    # >>>>>> Resource <<<<<<
	# -----------------------------
    # Cores available (used in MD(Amber) and QM(Gaussian) calculations)
    #
    n_cores = 24
    # -----------------------------
    # Per core memory in MB
    #
    max_core = 2000
    # -----------------------------
    # debug info level
    # 0: No info. 1: normal debug (warning) 2: verbose debug (running log)
    debug = 1
    # -----------------------------
    # command line name for cpu parallel computing
    # 
    PC_cmd = 'mpirun'
    @classmethod
    def get_PC_cmd(cls):
        if cls.PC_cmd == 'mpirun':
            return 'mpirun -np '+str(cls.n_cores) 
        return cls.PC_cmd

    
    # >>>>>> Software <<<<<<
    class Amber:
        # -----------------------------
        # Amber Home dir
        #
        AmberHome = '$AMBERHOME'
        # -----------------------------
        # User defined Amber CPU path (default: $AMBERHOME/bin/sander.MPI)
        #
        AmberEXE_CPU = None
        # -----------------------------
        # User defined Amber GPU path (default: $AMBERHOME/bin/pmemd.cuda)
        #
        AmberEXE_GPU = None
        # -----------------------------
        # Default configuration for MD
        #
        # -----------------------------
        # Water box type and size
        # availible type: box/oct
        box_type = 'oct'
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
        
        #==============================
        # Method to express Amber EXE path and PC_cmd
        #
        @classmethod
        def get_Amber_engine(cls, engine='Amber_GPU'):
            '''
            Give default value to Amber cpu/gpu engines
            engine: ['Amber_CPU', 'Amber_GPU'] 
            ---
            return (PC_cmd, AmberEXE path)
            '''
            #san check
            if engine not in ['Amber_CPU', 'Amber_GPU']:
                raise Exception("get_Amber_engin only support ['Amber_CPU', 'Amber_GPU']")

            if engine == 'Amber_CPU':
                if cls.AmberEXE_CPU == None:
                    engine_path = cls.AmberHome+'/bin/sander.MPI'
                else:
                    engine_path = cls.AmberEXE_CPU
                PC_cmd = Config.get_PC_cmd()

            if engine == 'Amber_GPU':
                if cls.AmberEXE_GPU == None:
                    engine_path = cls.AmberHome+'/bin/pmemd.cuda'
                else:
                    engine_path = cls.AmberEXE_GPU
                PC_cmd = ''

            return PC_cmd, engine_path



        class MMPBSA:
            # -----------------------------
            # User defined MMPBSA.py(.MPI) exe dir other than $AMBERHOME/bin/MMPBSA.py.MPI
            #
            MMPBSA_EXE = None

            #==============================
            # Method to express MMPBSA EXE path
            #
            @classmethod
            def get_MMPBSA_engine(cls):
                '''
                Give default value to MMPBSA engine
                Only support MPI version now.
                ---
                return engine_path
                '''
                if cls.MMPBSA_EXE == None:
                    engine_path = Config.Amber.AmberHome+'/bin/MMPBSA.py.MPI'
                else:
                    engine_path = cls.MMPBSA_EXE
                
                return engine_path
            
            # -----------------------------
            #           MMPBSA.in
            #
            # GB and PB calculation
            # &general
            #   startframe=1, interval=100,
            #   verbose=1, keep_files=0,
            # /
            # &gb
            #   igb=5, saltcon=0.150,
            # /
            # &pb
            #   istrng=0.15, fillratio=4.0
            # /
            conf_in = {
                'startframe': 1, 'endframe': None, 'interval': 1,
                'verbose': 1, 'keep_files': 0,
                'if_gb' : 1,
                'igb' : 5, 'saltcon' : 0.15,
                'if_pb' : 1,
                'istrng': 0.15, 'fillratio': 4.0
            }

            @classmethod
            def build_MMPBSA_in(cls, out_path=''):
                '''
                build MMPBSA.in in out_path
                '''
                if out_path == '':
                    mkdir('./tmp')
                    out_path = './tmp/MMPBSA.in'

                # make lines
                frame_line = '  '
                for i in ('startframe', 'endframe', 'interval'):
                    if cls.conf_in[i] != None:
                        frame_line = frame_line + i + '=' + str(cls.conf_in[i]) + ', '
                output_line = '  verbose='+ str(cls.conf_in['verbose']) +', keep_files='+ str(cls.conf_in['keep_files']) +','
                gb_line = '  igb='+str(cls.conf_in['igb']) + ', saltcon='+str(cls.conf_in['saltcon'])+','
                pb_line = '  istrng='+str(cls.conf_in['istrng'])+', fillratio='+str(cls.conf_in['fillratio'])

                with open(out_path, 'w') as of:
                    print('GB and PB calculation' , end=line_feed, file=of)
                    print('&general' , end=line_feed, file=of)
                    print(frame_line  , end=line_feed, file=of)
                    print(output_line  , end=line_feed, file=of)
                    print('/' , end=line_feed, file=of)
                    print('&gb' , end=line_feed, file=of)
                    print(gb_line  , end=line_feed, file=of)
                    print('/' , end=line_feed, file=of)
                    print('&pb' , end=line_feed, file=of)
                    print(pb_line  , end=line_feed, file=of)
                    print('/' , end=line_feed, file=of)





    class Gaussian:
        # -----------------------------
        #   >>>>>>>>ONIOM<<<<<<<<
        # -----------------------------

        # -----------------------------
        # Cores for gaussian job (higher pirority)
        # 
        # n_cores = n_cores
        # -----------------------------
        # Per core memory in MB for gaussian job (higher pirority)
        # 
        # max_core = max_core
        # -----------------------------
        # Executable g16 command for current environment
        # 
        g16_exe = 'g16'
        # -----------------------------
        # Executable g16 command for current environment
        # 
        g09_exe = 'g09'
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


    class Multiwfn:
        # -----------------------------
        # Cores for Multiwfn job (higher pirority)
        # 
        # n_cores = n_cores
        # -----------------------------
        # Per core memory in MB for Multiwfn job (higher pirority)
        # 
        # max_core = max_core
        # -----------------------------
        # Executable Multiwfn command for current environment
        # 
        exe = 'Multiwfn'
        # -----------------------------
        # Path of Multiwfn folder
        # 
        DIR = '$Multiwfnpath'
        

class Layer:
	'''
	set / use preset of oniom layer
	------------
	PDB: related PDB object
	layer: list of layer's list of atom indexes 
	'''
	def __init__(self, PDB_obj, atom_lists, if_set=0):
		'''
		general way to assign layer: list of atom indexes
		'''
		self.PDB = PDB_obj
		if not if_set:
			self.layer = []

			if len(atom_lists) not in [2,3]:
				raise Exception('Layer.__init__: only support 2 or 3 layers. Input layers: '+str(len(atom_lists)))
				
			for layer in atom_lists:
				atoms = []
				for atom_str in layer.split(','):
					if '-' in atom_str:
						if re.match(r'[0-9]+-(?:[0-9]+|L)',atom_str) == None:
							raise Exception('Wrong layer syntax of atom indexes. e.g.: 1,2,3-5')
						ID = atom_str.split('-')
						for i, j in enumerate(ID):
							if j == 'L':
								ID[i] = PDB_obj.get_last_A_id()
						atoms.extend(range(int(ID[0]),int(ID[1])+1))
					else:
						atoms.append(int(atom_str))
				self.layer.append(atoms)
		else:
			#for preset
			self.layer = atom_lists

		if Config.debug >= 1:
			print('current layer atoms: ')
			for i,atoms in enumerate(self.layer):
				print(str(i)+':',len(atoms))
		if Config.debug >= 2:
			print('current layer atoms: ' + repr(self.layer))
									
	
	@classmethod
	def preset(cls, PDB_obj, set_id, lig_list=[]):
		'''
		preset layer settings for ONIOM
                ------------
		set_id = 
			(two layers)
			1: Substrate only
			2: Key ligands (specify ligand index or name; all ligand by default)
			3: Substrate and key residues (manually assigned)
			4: Substrate and all residues/ligands within a assigned radius (Need to keep consistant molecular number)
			5: Based on some parameters to select the QM region
		lig_list: (set_id = 2) key ligand index
		'''
		layer_atoms = []
		PDB_obj.get_stru()
		stru = PDB_obj.stru

		if set_id not in [1,2,3,4,5]:
			raise Exception('Only support 1-5 set_id now. You are entering: ' +str(set_id))

		if set_id == 1:
			# TODO: connect with the database and recognize the substrate in the future
			pass

		if set_id == 2:
			if lig_list == []:
				h_atoms=[]
				lig_names = []
				for lig in stru.ligands:
					for atom in lig:
						h_atoms.append(atom.id)
					lig_names.append(repr((lig.name, lig.id)))
				layer_atoms.append(h_atoms)
				# log
				if Config.debug >= 1:
					print('Layer.preset(set_id=2): No id in lig_list, set all ligands to high layer by default')
					print(' '.join(lig_names))
			else:
				h_atoms=[]
				for lig_inp in lig_list:
					if type(lig_inp) == int:
						for lig in stru.ligands:
							if lig.id == lig_inp:
								for atom in lig:
									h_atoms.append(atom.id)
					if type(lig_inp) == str:
						for lig in stru.ligands:
							if lig.name == lig_inp:
								for atom in lig:
									h_atoms.append(atom.id)

			l_atoms = list(set(stru.get_atom_id()).difference(set(h_atoms)))

			#debug
			# a = stru.get_atom_id()
			# b = set(a)
			# for i in b:
			# 	count = 0
			# 	for j in a:
			# 		if j == i:
			# 			count +=1
			# 	if count > 1:
			# 		print(i, ':', count)
			
			layer_atoms = [h_atoms,l_atoms]
		
		if set_id == 3:
			pass
		
		return cls(PDB_obj, layer_atoms, if_set=1)


	'''
	====
	Special Method
	====
	'''

	def __getitem__(self, key: int):
		'''
		Chain_obj[int]: Chain_obj.residues[int]
		-----
		use residue index within the chain, start from 0
		'''
		return self.layer[key]


	def __len__(self):
		return len(self.layer)