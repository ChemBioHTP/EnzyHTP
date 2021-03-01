
__doc__='''
This module make configuration file and deploy input file for different next-stage software. (Amber, Gaussian(QM/ONIOM))
Take all job set up that not relate to the PDB itself
-------------------------------------------------------------------------------------
Class Conf
-------------------------------------------------------------------------------------
add_tag(self,job_tag='')
    Add a custom tag for filenames including: 
        all .in files;
-------------------------------------------------------------------------------------
Amber Classical MD Simulation
-------------------------------------------------------------------------------------
set_MD_min(self,maxcyc='20000',ncyc='10000',ntpr='1000',cut='10.0',ntc='1')
set_MD_heat(self,nstlim='20000',dt='0.001',tempi='0.0',temp0='300.0',A_istep2=nstlim*0.9,ntpr='100',ntwx='20000',cut='10.0',ntc='2',ntf='2')
set_MD_equi(self,nstlim='5000',dt='0.001',cut='10.0',temp0='300.0',ntpr='5000',ntwx='5000',ntc='2',ntf='2')
set_MD_prod(self,nstlim='100000',dt='0.001',cut='10.0',temp0='300.0',ntpr='5000',ntwx=nstlim/10000,ntc='2',ntf='2')
-------------------------------------------------------------------------------------
Gaussian ONIOM calculation
-------------------------------------------------------------------------------------

-------------------------------------------------------------------------------------
File Deployment
-------------------------------------------------------------------------------------
deploy(path='./',clean=True)       #deploy all files set in the deploy_list
                                   #clean the deploy_list after deployment
-------------------------------------------------------------------------------------
'''

class Config:
    # >>>>>> Resource <<<<<<
    # -----------------------------
    # Cores available (used in MD(Amber) and QM(Gaussian) calculations)
    #
    n_cores = 8
    # -----------------------------
    # Per core memory in MB
    #
    max_core = 2000
    # -----------------------------
    # debug info level
    # 0: No info. 1: normal debug 2: verbose debug
    debug = 1

    
    # >>>>>> Software <<<<<<
    class Amber:
        # -----------------------------
        # Amber Home dir
        #
        AmberHome = '$AMBERHOME'


class Conf():

    # A custom tag for filenames including: all .in files; 
    tag=''

    # a list for deployment
    deploy_list=[]
    deploy_path=''

    #Amber simulation
    min_conf=''
    min_fn=''
    heat_conf=''
    heat_fn=''
    equi_conf=''
    equi_fn=''
    prod_conf=''
    prod_fn=''

    def add_tag(self,job_tag=''):
        '''
        Add a custom tag for filenames including: 
            all .in files;
        '''

        if job_tag != '':
            self.tag='_'+job_tag



    def set_MD_min(self,maxcyc='20000',ncyc='10000',ntpr='1000',cut='10.0',ntc='2',ntf='2'):
        '''
        Set configuration for a minimization job
        default value: (self,maxcyc='20000',ncyc='10000',ntpr='1000',cut='10.0',ntc='2',ntf='2')
        '''

        # WARNING: support SHAKE or ntr restrain in the future
        self.min_conf='''Minimize
 &cntrl
  imin  = 1,  ntx   = 1,  irest = 0,
  ntc   = '''+ntc+''',    ntf = '''+ntf+''',
  cut   = '''+cut+''',
  maxcyc= '''+maxcyc+''', ncyc  = '''+ncyc+''',
  ntpr  = '''+ntpr+''',   ntwx  = 0,
  ntr   = 1,	restraint_wt = 2.0,
  restraintmask = '@C,CA,N',
 /
 '''

        self.min_fn='min'+self.tag+'.in'
        self.deploy_list=self.deploy_list+[(self.min_fn,self.min_conf),]


    def set_MD_heat(self,nstlim='20000',dt='0.001',tempi='0.0',temp0='300.0',A_istep2='',ntpr='100',ntwx='20000',cut='10.0',ntc='2',ntf='2'):
        '''
        Set configuration for a heat job
        default value: (self,nstlim='20000',dt='0.001',tempi='0.0',temp0='300.0',A_istep2=nstlim*0.9,ntpr='100',ntwx='20000',cut='10.0',ntc='2',ntf='2')
        '''

        if A_istep2 == '':
            A_istep2=str(int(float(nstlim)*0.9))
        B_istep1=str(int(A_istep2)+1)

        self.heat_conf='''Heat
 &cntrl
  imin  = 0,  ntx = 1, irest = 0,
  ntc   = '''+ntc+''', ntf = '''+ntf+''',
  cut   = '''+cut+''',
  nstlim= '''+nstlim+''', dt= '''+dt+''',
  tempi = '''+tempi+''',  temp0='''+temp0+''',  
  ntpr  = '''+ntpr+''',  ntwx='''+ntwx+''',
  ntt   = 3, gamma_ln = 5.0,
  ntb   = 1,  ntp = 0,
  iwrap = 1,
  nmropt= 1,
  ig    = -1,
  ntr   = 1,	restraint_wt = 2.0,
  restraintmask = '@C,CA,N',
 /
 &wt
  type  = 'TEMP0',
  istep1= 0, istep2='''+A_istep2+''',
  value1= '''+tempi+''', value2='''+temp0+''',
 /
 &wt
  type  = 'TEMP0',
  istep1= '''+B_istep1+''', istep2='''+nstlim+''',
  value1= '''+temp0+''', value2='''+temp0+''',
 /
 &wt
  type  = 'END',
 /
 '''

        self.heat_fn='heat'+self.tag+'.in'
        self.deploy_list=self.deploy_list+[(self.heat_fn,self.heat_conf),]


    def set_MD_equi(self,nstlim='500000',dt='0.002',cut='10.0',temp0='300.0',ntpr='10000',ntwx='10000',ntc='2',ntf='2'):
        '''
        Set configuration for a equilibration job
        default value: (self,nstlim='500000',dt='0.002',cut='10.0',temp0='300.0',ntpr='5000',ntwx='5000',ntc='2',ntf='2')
        '''

        self.equi_conf='''Equilibration:constant pressure
 &cntrl
  imin  = 0,  ntx = 1,  irest = 1,
  ntf   = '''+ntf+''',  ntc = '''+ntc+''',
  nstlim= '''+nstlim+''', dt= '''+dt+''',
  cut   = '''+cut+''',
  temp0 = '''+temp0+''',
  ntpr  = '''+ntpr+''', ntwx = '''+ntwx+''',
  ntt   = 3, gamma_ln = 5.0,
  ntb   = 2,  ntp = 1,
  iwrap = 1,
  ig    = -1,
  ntr   = 1,	restraint_wt = 2.0,
  restraintmask = '@C,CA,N',
 /
 '''

        self.equi_fn='equi'+self.tag+'.in'
        self.deploy_list=self.deploy_list+[(self.equi_fn,self.equi_conf),]

    def set_MD_prod(self,nstlim='50000000',dt='0.002',cut='10.0',temp0='300.0',ntpr='25000',ntwx='',ntc='2',ntf='2',irest='1'):
        '''
        Set configuration for a production job
        default value: (self,nstlim='50000000',dt='0.002',cut='10.0',temp0='300.0',ntpr='5000',ntwx=nstlim/10000,ntc='2',ntf='2')
        '''

        #default 10k frames
        if ntwx == '':
            ntwx=str(int(float(nstlim)/10000))

        self.prod_conf='''Equilibration:constant pressure
 &cntrl
  imin  = 0, ntx = 1, irest = '''+irest+''',
  ntf   = '''+ntf+''',  ntc = '''+ntc+''',
  nstlim= '''+nstlim+''', dt= '''+dt+''',
  cut   = '''+cut+''',
  temp0 = '''+temp0+''',
  ntpr  = '''+ntpr+''', ntwx = '''+ntwx+''',
  ntt   = 3, gamma_ln = 5.0,
  ntb   = 2,  ntp = 1,
  iwrap = 1,
  ig    = -1,
 /
 '''
        
        self.prod_fn='prod'+self.tag+'.in'
        self.deploy_list=self.deploy_list+[(self.prod_fn,self.prod_conf),]


    def deploy(self,path='.',clean=True):

        self.deploy_path=path
        for filename,fileStr in self.deploy_list:
            conf_file=open(self.deploy_path+'/'+filename,'w')
            conf_file.write(fileStr)
            conf_file.close()
        
        if clean:
            self.deploy_list=[]

#TestOnly
# a=Conf()
# a.set_MD_min()
# a.set_MD_heat()
# a.set_MD_equi()
# a.set_MD_prod()
# a.deploy()


