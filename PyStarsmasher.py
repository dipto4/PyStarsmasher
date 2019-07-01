import numpy as np
import os
from optparse import OptionParser
from astropy import units
from errors import *
from ctypes import cdll, c_int, c_double, byref, c_bool
from mpi4py import MPI
from Star import Star
from readNbody import *
import fnmatch
#from multiprocessing import Process

class Starsmasher(object):
    __instance__ = None

    def __init__(self,num_of_workers=1):
        Starsmasher.__instance__ = self

        #self.library = cdll.LoadLibrary('./libstarsmasher.so')
        self.toolsLibrary = cdll.LoadLibrary('./tools/tools.so')
        # The following variables are present in
        # sph.input

        self.num_of_workers = num_of_workers


        self.ndisplace = 0
        self.displacex = 0.0
        self.displacey = 0.0
        self.displacez = 0.0

        self.semimajoraxis = 0.0
        self.bimpact = -1.e30
        self.e0 = -1.e30
        self.vinf2 = 1.e30

        self.tf = 50000
        self.dtout = 100

        self.n = 100000

        self.gflag = 1
        self.nnopt = 22 + self.gflag
        self.nav = 3
        self.alpha = 1
        self.beta = 2
        self.ngr = 3
        self.hco = 1.0
        self.hfloor = 0.0
        self.nrelax = 1
        self.trelax =1.e30

        self.sep0 = 200
        self.equalmass=0
        self.treloff=0
        self.tresplintmuoff=0.0

        self.nitpot=1
        self.tscanon = 0
        self.sepfinal = 1.e30
        self.nintvar = 2

        self.ngravprocs=-2
        self.qthreads = 0

        self.mbh = 20.0
        self.runit = 6.9599e10
        self.munit = 1.9891e33

        self.cn1 = 0.5e0
        self.cn2 = 0.06e0
        self.cn3 = 0.06e0
        self.cn4 = 1.e30
        self.cn5 = 0.02e0
        self.cn6 = 0.02e0
        self.cn7 = 4.0e0

        self.computeexclusivemode = 0

        self.omega_spin=0.0e0
        #check if this is required
        self.ppn = 12
        self.neos = 1
        self.nselfgravity = 1
        self.gam = 5.0/3.0
        self.reat = -1.0e0
        self.npoly = 1.50

        #starmass and starradius are in units of
        #Msun and Rsun

        self.starmass = 1.0e0
        self.starradius = 1.0e0
        self.ncooling = 0
        self.nkernel = 0
        self.teq = 100.0e0

        self.tjumpahead=1.0e30
        self.startfile1='sph.start1u'
        self.startfile2='sph.start2u'

        self.eosfile='sph.eos'
        self.opacityfile='sph.opacity'

        self.profilefile='eg.last1.muse_s2mm'
        self.throwaway = False

        self.stellarevolutioncodetype=1

        # The following parameters are present in sph.init

        self.simulationtype = '1es'
        self.dirname = 'data'

        self.star1 = Star()
        self.star2 = Star()


    def __errorHandler(self):
        # 1 = Error in Inputs
        # 2 = Internal Error. Check log / inputs
        # 2 is triggered when control returns to python and tfinal
        # in starsmasher is less than tfinal that the user entered

        #Raise Error if simulation type mismatch
        try:
            simulationtypes = {"1es" : 1,
                    "1mc" : 2,
                    "2cr" : 3,
                    "bps" : 4,
                    "bph" : 5,
                    "hyp" : 6,
                    "hbs" : 7,
                    "erg" : 8,
                    "tri" : 9,
                    "bhe" : 10,
                    "rin" : 11,
                    "grs" : 12,
                    "txt" : 13,
                    "dbl" : 14}

            value = simulationtypes.get(self.simulationtype, None)

            print value
            if (value is None):
                raise SimulationTypeError

            if(self.n <= 0):
                print("n must be a positive number")
                raise BadValueError
            if(self.starmass <= 0):
                print("starmass must be a positive number")
                raise BadValueError
            if(self.starradius <= 0):
                print("starradius must be a positive number")
                raise BadValueError
            if(self.npoly <= 0):
                print("npoly has to be greater than zero")
                raise BadValueError
            if(self.num_of_workers <=0 or isinstance(self.num_of_workers,float)):
                print("Something wrong with num_of_workers")
                raise BadValueError

            if(self.tf <= 0 ):
                print("tf cannot be less than 0")
                raise BadValueError
        except SimulationTypeError:
            print("Unexpected Simulation Type")
        except BadValueError:
            print("Unexpected value found!")
        pass
    #TODO: implement the mpi based setValuesStarsmasher_mpi function
    #TODO: throwaway causes issues. CHange it

    def __setupValuesStarsmasher_mpi(self):

        self.__getGamForStar()

        dtype_input = [('ndisplace',np.int32),('displacex',np.float64),('displacey',np.float64),
                ('displacez',np.float64),('semimajoraxis',np.float64),('e0',np.float64),
                ('bimpact',np.float64),('trelax',np.float64),('vinf2',np.float64),('tf',np.float64),
                ('dtout',np.float64),('equalmass',np.float64),('n',np.int32),('gflag',np.int32),
                ('nnopt',np.int32),('nav',np.int32),('ngr',np.int32),('nrelax',np.int32),
                ('alpha',np.float64),('beta',np.float64),('cn1',np.float64),('cn2',np.float64),
                ('cn3',np.float64),('cn4',np.float64),('cn5',np.float64),('cn6',np.float64),
                ('cn7',np.float64),('hco',np.float64),('hfloor',np.float64),('tscanon',np.float64),
                ('sepfinal',np.float64),('sep0',np.float64),('treloff',np.float64),('tresplintmuoff',np.float64),
                ('nitpot',np.int32),('nintvar',np.int32),('ngravprocs',np.int32),('qthreads',np.int32),
                ('mbh',np.float64),('runit',np.float64),('munit',np.float64),
                ('computeexclusivemode',np.int32),('ppn',np.int32),('omega_spin',np.float64),
                ('neos',np.int32),('nselfgravity',np.int32),('ncooling',np.int32),('nkernel',np.int32),
                ('gam',np.float64),('reat',np.float64),('teq',np.float64),('tjumpahead',np.float64),
                ('starmass',np.float64),('starradius',np.float64),
                ('stellarevolutioncodetype',np.int32),('npoly',np.float64)]


        dtype_input_aligned = np.dtype(dtype_input,align=True)


        inputs = np.array([(self.ndisplace,self.displacex,self.displacey,self.displacez,
            self.semimajoraxis,self.e0,self.bimpact,self.trelax,self.vinf2,self.tf,
            self.dtout,self.equalmass,self.n,self.gflag,self.nnopt,self.nav,self.ngr,
            self.nrelax,self.alpha,self.beta,self.cn1,self.cn2,self.cn3,self.cn4,self.cn5,
            self.cn6,self.cn7,self.hco,self.hfloor,self.tscanon,self.sepfinal,self.sep0,
            self.treloff,self.tresplintmuoff,self.nitpot,self.nintvar,self.ngravprocs,
            self.qthreads,self.mbh,self.runit,self.munit,self.computeexclusivemode,
            self.ppn,self.omega_spin,self.neos,self.nselfgravity,self.ncooling,
            self.nkernel,self.gam,self.reat,self.teq,self.tjumpahead,self.starmass,
            self.starradius,self.stellarevolutioncodetype,self.npoly)],dtype=dtype_input_aligned)

        dtype_input_double = None

        inputs_double_s1 = None
        inputs_double_s2 = None

        if(self.simulationtype == 'dbl'):
            dtype_input_double = [('x',np.float64),('y',np.float64),('z',np.float64),
                    ('vx',np.float64),('vy',np.float64),('vz',np.float64),('m',np.float64)]


            inputs_double_s1 = np.array([(self.star1.x,self.star1.y,self.star1.z,
                self.star1.vx,self.star1.vy,self.star1.vz,self.star1.m)],dtype=dtype_input_double)
            inputs_double_s2 = np.array([(self.star2.x,self.star2.y,self.star2.z,
                self.star2.vx,self.star2.vy,self.star2.vz,self.star2.m)],dtype=dtype_input_double)



        if(self.simulationtype == 'dbl'):
            return inputs, dtype_input, inputs_double_s1, inputs_double_s2, dtype_input_double
        else:
            return inputs, dtype_input

    def __setValuesStarsmasher(self):

        if(self.simulationtype == 'dbl'):
            inputs, dtype_input, inputs_double_s1, inputs_double_s2, dtype_input_double = self.__setupValuesStarsmasher_mpi()
        else:
            inputs, dtype_input = self.__setupValuesStarsmasher_mpi()

        input_field_name = [lis[0] for lis in dtype_input]
        input_field_type = [lis[1] for lis in dtype_input]

        input_struct_size = inputs.nbytes

        print input_struct_size



        mpitype_dict = {np.int32:MPI.INT, np.float64:MPI.DOUBLE}

        input_offsets = [inputs.dtype.fields[field][1] for field in input_field_name]
        input_field_mpitypes = [mpitype_dict[dtype] for dtype in input_field_type]

        #print input_field_mpitypes
        #print len(input_offsets)



        input_structtype = MPI.Datatype.Create_struct([1]*len(input_field_name),input_offsets,input_field_mpitypes)
        input_structtype = input_structtype.Create_resized(0,input_struct_size)
        input_structtype.Commit()

        input_double_structtype = None


        if(self.simulationtype == 'dbl'):
            input_double_field_name = [lis[0] for lis in dtype_input_double]
            input_double_field_type = [lis[1] for lis in dtype_input_double]

            input_double_struct_size = inputs_double_s1.nbytes

            input_double_offsets = [inputs_double_s1.dtype.fields[field][1] for field in input_double_field_name]
            input_double_field_mpitypes = [mpitype_dict[dtype] for dtype in input_double_field_type]

            input_double_structtype = MPI.Datatype.Create_struct([1]*len(input_double_field_name),input_double_offsets,input_double_field_mpitypes)

            input_double_structtype = input_double_structtype.Create_resized(0,input_double_struct_size)
            input_double_structtype.Commit()


            return inputs, input_structtype, inputs_double_s1, inputs_double_s2, input_double_structtype
        else:
            return inputs, input_structtype



        pass

    def __getErrorFromStarsmasher():
        pass

    def __getTimeFromStarsmasher():
        pass

    def getDataFromNBODY6(self,folder):
        collnums, namei1, namei2, m1,m2,mnew,time = read_collision_file(folder)
        main_collnums, main_names , main_masses,main_time, rco_name = get_main_collisions(collnums,namei1,namei2,m1,m2,mnew,time)
        collision_track = get_collision_data(folder,main_collnums,main_names,main_masses,main_time)
        return rco_name,collision_track

    def getLastOutFile(self,folder):
        fls = fnmatch.filter(os.listdir(folder),'out*.sph')
        fls.sort()
        return fls[-1]

    def getBoundMass(self,filename):
        Pfilename = filename.ljust(32)
        mass1 = c_double(0)
        mass2 = c_double(0)
        mass3 = c_double(0)
        self.toolsLibrary.getmergerproduct_(Pfilename,byref(mass1),byref(mass2),byref(mass3))


        # mass1 > mass2 > mass3

        if(abs(mass1.value) <= 1.e-3):
            mass1 = None
        else:
            mass1 = mass1.value
        if(abs(mass2.value) <= 1.e-3):
            mass2 = None
        else:
            mass2 = mass2.value
        if(abs(mass3.value) <= 1.e-3):
            mass3 = None
        else:
            mass3 = mass3.value

        return mass1,mass2,mass3


    def __getGamForStar(self):
        if (self.starmass >= 1.0):
            mu_m = 0.617534  # meam molecular weight in proton masses
            m = self.starmass
            m_c = 2.28      # obtained from equation 2.3(b) "Collisions and Close Encounters Between"
                            # "Massive Main Sequence Stars" by Lai Et al. 1993

            alpha_gam = mu_m*mu_m*m/m_c     # eqn 2.3(a) from Lai Et al. 1993

            # eqn 2.4 from Lai Et al. 1993 follows

            c = 7.89 / alpha_gam

            beta_gam_sq = (-(c**2) + (c**4 + 4*c**2) ** 0.5) /2
            beta_gam = beta_gam_sq ** 0.5

            gam_dyn = (32 - 24*beta_gam - 3*(beta_gam**2)) / (24 - 21*beta_gam)

            self.gam = gam_dyn
        else:
            self.gam = (5.0)/(3.0)

    def __setAndBcast(self,comm):

        if (not os.path.isdir(self.dirname)):
            os.system("mkdir " + self.dirname)


        if(self.simulationtype == 'dbl'):
            inputs, input_structtype, inputs_double_s1, inputs_double_s2, input_double_structtype = self.__setValuesStarsmasher()
        else:
            inputs, input_structtype = self.__setValuesStarsmasher()

        Pstartfile1 = self.startfile1.ljust(255)
        Pstartfile2 = self.startfile2.ljust(255)

        Peosfile = self.eosfile.ljust(255)
        Popacityfile = self.opacityfile.ljust(255)

        Pprofilefile = self.profilefile.ljust(255)

        Psimulationtype = self.simulationtype.ljust(3)
        Pdirname = self.dirname.ljust(32)

        comm.Bcast([inputs,input_structtype],root=MPI.ROOT)
        #comm.Send([inputs,input_structtype],dest=0,tag=1)

        comm.Bcast([Pstartfile1,MPI.CHAR],root=MPI.ROOT)
        comm.Bcast([Pstartfile2,MPI.CHAR],root=MPI.ROOT)
        comm.Bcast([Peosfile,MPI.CHAR],root=MPI.ROOT)
        comm.Bcast([Popacityfile,MPI.CHAR],root=MPI.ROOT)
        comm.Bcast([Pprofilefile,MPI.CHAR],root=MPI.ROOT)
        comm.Bcast([Psimulationtype,MPI.CHAR],root=MPI.ROOT)
        comm.Bcast([Pdirname,MPI.CHAR],root=MPI.ROOT)

        if (self.simulationtype == 'dbl'):
            comm.Bcast([inputs_double_s1,input_double_structtype],root=MPI.ROOT)
            comm.Bcast([inputs_double_s2,input_double_structtype],root=MPI.ROOT)


    def run(self):
        self.__errorHandler()

        worker = os.path.abspath("src/worker_code")
        comm = MPI.COMM_SELF.Spawn(worker,args=[],maxprocs=self.num_of_workers)
        self.__setAndBcast(comm)
        comm.barrier()
        MPI.Finalize()

