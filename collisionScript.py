from PyStarsmasher import Starsmasher
import os
import sys

def relax_star(code, name,mass, cnum=None):
    code.nrelax = 1

    code.treloff = 100 * (code.runit ** 3/code.munit) ** 0.5
    code.trelax = (code.runit **3 / code.munit) ** 0.5

    #if (mass >= 20.0):
    #    code.gam = 1.5
    #else:
    #    code.gam = 5.0/3.0

    code.tf = 300
    rad = mass ** 0.8

    #if (mass > 50):
    #    code.n = 7000
    #else:
    #    code.n = 1000

    code.n = int((14000.0/32.0)*rad)

    #code.n = 7000
    if(name != -1):
        code.dirname = 'relax_'+str(name)
    else:
        code.dirname = 'relax_coll_'+str(cnum)
    code.simulationtype = '1es'
    code.starradius = mass ** 0.8
    if (mass < 0.6):
        code.npoly = 1.5
    else:
        code.npoly = 3
    code.starmass = mass
    #code.setParams()
    code.run(num_of_wokers=16)

def collide_stars(cnum,code,file1,file2,s1,s2):
    code.tf = 600

    #code.treloff = 100 * (code.runit ** 3/code.munit) ** 0.5
    #code.trelax = (code.runit **3 / code.munit) ** 0.5


    code.nrelax = 0
    code.simulationtype = 'dbl'

    code.dirname = 'coll_'+str(cnum)

    code.startfile1 = file1
    code.startfile2 = file2

    code.star1.x = s1[2][0]
    code.star1.y = s1[2][1]
    code.star1.z = s1[2][2]

    code.star1.vx = s1[3][0]
    code.star1.vy = s1[3][1]
    code.star1.vz = s1[3][2]

    code.star1.m = s1[4]

    code.star2.x = s2[2][0]
    code.star2.y = s2[2][1]
    code.star2.z = s2[2][2]

    code.star2.vx = s2[3][0]
    code.star2.vy = s2[3][1]
    code.star2.vz = s2[3][2]

    code.star2.m = s2[4]


    # approximate how long the code needs to run for
    # a collision
    rel_pos2 = (code.star2.x - code.star1.x) ** 2 +(code.star2.y - code.star1.y) ** 2 +(code.star2.z - code.star1.z) ** 2

    rel_pos = (rel_pos2 ** 0.5) * 695510.0

    rel_vel2 = (code.star2.vx - code.star1.vx) ** 2 +(code.star2.vy - code.star1.vy) ** 2 +(code.star2.vz - code.star1.vz) ** 2

    rel_vel = rel_vel2 ** 0.5

    tunit = 1592.48283336 #in seconds

    n_tunits = int((rel_pos/rel_vel)/tunit)

    code.tf = 4000

    #code.setParams()
    code.run(num_of_workers=16)


def run():
    folder = '/data/mukherjeed2/nbody6_starsmasher_test/N_8000_S_0.8_chain'
    #file_out = open("output.txt",'a')
    code = Starsmasher()

    code.dtout = 1
    code.n = 7000
    code.nrelax = 1
    code.ngravprocs = 4
    code.computeexclusivemode = 0
    code.ppn = 24
    code.npoly = 3

    rco_name , collision_track = code.getDataFromNBODY6(folder)

    coll_count = 0
    rco_has_been_relaxed = False

    init_num = 0
    restartFileExists = os.path.isfile('restart.starsmasher')

    attemptRestart = False


    #if(restartFileExists):
    #    frestart = open('restart.starsmasher','r')
    #    cs = frestart.readlines()
    #    lastColl = int(cs[-1])
    #    init_num = lastColl*2
    #    coll_count = lastColl
    #    frestart.close()
        #rco_has_been_relaxed = True
    #    attemptRestart = True

    #print(attemptRestart)
    #sys.exit()

    for i in xrange(init_num,len(collision_track), 2):
        frestart = open('restart.starsmasher','w')
        frestart.write(str(coll_count)+'\n')
        frestart.close()



        s1 = collision_track[i]
        s2 = collision_track[i+1]

        # Tracks which star is new and has to be relaxed
        ns = None
        merger = None
        rco = None
        fs = None

        # First time. Both stars have to be relaxed before collision
        if (not rco_has_been_relaxed):
            #file_out.write("relaxing first set of stars ...\n")
            #file_out.write("relax: "+str(s1[1])+" "+str(s2[1]))
            #if (not attemptRestart):
            relax_star(code,s1[1],s1[4])
            relax_star(code,s2[1],s2[4])

            rco_has_been_relaxed = True
            #file_out.write("stars have been relaxed properly ...\n")

            if (s1[1] == rco_name):
                rco = s1
                fs = s2
            else:
                rco = s2
                fs = s1
        if(rco_has_been_relaxed):
            if(s1[1] != rco_name):
                #file_out.write("relaxing new star")
                #file_out.write("name = " + str(s1[1]))
                #if(not attemptRestart):
                relax_star(code,s1[1],s1[4])
                ns = s1
                merger = s2
            else:
                #file_out.write("relaxing new star")
                #file_out.write("name = " + str(s1[1]))
                #if(not attemptRestart):
                relax_star(code,s2[1],s2[4])
                ns = s2
                merger = s1

        # First collision is a little different
        # It relaxes both stars and then collides them
        if(coll_count == 0):
            #file_out.write("colliding first set of stars ...\n")
            #file_out.write("coll_count = "+str(coll_count)+'\n')

            file1 = 'relax_'+str(rco[1])+'/out0300.sph'
            file2 = 'relax_'+str(fs[1])+'/out0300.sph'

            #file_out.write("stars are = "+file1 + " "+ file2+'\n')

            collide_stars(coll_count, code, file1,file2, rco, fs)

        elif(coll_count > 0 and (ns[5]-collision_track[i-1][5]) <= 1.0e-14):
            #file_out.write("colliding stars ... Time difference less than 1 time unit\n")
            #file_out.write("coll_coount = "+str(coll_count)+'\n')

            out_folder = 'coll_'+str(coll_count-1)
            fil = code.getLastOutFile(out_folder)
            file1 = out_folder+'/'+fil
            file2 = 'relax_'+str(ns[1])+'/out0300.sph'

            #file_out.write("stars are = "+file1 + " "+ file2+'\n')

            collide_stars(coll_count, code, file1, file2, merger, ns)
        elif(coll_count > 0 and (ns[5]-collision_track[i-1][5]) > 1.0e-14):
            #create new star from the bound mass remaining
            #file_out.write("colliding stars ... Time difference more than 1 time unit\n")
            #file_out.write("coll_coount = "+str(coll_count)+'\n')

            out_folder = 'coll_'+str(coll_count-1)
            fil = code.getLastOutFile(out_folder)
            fmerger = out_folder+'/'+fil

            boundMass1, boundMass2, boundMass3 = code.getBoundMass(fmerger)
            relax_star(code,-1, boundMass1, coll_count)

            file1 = 'relax_coll_'+str(coll_count)+'/out0300.sph'
            file2 = 'relax_'+str(ns[1])+'/out0300.sph'
            #file_out.write("stars are = "+file1 + " "+ file2+'\n')

            collide_stars(coll_count,code,file1,file2, merger, ns)

        coll_count+=1
    #file_out.close()

if __name__ == '__main__':
    run()

