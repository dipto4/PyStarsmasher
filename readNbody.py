from __future__ import print_function
import os
import numpy as np

def check_distance_between_stars(star_data,flag):
    star1 = star_data[0]
    star2 = star_data[1]

    pos1 = star1[2]
    pos2 = star2[2]

    sep2 = (pos1[0]-pos2[0])**2 + (pos1[1]-pos2[1])**2 + (pos1[2]-pos2[2])**2
    sep = sep2 ** 0.5
    #sep = (sep * 3.086e13) / 695510.0
    m1 = star1[4]
    m2 = star2[4]

    r1r2 = (m1 ** 0.8) + (m2 ** 0.8)
    #print("sep = ",sep)
    #print("r1r2 = ",r1r2)

    if (flag!= 1):
        if(sep <= 10*r1r2 and sep >= 2*r1r2):
            return True
        else:
            return False
    else:
        if(sep <= 2*r1r2):
            return True
        else:
            return False



def read_collision_file(folder):
    # these get returned
    collnums = []
    namei1 = []
    namei2 = []
    m1 = []
    m2 = []
    mnew = []
    time = []

    fname = '{}/coldat.69'.format(folder)
    print(fname)
    fexist = os.path.isfile(fname)

    if(not fexist):
        print("coldat.69 does not exist!")
        return

    f = open(fname)
    linecount = 0
    for line in f:
        if(linecount <= 1):
            linecount+=1
            continue

        # read in data
        # COLLNUM, TIME, NAME(I1), NAME(I2), M(I1)[M*], M(I2)[M*], X1(I1) X2(I1) X3(I1), X1(I2) X2(I2) X3(I2), V1(I1) V2(I1) V3(I1), V1(I2) V2(I2) V3(I2)
        # RS(I1)[R*]  RS(I2)[R*]

        params = line.split()
        #print(params)
        collnum = int(params[0])
        collnums.append(collnum)
        time.append(float(params[1]))
        i1 = int(params[2])
        i2 = int(params[3])
        namei1.append(i1)
        namei2.append(i2)

        mcoll = float(params[4]) + float(params[5])
        mnew.append(mcoll)

        m1.append(float(params[4]))
        m2.append(float(params[5]))


        linecount+=1

    if(linecount <= 1):
        print("corrupted file!")
        return

    f.close()
    return collnums, namei1, namei2, m1,m2, mnew, time


def get_main_collisions(collnums,namei1,namei2,m1,m2,mnew,time):
    # these get returned
    main_collnums = []
    main_names = []
    main_masses = []
    main_time = []

    rco_index = mnew.index(max(mnew))

    temp_mass_i1 = m1[rco_index]
    temp_mass_i2 = m2[rco_index]

    rco_name = None

    if(temp_mass_i1 > temp_mass_i2):
        rco_name = namei1[rco_index]
    else:
        rco_name = namei2[rco_index]

    np_namei1 = np.asarray(namei1)
    np_namei2 = np.asarray(namei2)

    indices_as_i1 = np.where(rco_name == np_namei1)[0]
    indices_as_i2 = np.where(rco_name == np_namei2)[0]

    for i in indices_as_i1:
        c = collnums[i]
        t = time[i]
        main_collnums.append(c)
        main_time.append(t)

    for i in indices_as_i2:
        c = collnums[i]
        t = time[i]
        main_collnums.append(c)
        main_time.append(t)

    main_collnums.sort()
    main_time.sort()

    for c in main_collnums:
        i = collnums.index(c)
        main_names.append((namei1[i],namei2[i]))
        if(rco_name == namei1[i]):
            if(m2[i] > m1[i]):
                main_masses.append((m2[i],m1[i]))
            else:
                main_masses.append((m1[i],m2[i]))
        else:
            if(m1[i] > m2[i]):
                main_masses.append((m2[i],m1[i]))
            else:
                main_masses.append((m1[i],m2[i]))

    return main_collnums, main_names, main_masses, main_time, rco_name

def get_collision_data(folder,main_collnums,main_names,main_masses, main_time):
    #get appropriate collision_track.469_<collnum>_number file
    track_collision = []

    count = 0
    for c in main_collnums:
        collision_file = "collision_track.469_{}_".format(str(c).zfill(2))
        collision_file = folder + "/" + collision_file
        t = main_time[count]
        for i in xrange(5,6):

            collision_filei = collision_file+str(i).zfill(1)
            flag = 0
            f = open(collision_filei,'r')
            pos_temp=[]
            nameCount = 0
            for line in f:
                params = line.split()
                name = int(params[0])
                x = float(params[1]) * 3.086e13 / 695510.0
                y = float(params[2]) * 3.086e13 / 695510.0
                z = float(params[3]) * 3.086e13 / 695510.0
                vx = float(params[4])
                vy = float(params[5])
                vz = float(params[6])
                m = float(params[7])

                # if there is a chain present
                if(name == 0):
                    m1 = main_masses[count][0]
                    m2 = main_masses[count][1]
                    if(abs(m-m1) <= 1.e-6):
                        name = main_names[count][0]
                        pos_temp.append((c,name,[x,y,z],[vx,vy,vz],m,t))
                        flag =2
                    elif(abs(m-m2) <= 1.e-6):
                        name = main_names[count][1]
                        pos_temp.append((c,name,[x,y,z],[vx,vy,vz],m,t))
                        flag = 3
                    else:
                        #print("here lies the problem")
                        pos_temp = []
                        flag = 4

                    break

                if(name == main_names[count][0]):
                    m = main_masses[count][0]
                    pos_temp.append((c,name,[x,y,z],[vx,vy,vz],m,t))
                    nameCount+=1
                elif (name == main_names[count][1]):
                    m = main_masses[count][1]
                    pos_temp.append((c,name,[x,y,z],[vx,vy,vz],m,t))
                    nameCount+=1


                if (nameCount == 2):
                    #print("do i get here?")
                    flag = 1
                    break
            #print(pos_temp)
            if(flag == 1):
                #if(check_distance_between_stars(pos_temp,flag) == True):
                if(True):
                    #print(pos_temp)
                    track_collision.append(pos_temp[0])
                    track_collision.append(pos_temp[1])
                    #print("found!",c,i)
                    break

            # flag = 2, flag = 3if particle 0 corresponds to rco
            if (flag == 2 or flag == 3 or flag == 4):
                file_temp = open(collision_filei,'r')
                lines_temp = file_temp.readlines()

                #file_size = len(lines_temp)

                part_lines = lines_temp[-3:]

                for p in part_lines:
                    params = p.split()
                    name = int(params[0])
                    x = float(params[1]) * 3.086e13 / 695510.0
                    y = float(params[2]) * 3.086e13 / 695510.0
                    z = float(params[3]) * 3.086e13 / 695510.0
                    vx = float(params[4])
                    vy = float(params[5])
                    vz = float(params[6])

                    if (flag == 2):
                        if(name == main_names[count][1]):
                            m = main_masses[count][1]
                            pos_temp.append((c,name,[x,y,z],[vx,vy,vz],m,t))
                    elif (flag == 3):
                        if(name == main_names[count][0]):
                            m = main_masses[count][0]
                            pos_temp.append((c,name,[x,y,z],[vx,vy,vz],m,t))
                    elif (flag == 4):
                        #print("problem zone 1")
                        #print(name)
                        #print(main_names[count][1])
                        if (name == main_names[count][0]):
                            #print("problem zone 2")
                            m = main_masses[count][0]
                            pos_temp.append((c,name,[x,y,z],[vx,vy,vz],m,t))
                            #print(pos_temp)
                        if (name == main_names[count][1]):
                            #print("problem zone 3")
                            m = main_masses[count][1]
                            pos_temp.append((c,name,[x,y,z],[vx,vy,vz],m,t))
                    else:
                        print("Something has seriously gone wrong!")

                file_temp.close()

                #if(check_distance_between_stars(pos_temp,flag) == True):
                if(True):
                    #print("chain")
                    #print(pos_temp)
                    #print(i)
                    track_collision.append(pos_temp[0])
                    track_collision.append(pos_temp[1])
                    break

            # flag = 4, particle 0 does not correspond to rco
        count+=1
    return track_collision





