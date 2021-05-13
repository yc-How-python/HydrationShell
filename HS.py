#import
# =============================================================================
# =============================================================================
from math import acos,fabs,radians,trunc,degrees,cos
import os
import sys
import numpy as np
import ctypes



import multiprocessing as mp
import time


def distat(atomi,atomj):
    disx=fabs(atomi.x-atomj.x)
    disy=fabs(atomi.y-atomj.y)
    disz=fabs(atomi.z-atomj.z)
    dis=pow(disx**2+disy**2+disz**2,1/2)
    return dis


def pbc(x,y,z):


    (x1,x2)=x.split(maxsplit=2)

    (y1,y2)=y.split(maxsplit=2)

    (z1,z2)=z.split(maxsplit=2)
    #read boundary
    LX=float(x2)-float(x1)
    LY=float(y2)-float(y1)
    LZ=(float(z2)-float(z1))
    return (LX, LY, LZ, x2, x1, y2, y1, z2, z1)





def pro_check(matrix,box):#[LX, LY, LZ, x2, x1, y2, y1, z2, z1]
    [LX, LY, LZ, x2, x1, y2, y1, z2, z1] = box
    listcopy = []
    dO_O = 3.5
    for atomsarray in matrix:
        dx = fabs(atomsarray[0][1] - float(x2))
        dy = fabs(atomsarray[0][2] - float(y1))
        dz = fabs(atomsarray[0][3] - float(z1))

        if (dx <= dO_O or dy <= dO_O or dz <= dO_O):
            if dx <= dO_O:
                movex = -LX

                k = [atomsarray[0][0], atomsarray[0][1] + movex, atomsarray[0][2], atomsarray[0][3]]
                kh1 = [atomsarray[1][0], atomsarray[1][1] + movex, atomsarray[1][2], atomsarray[1][3]]
                kh2 = [atomsarray[2][0], atomsarray[2][1] + movex, atomsarray[2][2], atomsarray[2][3]]

                listcopy.append([k, kh1, kh2])
            if dy <= dO_O:
                movey = LY
                k = [atomsarray[0][0], atomsarray[0][1], atomsarray[0][2] + movey, atomsarray[0][3]]

                kh1 = [atomsarray[1][0], atomsarray[1][1], atomsarray[1][2] + movey, atomsarray[1][3]]

                kh2 = [atomsarray[2][0], atomsarray[2][1], atomsarray[2][2] + movey, atomsarray[2][3]]

                listcopy.append([k, kh1, kh2])
            if dz <= dO_O:
                movez = LZ
                k = [atomsarray[0][0], atomsarray[0][1], atomsarray[0][2], atomsarray[0][3] + movez]

                kh1 = [atomsarray[1][0], atomsarray[1][1], atomsarray[1][2], atomsarray[1][3] + movez]

                kh2 = [atomsarray[2][0], atomsarray[2][1], atomsarray[2][2], atomsarray[2][3] + movez]

                listcopy.append([k, kh1, kh2])
            if dx <= dO_O and dy <= dO_O:
                movex = -LX
                movey = LY
                movez = 0
                k = [atomsarray[0][0], atomsarray[0][1] + movex, atomsarray[0][2] + movey, atomsarray[0][3] + movez]

                kh1 = [atomsarray[1][0], atomsarray[1][1] + movex, atomsarray[1][2] + movey, atomsarray[1][3] + movez]

                kh2 = [atomsarray[2][0], atomsarray[2][1] + movex, atomsarray[2][2] + movey, atomsarray[2][3] + movez]

                listcopy.append([k, kh1, kh2])
            if dx <= dO_O and dz <= dO_O:
                movex = -LX
                movez = LZ
                movey = 0
                k = [atomsarray[0][0], atomsarray[0][1] + movex, atomsarray[0][2] + movey, atomsarray[0][3] + movez]
                kh1 = [atomsarray[1][0], atomsarray[1][1] + movex, atomsarray[1][2] + movey, atomsarray[1][3] + movez]

                kh2 = [atomsarray[2][0], atomsarray[2][1] + movex, atomsarray[2][2] + movey, atomsarray[2][3] + movez]

                listcopy.append([k, kh1, kh2])
            if dy <= dO_O and dz <= dO_O:
                movex = 0
                movey = LY
                movez = LZ
                k = [atomsarray[0][0], atomsarray[0][1] + movex, atomsarray[0][2] + movey, atomsarray[0][3] + movez]

                kh1 = [atomsarray[1][0], atomsarray[1][1] + movex, atomsarray[1][2] + movey, atomsarray[1][3] + movez]

                kh2 = [atomsarray[2][0], atomsarray[2][1] + movex, atomsarray[2][2] + movey, atomsarray[2][3] + movez]

                listcopy.append([k, kh1, kh2])
            if dx <= dO_O and dy <= dO_O and dz <= dO_O:
                movex = -LX
                movey = LY
                movez = LZ
                k = [atomsarray[0][0], atomsarray[0][1] + movex, atomsarray[0][2] + movey, atomsarray[0][3] + movez]

                kh1 = [atomsarray[1][0], atomsarray[1][1] + movex, atomsarray[1][2] + movey, atomsarray[1][3] + movez]

                kh2 = [atomsarray[2][0], atomsarray[2][1] + movex, atomsarray[2][2] + movey, atomsarray[2][3] + movez]

                listcopy.append([k, kh1, kh2])


        else:
            pass

    return np.array(listcopy)
def pro_checkC(matrix,box):#[LX, LY, LZ, x2, x1, y2, y1, z2, z1]
    [LX, LY, LZ, x2, x1, y2, y1, z2, z1] = box
    listcopy = []
    dO_O = 3.5
    for atomsarray in matrix:
        dx = fabs(atomsarray[1] - float(x2))
        dy = fabs(atomsarray[2] - float(y1))
        dz = fabs(atomsarray[3] - float(z1))

        if (dx <= dO_O or dy <= dO_O or dz <= dO_O):
            if dx <= dO_O:
                movex = -LX

                k = [atomsarray[0], atomsarray[1] + movex, atomsarray[2], atomsarray[3]]

                listcopy.append(k)
            if dy <= dO_O:
                movey = LY
                k = [atomsarray[0], atomsarray[1], atomsarray[2] + movey, atomsarray[3]]

                listcopy.append(k)
            if dz <= dO_O:
                movez = LZ
                k = [atomsarray[0], atomsarray[1], atomsarray[2], atomsarray[3] + movez]

                listcopy.append(k)
            if dx <= dO_O and dy <= dO_O:
                movex = -LX
                movey = LY
                movez = 0
                k = [atomsarray[0], atomsarray[1] + movex, atomsarray[2] + movey, atomsarray[3] + movez]

                listcopy.append(k)
            if dx <= dO_O and dz <= dO_O:
                movex = -LX
                movez = LZ
                movey = 0
                k = [atomsarray[0], atomsarray[1] + movex, atomsarray[2] + movey, atomsarray[3] + movez]
                listcopy.append(k)
            if dy <= dO_O and dz <= dO_O:
                movex = 0
                movey = LY
                movez = LZ
                k = [atomsarray[0], atomsarray[1] + movex, atomsarray[2] + movey, atomsarray[3] + movez]

                listcopy.append(k)
            if dx <= dO_O and dy <= dO_O and dz <= dO_O:
                movex = -LX
                movey = LY
                movez = LZ
                k = [atomsarray[0], atomsarray[1] + movex, atomsarray[2] + movey, atomsarray[3] + movez]

                listcopy.append(k)


        else:
            pass

    return np.array(listcopy)
def data_extractor(line):
    (aid, atype, ax, ay, az) = line.split(maxsplit=4)
    return [float(aid), float(ax), float(ay), float(az)]
def hydrogen_reset(H1, O, LX, LY, LZ):
    dO_O = 3.5
    if ((O[1] - 3.5 <= H1[1] <= O[1] + 3.5) and (O[2] - 3.5 <= H1[2] <= O[2] + 3.5) and (
            O[3] - 3.5 <= H1[3] <= O[3] + 3.5)) == False:
        # print('found situation1')#将H移动至O附近
        distH1 = [H1[1] - O[1], H1[2] - O[2], H1[3] - O[3]]
        if fabs(distH1[0]) > dO_O:
            if distH1[0] > 0:
                movex = -LX
            else:
                movex = LX
        else:
            movex = 0
        if fabs(distH1[1]) > dO_O:
            if distH1[1] > 0:
                movey = -LY
            else:
                movey = LY
        else:
            movey = 0
        if fabs(distH1[2]) > dO_O:
            if distH1[2] > 0:
                movez = -LZ
            else:
                movez = LZ
        else:
            movez = 0
        H1K = [H1[0], H1[1] + movex, H1[2] + movey, H1[3] + movez]
        return H1K
    else:
        return H1
def Loading_atom(target,datafile):
    listatom = []
    listC=[]
    while 1:

        line = datafile.readline()
        if line == 'ITEM: NUMBER OF ATOMS\n':  # 确认原子数
            Noatom_str = datafile.readline()
            Noatom = int(Noatom_str)
            # print('Number of atoms :  ', Noatom)
        elif line == 'ITEM: TIMESTEP\n':  # timestep
            title = datafile.readline()
            # print('#' * 5 ** 2, '\n')
            # print('timestep:  ', title)
        elif line == 'ITEM: BOX BOUNDS pp pp pp\n':
            (LX, LY, LZ, x2, x1, y2, y1, z2, z1) = pbc(datafile.readline(), datafile.readline(), datafile.readline())
        elif line == 'ITEM: ATOMS id type xu yu zu\n':  # atomdata
            roundnum = 1
            while roundnum <= Noatom:

                atomstr = datafile.readline()

                if atomstr == '':
                    break
                else:
                    # print(line)
                    (aid, atype, ax, ay, az) = atomstr.split(maxsplit=4)
                if atype == target :
                    O = [float(aid), float(ax), float(ay), float(az)]
                    # print('发现'+ax+' '+ay+' '+az+' '+atype+' \n')#check
                    H1 = data_extractor(datafile.readline())
                    H2 = data_extractor(datafile.readline())
                    H1 = hydrogen_reset(H1, O, LX, LY, LZ)
                    H2 = hydrogen_reset(H2, O, LX, LY, LZ)
                    listatom.append([O, H1, H2])
                    roundnum = roundnum + 3
                    continue
                elif atype == '4':
                    listC.append([float(aid), float(ax), float(ay), float(az)])
                    roundnum += 1
                    continue
                else:
                    roundnum += 1
                    continue
            # print('共', len(listatom), '个原子')
            matrix = np.array(listatom)
            matrixC=np.array(listC)
            return [matrix, title, [LX, LY, LZ, x2, x1, y2, y1, z2, z1], matrixC]

            break  # while1
#torsion angle calculator



def PHI(wateri, waterj,OH):  # 原子i与其邻num之间二面角，参数为，i的第ih个oh和num的第neih个oh

    vij=OH
    costheta=np.dot( vij, [0,0,1])/((vij[0] ** 2 + vij[1] ** 2 + vij[2] ** 2) ** (0.5))
    # vz=[0,0,1]
    deg=np.rad2deg(np.arccos(costheta))
    if deg>90:
        return 180-deg
    return deg


def vector(atomi, atomj):
    disx = atomi[1] - atomj[1]
    disy = atomi[2] - atomj[2]
    disz = atomi[3] - atomj[3]
    return (disx, disy, disz)
def anglecheck(wateri,waterj,Onorm):
    #vector
    vO_O = vector(waterj[0], wateri[0])  # j-i
    # print(Onorm)
    v1 = vector(wateri[1], wateri[0])

    v2 = vector(wateri[2], wateri[0])
    v3 = vector(waterj[0], waterj[1])
    v4 = vector(waterj[0], waterj[2])  # ) [H1I-OI,H2I-OI,OJ-H1J,OJ-H2J]
    dotarray=np.dot([v1,v2,v3,v4],vO_O)#[H1I*O_O,H2I*O_O......]
    normarray=np.linalg.norm([v1,v2,v3,v4],axis=1)#沿着y轴进行求模
    O_Onorm=np.linalg.norm(vO_O)#O_O模
    costheta_array=(dotarray/normarray)/O_Onorm #numpy 向量除法得向量
    theta_array=np.arccos(np.round(costheta_array,3))
    thetalist=list(theta_array)
    minvalue=theta_array.min()
    HBvectors=[vector(wateri[1], waterj[0]),vector(wateri[2], waterj[0]),vector(waterj[1], wateri[0]),vector(waterj[2], wateri[0])]
    HB=HBvectors[thetalist.index(minvalue)] #返回氢键所代表的OH

    return [minvalue,HB]



def search_atom(matrix,matrixC,listcopy):
    listtheta=[]
    listalphaLeft=[]
    listalphaRight = []
    dO_O=3.5
    dO_O2=3.5*3.5
    anglecut=np.radians(30)

    for i in range(len(matrixC)): #原子i


        idi = matrixC[i][0]
        Ci=matrixC[i]
        listj = [] #主原子本身的邻

        count=0

        j=np.where(np.square(matrix[:,0,1:] - Ci[1:]).sum(axis=1) < 5.5 **2)[0]

        for order in j:
            OW=matrix[order][0]
            HW1 = matrix[order][1]
            HW2= matrix[order][2]

            vij = np.sum((vector(HW1, OW), vector(HW2, OW)), axis=0)
            H1= vector(HW1, OW)
            H2 = vector(HW2, OW)
            MeO = vector(Ci, OW)
            costheta = np.dot(vij, MeO) / np.linalg.norm( vij)/ (np.linalg.norm(MeO))
            theta=np.rad2deg(round(np.arccos(costheta), 3))
            listtheta.append(theta)
            cosalpha1 = np.dot(H1, MeO) / np.linalg.norm(H1) / (np.linalg.norm(MeO))
            cosalpha2 = np.dot(H2, MeO) / np.linalg.norm(H2) / (np.linalg.norm(MeO))
            alpha1 = np.rad2deg(round(np.arccos(cosalpha1), 3))
            alpha2 = np.rad2deg(round(np.arccos(cosalpha2), 3))
            if OW[-1]>=Ci[-1]:

                listalphaRight.append(alpha1 )
                listalphaRight.append(alpha2)
            else:
                listalphaLeft.append(alpha1 )
                listalphaLeft.append(alpha2)

        if listcopy.shape[0]>1:
            whereCopy = np.where(listcopy[:, 0] == idi)[0]
            for i in whereCopy:

                j = np.where(np.square(matrix[:, 0, 1:] - listcopy[i][1:]).sum(axis=1) < 5.5 ** 2)[0]
                for order in j:
                    OW = matrix[order][0]
                    HW1 = matrix[order][1]
                    HW2 = matrix[order][2]

                    vij = np.sum((vector(HW1, OW), vector(HW2, OW)), axis=0)
                    MeO = vector(listcopy[i], OW)
                    costheta = np.dot(vij, MeO) / np.linalg.norm(vij) / (np.linalg.norm(MeO))
                    theta = np.rad2deg(round(np.arccos(costheta), 3))

                    listtheta.append(theta)
                    cosalpha1 = np.dot(H1, MeO) / np.linalg.norm(H1) / (np.linalg.norm(MeO))
                    cosalpha2 = np.dot(H2, MeO) / np.linalg.norm(H2) / (np.linalg.norm(MeO))
                    alpha1 = np.rad2deg(round(np.arccos(cosalpha1), 3))
                    alpha2 = np.rad2deg(round(np.arccos(cosalpha2), 3))
                    if OW[-1] >= listcopy[i][-1]:

                        listalphaRight.append(alpha1)
                        listalphaRight.append(alpha2)
                    else:
                        listalphaLeft.append(alpha1)
                        listalphaLeft.append(alpha2)
        elif listcopy.shape[0]==1 :
            whereCopy = np.where(listcopy[0] == idi)[0]
            for i in whereCopy:

                j = np.where(np.square(matrix[:, 0, 1:] - listcopy[i][1:]).sum(axis=1) < 5.5 ** 2)[0]
                for order in j:
                    OW = matrix[order][0]
                    HW1 = matrix[order][1]
                    HW2 = matrix[order][2]

                    vij = np.sum((vector(HW1, OW), vector(HW2, OW)), axis=0)
                    MeO = vector(listcopy[i], OW)
                    costheta = np.dot(vij, MeO) / np.linalg.norm(vij) / (np.linalg.norm(MeO))
                    theta = np.rad2deg(round(np.arccos(costheta), 3))

                    listtheta.append(theta)
                    cosalpha1 = np.dot(H1, MeO) / np.linalg.norm(H1) / (np.linalg.norm(MeO))
                    cosalpha2 = np.dot(H2, MeO) / np.linalg.norm(H2) / (np.linalg.norm(MeO))
                    alpha1 = np.rad2deg(round(np.arccos(cosalpha1), 3))
                    alpha2 = np.rad2deg(round(np.arccos(cosalpha2), 3))
                    if OW[-1] >= listcopy[i][-1]:

                        listalphaRight.append(alpha1)
                        listalphaRight.append(alpha2)
                    else:
                        listalphaLeft.append(alpha1)
                        listalphaLeft.append(alpha2)
        else:
            pass


         # costheta = np.dot(vij, MeO) / ((vij[0] ** 2 + vij[1] ** 2 + vij[2] ** 2) ** (0.5)) / (np.linalg.norm(MeO))


        #print('counttol   ',count+mcount)


    print(len(listalphaLeft))
    print(len(listalphaRight))
    allr=[listalphaLeft, listalphaRight]
    return allr
def is_suffix_lammpstrj(suffix: str):
    if suffix == 'lammpstrj':
        return True
    return True




import argparse
import matplotlib.pyplot as plt

def hbond(listdata):
    listl=[]
    listr=[]
    for data in listdata:
        listCcopy = pro_checkC(data[-1], data[2])
        # [matrix, title, [LX, LY, LZ, x2, x1, y2, y1, z2, z1]]
        allr=search_atom(data[0],data[-1], listCcopy)
        listl.extend(allr[0])
        listr.extend(allr[1])
    return [listl,listr]





if __name__ == '__main__':
    fig=plt.figure(figsize=(3,3))
    ax=fig.gca()
    fig2 = plt.figure(figsize=(3, 3))
    ax2 = fig2.gca()
    from matplotlib import cm
    C1 = cm.get_cmap("Greens", 50)
    C2 = cm.get_cmap("Blues", 50)
    C3 = cm.get_cmap("Purples", 50)
    C4 = cm.get_cmap("Oranges", 50)
    cmap = np.vstack((C1(np.linspace(0.6, 0.9, 6)), C2(np.linspace(0.6, 0.9, 5)), C3(np.linspace(0.6, 0.9, 5)),
                      C4(np.linspace(0.6, 0.9, 5))))
    parser = argparse.ArgumentParser(

        description='Notice:\n' + '\n The program and folder that contains trjs to be calculated should to be in ' + 'the same directory.\n ')
    parser.add_argument('-i', type=str, default=r'F:\TrjData\run2', help="The name of folder that contains input trjs")
    parser.add_argument('-o', type=str, default='Default_F4.txt', help="The name of output file")
    parser.add_argument('-c', type=int, default=mp.cpu_count(), help="Number of cores")
    args = parser.parse_args()
    # Constant
    foldername=args.i
    output=args.o
    CORE=args.c
    pool = mp.Pool(2)

    #
    t0=time.time()

    # Batch

    PROJECT_DIR_PATH = os.path.dirname(os.path.abspath(os.path.abspath(__file__)))
    print('当前路径',PROJECT_DIR_PATH)
    DIR_PATH = os.path.join(PROJECT_DIR_PATH, foldername)
    files = os.listdir(DIR_PATH)
    print(files)
    listR=[]
    listN=[]
    listA=[]
    for filename in files:
        foldername = args.i
        output = args.o

        dO_O = 3.5
        target = '1'
        anglecut = np.radians(30)

        PATH = os.path.abspath((os.path.abspath(filename)))  # 改双斜杠绝对路径

        DIR = os.path.dirname(PATH)  # 获得文件夹路径DIR
        name, suffix = os.path.splitext(PATH)  # 获得后缀与前缀
        name = name.split('/')[-1]  # 去除斜杠

        if is_suffix_lammpstrj(suffix):
            filepath=os.path.join(DIR_PATH, filename)
            datafile = open(filepath, 'r')
            ending = datafile.seek(0, 2)
            datafile.seek(0, 0)
            data = []
            STEPLIST = []
            listdeg = []
            listalphaL=[]
            listalphaR=[]
            FINALRESULT = []

            while 1:
                progress = datafile.tell()
                if ending == progress:
                    print('End')
                    break
                else:
                    # pbar.setValue(progress/ending*100)
                    pass
                data.append(Loading_atom(target, datafile))  # atomarray, title, box

            # 启动多个进程
            listdata=data

            multi_process = []
            totalTasks = len(listdata)  # 任务分配
            averTasks = totalTasks // CORE
            for i in range(CORE):
                if i == CORE - 1:
                    remain = totalTasks % CORE + 1
                else:
                    remain = 0
                multi_process.append(
                    pool.apply_async(hbond, (listdata[i * averTasks:(i + 1) * averTasks + remain],)))

            for each_process in multi_process:
                if each_process.get() == []:
                    r = []



                else:
                    r = each_process.get()

                listalphaL.extend(r[0])
                listalphaR.extend(r[1])

            lb=float(os.path.splitext(filename)[0].replace('_', '.'))
            print('Initiated F4_proc. ')
            listN.append(lb)


            listR.append([listalphaL,listalphaR])

            print('Finished F4_proc. ')

            # outdata.write(str(res)+' \n')
    t1 = time.time()
    Name=np.array(listN)
    Results=np.array(listR)
    sorder=Name.argsort()
    Results = np.array(listR)[sorder]

    Name = Name[sorder]







    from matplotlib.pyplot import MultipleLocator
    def drawsub(xlist, ylist, namelist,bins):
        b = 0
        fig, ax = plt.subplots(5, 1, facecolor='w', figsize=(1, 2), dpi=600, sharex=True, )
        fig2, ax2 = plt.subplots(4, 1, facecolor='w', figsize=(1, 2), dpi=600, sharex=True)
        plt.rc('font', family='Times New Roman', )
        C1 = cm.get_cmap("Greens", 50)
        C2 = cm.get_cmap("Blues", 50)
        C3 = cm.get_cmap("Purples", 50)
        C4 = cm.get_cmap("Oranges", 50)
        markerlist = ['*', '^', '>', 'v', '<', '^', '>', 'v', '<', '^', '>', 'v', '<', '^', '>', 'v', '<', '^', '>',
                      'v',
                      '<']
        cmap = np.vstack((C1(np.linspace(0.6, 0.9, 6)), C2(np.linspace(0.6, 0.9, 5)), C3(np.linspace(0.6, 0.9, 5)),
                          C4(np.linspace(0.6, 0.9, 5))))
        x_major_locator = MultipleLocator(30)

        # ax.xaxis.set_major_locator(x_major_locator)
        # ax2.xaxis.set_major_locator(x_major_locator)

        for order in list(range(len(ylist))):
            if order == 0:
                # s = fig2.add_subplot(510 + order+1,sharex=True,)
                d = np.histogram(ylist[order][0], bins=bins, density=True)
                ax2[order].plot(d[1][:-1], d[0] * 100, c='black', lw=0.7, label='E = ' + str(namelist[order]) + 'V/nm')
                ax2[order].tick_params(labelsize=6, direction='in', width=0.5)
                ax2[order].xaxis.set_major_locator(x_major_locator)
                ax2[order].set_ylim([0, 2])
                d = np.histogram(ylist[order][1], bins=bins, density=True)
                ax2[order].plot(d[1][:-1], d[0] * 100, c='r',lw=0.7, label='E = ' + str(namelist[order]) + 'V/nm')
                ax2[order].tick_params(labelsize=6, direction='in', width=0.5)
                ax2[order].xaxis.set_major_locator(x_major_locator)
                ax2[order].set_ylim([0, 2])
            elif order <= 8:
                d = np.histogram(ylist[order][0], bins=bins, density=True)
           # 计算均值

                num_bins = bins  # 直方图柱子的数量
                # n, bins, patches = ax2.hist(ylist[order], num_bins, density=1, alpha=0.75)
                # # 直方图函数，x为x轴的值，normed=1表示为概率密度，即和为一，绿色方块，色深参数0.5.返回n个概率，直方块左边线的x值，及各个方块对象
                # y = norm.pdf(bins, mu, sigma)
                # print('%s ns: mu = %d sigma = %d'% (str(namelist[order]),mu, sigma))
                # radlimit = np.rad2deg(np.arccos(
                #     np.cos(0) - (3 / 2 * 1.38 * (10 ** (-23)) * 290 / (
                #                 order / 10 * 2.426 * 3.335 * (10 ** (-30)) * (10 ** 9)))))

                if order < 4:
                    # s = fig2.add_subplot(510 + order+1,sharex=True)
                    ax2[order].plot(d[1][:-1], d[0] * 100, lw=0.7, label='E = ' + str(namelist[order - 1]) + 'V/nm',
                                    c=cmap[order],
                                    )
                    ax2[order].tick_params(labelsize=6, direction='in', width=0.5)
                    ax2[order].xaxis.set_major_locator(x_major_locator)
                    ax2[order].set_ylim([0, 2])
                    d = np.histogram(ylist[order][1], bins=bins, density=True)
                    ax2[order].plot(d[1][:-1], d[0] * 100, lw=0.7, label='E = ' + str(namelist[order - 1]) + 'V/nm',
                                    c='r'
                                    )
                    ax2[order].tick_params(labelsize=6, direction='in', width=0.5)
                    ax2[order].xaxis.set_major_locator(x_major_locator)
                    ax2[order].set_ylim([0, 2])
                else:
                    # s = fig.add_subplot(510 + order-3,sharex=True)
                    # ax2.plot(range(len(y)),y*100,label=str(namelist[order]) + 'V/nm')

                    ax[order - 4].plot(d[1][:-1], d[0] * 100, lw=0.7, label='E = ' + str(namelist[order - 1]) + 'V/nm',
                                      c=cmap[order])
                    ax[order - 4].tick_params(labelsize=6, direction='in', width=0.5)
                    ax[order - 4].xaxis.set_major_locator(x_major_locator)
                    ax[order - 4].set_ylim([0, 2])
                    d = np.histogram(ylist[order][1], bins=bins, density=True)
                    ax[order - 4].plot(d[1][:-1], d[0] * 100, lw=0.7, label='E = ' + str(namelist[order - 1]) + 'V/nm',
                                       c='r')
                    ax[order - 4].tick_params(labelsize=6, direction='in', width=0.5)
                    ax[order - 4].xaxis.set_major_locator(x_major_locator)
                    ax[order - 4].set_ylim([0, 2])
                    # ax[order - 4].annotate(round(radlimit, 1), xy=(radlimit, 0.5), xytext=(radlimit, 0.5),
                    #                        arrowprops=dict(facecolor='red', shrink=1, width=0.1), fontsize=4)
        # ax.legend(fontsize=6, )
        # ax2.legend(fontsize=6,)
        # ax.tight_layout(pad=0)
        # ax2.tight_layout(pad=0)

        fig.tight_layout(pad=0)
        fig2.tight_layout(pad=0)


        fig.tight_layout(pad=0)
        fig2.tight_layout(pad=0)
        fig.subplots_adjust(wspace=0, hspace=0)
        fig2.subplots_adjust(wspace=0, hspace=0)
        fig.savefig('alpha1.tiff',dpi=600)
        fig2.savefig('alpha2.tiff', dpi=600)
    drawsub([], Results, Name, 181)


    def drawsum(xlist, ylist, namelist, bins):
        b = 0.75
        fig = plt.figure(figsize=(2, 2), dpi=300)
        ax = fig.gca()
        fig2 = plt.figure(figsize=(2, 2), dpi=300)
        ax2 = fig2.gca()
        plt.rc('font', family='Times New Roman', )

        C1 = cm.get_cmap("Greens", 50)
        C2 = cm.get_cmap("Blues", 50)
        C3 = cm.get_cmap("Purples", 50)
        C4 = cm.get_cmap("Oranges", 50)
        markerlist = ['*', '^', '>', 'v', '<', '^', '>', 'v', '<', '^', '>', 'v', '<', '^', '>', 'v', '<', '^', '>',
                      'v',
                      '<']
        cmap = np.vstack((C1(np.linspace(0.6, 0.9, 6)), C2(np.linspace(0.6, 0.9, 5)), C3(np.linspace(0.6, 0.9, 5)),
                          C4(np.linspace(0.6, 0.9, 5))))

        plt.xlabel('Time(ns)', fontsize=12)
        plt.ylabel('Angle', fontsize=12)

        ax.tick_params(labelsize=8, direction='in', width=0.5)
        ax.spines['bottom'].set_linewidth(b)
        ax.spines['left'].set_linewidth(b)
        ax.spines['top'].set_linewidth(b)
        ax.spines['right'].set_linewidth(b)
        plt.gca().yaxis.set_ticks_position('left')

        for order in list(range(len(ylist))):
            if order > 10:
                continue
            if order == 0:

                d = np.histogram(np.hstack([ylist[order][0], ylist[order][1]]), bins=bins, density=True)
                ax.plot(d[1][:-1], d[0] * 100, c='black', lw=0.5, label='E = ' + str(namelist[order]) + 'V/nm')
                ax.tick_params(labelsize=6, direction='in', width=0.5)

            else:
                # s = fig.add_subplot(510 + order-3,sharex=True)
                # ax2.plot(range(len(y)),y*100,label=str(namelist[order]) + 'V/nm')
                d = np.histogram(np.hstack([ylist[order][0], ylist[order][1]]), bins=bins, density=True)
                ax.plot(d[1][:-1], d[0] * 100, lw=0.5, label='E = ' + str(namelist[order - 1]) + 'V/nm',
                        c=cmap[order], markersize=4, marker=markerlist[order], markevery=90)
        d = np.histogram(np.hstack([ylist[0][0], ylist[0][1]]), bins=bins, density=True)
        ax.plot(d[1][:-1], d[0] * 100, c='black', lw=0.5, label='E = ' + str(namelist[0]) + 'V/nm')
        ax.tick_params(labelsize=6, direction='in', width=0.5)
        # ax2.legend(fontsize=6,)
        # ax.tight_layout(pad=0)
        # ax2.tight_layout(pad=0)

        fig.tight_layout(pad=0)
    drawsum([], Results, Name, 90)
