__author__ = 'wilsonp'
import VGenesSeq
from operator import itemgetter
from math import ceil
import os
def CloneCaller(DataList, Duplicates, signal, poolID):


    # 'SeqName', 'VLocus', 'JLocus', 'CDR3Length', 'CDR3DNA', 'Mutations', 'Vbeg', 'Vend', 'Sequence', 'ID', 'GVend', 'GJbeg', 'GD1beg', 'GD1end'
    # 0             1        2            3           4          5           6       7         8        9       10       11       12        13

    DataList.sort(key=itemgetter(1,2,3)) #sort data list by Vlocus, Jlocus, and CDR3 length
    i = 0
    r = 0
    StartNew = False
    NewList = []
    PrePools = []
    NewPools = []
    for record in DataList:
        # Sequence = record[1]
        CDR3 = record[4]
        # CDR1from = int(record[7])
        # if len(CDR3) > 0:
        if r+1 < len(DataList):
            if StartNew == False:
                NewList.clear()
                NewList.append(record)
                # test1 = record[3]
                # test2 = DataList[i+1][3]
            SeNm = record[0]
            CSeNm = DataList[r + 1][0]
            VGene = record[1]
            CVgene = DataList[r+1][1]
            JGene = record[2]
            CJgene = DataList[r + 1][2]
            CDR3leng = record[3]
            CCDR3leng = DataList[r + 1][3]



            if record[1] == DataList[r+1][1] and record[2] == DataList[r+1][2] and record[3] == DataList[r+1][3]: #has same V, J and CDR3 length and CDR3 is >2
                StartNew = True
                NewList.append(DataList[r+1])
            else:
                if len(NewList) > 1:

                    PrePools.append(tuple(NewList))
                    NewList.clear()
                StartNew = False

        r += 1
        # else:
        #     r += 1

    if len(NewList) > 1:
        PrePools.append(tuple(NewList))

    SeqDict = []
    i =0
    MutFreqDict = {}

    ClonalPools = []

    # 'SeqName', 'VLocus', 'JLocus', 'CDR3Length', 'CDR3DNA', 'Mutations', 'Vbeg', 'Vend', 'Sequence', 'ID'
    # 0             1        2            3           4          5           6       7         8        9
    progress = 1
    for pool in PrePools:
        label = "Processing " + str(poolID) + " group: " + str(progress) + '/' + str(len(PrePools))
        pct = int(progress / len(PrePools) * 100)
        signal.emit(pct, label)
        progress += 1

        SeqDict.clear()
        NearPool = []
        for record in pool:
            Seqname = record[0]
            Sequence = record[8]
            NewList.clear()


            # CDR3 = Sequence[int(record[11]):int(record[12])]
            CDR3 = record[4]


            NewList.append(Seqname)
            NewList.append(Sequence)

            SeqDict.append(tuple(NewList))

            Mutations = record[5]
            MutList = Mutations.split(',')
            Vend = int(record[7])
            Vbeg = int(record[6])
            VLen = (Vend+1)-Vbeg
            Vmuts  = 0
            for mutation in MutList:
                mut = str(mutation)
                MutParse = mut.split('-')
                try:
                    if int(MutParse[1]) < (Vend + 1):
                        Vmuts += 1
                except:
                    Vmuts = Vmuts

                if VLen > 0:
                    VMutFreq = Vmuts/ VLen
                else:
                    VMutFreq = 0
                    print('Problem with Vlen = 0 for ' + Seqname)



            # VMuts = 0
            # for i in range(0, len(CDRSeq)-1):
            #     try:
            #         if CDRSeq[i] != GCDRSeq[i]:
            #             CDRMuts +=1
            #     except:
            #         print('miss')
            # if CDRMuts > 0:
            #     CDRMutFreq =  CDRMuts/len(CDRSeq)#+1
            # else:
            #     CDRMutFreq = 0.0
            # # the +1 is for seq errors
            # Mutations = record[13]
            # MutList = Mutations.split(',')
            # Vend = int(record[14])

            MutFreqDict[Seqname] = (VMutFreq, CDR3, MutList, Vend)

        outfilename = VGenesSeq.ClustalO(SeqDict, 1000, False)

        # outfilename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes','ClustalOmega', 'my-out-seqs.fa')
        Aligned = VGenesSeq.readClustalOutput(outfilename)
        os.remove(outfilename)

        Similars = {}

        SimPool = []
        Matches = []
        Different = []

        i = 0
        LastSeqName = ''

        for seq in range(0,len(Aligned)-1):
            try:
                SeqName = Aligned[seq][0]
                if SeqName == '':
                    print('bad')

                VMutFreq = MutFreqDict[SeqName][0]
                CDR3 = MutFreqDict[SeqName][1]
                CDR3dif = VMutFreq * len(CDR3)
            except:
                CDR3dif = 0

            CCDR3 = ''
            SimPool.clear()
            SimPool.append(SeqName)
            # DupPool.clear()
            # DupPool.append(SeqName)
            for comp in range(seq+1, len(Aligned)):
                try:
                    CSeqName = Aligned[comp][0] #comparison SeqName, etc
                    CVMutFreq = MutFreqDict[CSeqName][0]
                    CCDR3 = MutFreqDict[CSeqName][1]
                except:
                    CVMutFreq = 0
                if Duplicates == False:
                    CCDR3dif = CVMutFreq * len(CCDR3)
                    SeventyFivePercent = len(CDR3)- int(0.75*len(CDR3))

                    differences = 0
                    AllowedDif = CDR3dif + 1 #+ CCDR3dif
                    # if AllowedDif < 1 and CDR3dif >0 and CCDR3dif >0: AllowedDif = 1
                    # AllowedDif+=1 # might need 2...but this is eror rate
                    AllowedDif = ceil(AllowedDif)

                    for j in range(0, len(CDR3)-1):
                        try:
                            if CDR3[j] != CCDR3[j]:
                                differences += 1
                        except:
                            print('stop')
                else:

                    try:
                        if CDR3 == CCDR3:

                            differences = 1

                        else:
                            differences = 0
                    except:
                        differences = 0


                if Duplicates == False:
                    try:
                        if differences <= AllowedDif:
                            # SimPool.append(SeqName)
                            SimPool.append(CSeqName)
                        elif differences <= SeventyFivePercent:  # if 75% nucleotides are shared in CDR3 and also share 1% base exchanges (3/~300), then counted as clonal

                            MutList = MutFreqDict[SeqName][2]
                            # NumMuts = len(MutList)
                            # Vend = MutFreqDict[SeqName][3]
                            MutListComp = MutFreqDict[CSeqName][2]
                            # NumMutsComp = len(MutListComp)
                            # VendComp = MutFreqDict[CSeqName][3]
                            Matches.clear()
                            Matches = [element for element in MutList if element in MutListComp]
                            if len(Matches)>3:  #need get average shared + 2SD from clones for this instead of arbitrary 3
                                SimPool.append(CSeqName)
                        else:

                            if differences == 0:
                                MutList = MutFreqDict[SeqName][2]
                                MutListComp = MutFreqDict[CSeqName][2]
                                if MutList == MutListComp:
                                    SimPool.append(CSeqName)
                    except:
                        SimPool.append(CSeqName)





                # i+=1
            if len(SimPool) > 1:

                # because aligned now can check if this one has everything in previous
                # already and only add to Similars dict if something differs
                # only need to find when stuff is different from list that already exists
                if len(Similars) >0:
                    SetLast = list(Similars[LastSeqName]) #set(Similars[LastSeqName])
                    Matches.clear()
                    # Different.clear()
                    Matches = [element for element in SimPool if element in SetLast]
                    if len(Matches) > 0: #if some match
                        SetLast += [element for element in SimPool if element not in SetLast]
                        Entry = tuple(SetLast)
                        Similars[LastSeqName] = Entry  #update

                    else:

                        Entry = tuple(SimPool)  #Start  new one
                        Similars[SeqName] = Entry
                        LastSeqName = SeqName

                elif len(Similars) == 0:
                    Entry = tuple(SimPool)
                    Similars[SeqName] = Entry
                    LastSeqName = SeqName


        Different.clear()

        Pools = list(Similars.keys())
        # i=0
        KeepList = []
        if len(Similars) >1:

            for i in range(0,len(Pools)-1):
                SimPoolName = Pools[i]
                SimPool = list(Similars[SimPoolName])
                if len(SimPool) != 0:
                    for j in range(i+1,len(Pools)):
                        NextName = Pools[j]
                        nextset = list(Similars[NextName])
                        Matches = [element for element in SimPool if element in nextset]
                        if len(Matches) > 0:
                            Different = [element for element in nextset if element not in SimPool]

                            SimPool = SimPool + Different
                            Similars[SimPoolName] = tuple(SimPool)
                            Similars[NextName] = ''

            # could marge  test name and make rest merged empty so no further matches
            # in end iterate through and make list of tuples of all final merged lists
            # Different.append(list(set))

        NewPools.clear()
        for i in range(0,len(Pools)):  #code to remove empty, merged pools and convert from dictionary to lists
            SimPoolName = Pools[i]
            SimPool = list(Similars[SimPoolName])
            if len(SimPool) != 0:
                NewPools.append(SimPool)


        for pool in NewPools:  #removes duplicates from list, faster then finding records and remarking as clonal later
            NearPool = pool
            NearPool.sort()
            Lennear = len(NearPool)-1
            for k in range (0, Lennear):
                try:
                    if k < (len(NearPool)-1):
                        if k > 0 and NearPool[k] == NearPool[k-1]:
                            NearPool.remove(NearPool[k])
                except:
                    print('Stop')

            if len(NearPool) >1: ClonalPools.append(tuple(NearPool))  #after dups removed if still pool add to final list of clonal pools



    return ClonalPools

# calculates CDR mutation frequency and uses that as the cutoff for
# mutation frequency allowed in the CDR3
#todo need to calculat standard deviation as well and allow within
# it then orders sequences by clustalO and scores CDR3 differences
# in order to see if different from CDR frequnecies

