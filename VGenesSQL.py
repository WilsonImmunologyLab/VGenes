__author__ = 'wilsonp'
import sqlite3 as db
import os
# first need connect to a database
from VGenesDialogues import openFile, openFiles, newFile, questionMessage, setText

def creatnewDB(DBpathname):

    (dirname, filename) = os.path.split(DBpathname)

    os.chdir(dirname)
    # print(dirname)


    conn = db.connect(DBpathname)
    cursor = conn.cursor()
    cursor.execute('drop table if exists vgenesdb')

    cursor.execute("create table vgenesDB(SeqName text, SeqLen, GeneType text, V1 text, V2 text, V3 text, D1 text, D2 text, D3 text, J1 text, J2 text, J3 text, StopCodon text, ReadingFrame text, productive text, Strand text, VSeqend text, VDJunction text, Dregion text, DJJunction text, begJ text, VJunction text, FR1From text, FR1To text, FR1length text, FR1matches text, FR1mis text, FR1gaps text, FR1PercentIdentity text, CDR1From text, CDR1to text, CDR1length text, CDR1matches text, CDR1mis text, CDR1gaps text, CDR1PercentIdentity text, FR2From text, FR2To text, FR2length text, FR2matches text, FR2mis text, FR2gaps text, FR2PercentIdentity text, CDR2From text, CDR2to text, CDR2length text, CDR2matches text, CDR2mis text, CDR2gaps text, CDR2PercentIdentity text, FR3From text, FR3To text, FR3length text, FR3matches text, FR3mis text, FR3gaps text, FR3PercentIdentity text, TotMut text, SeqAlignment text, GVbeg text, GVend text, GD1beg text, GD1end text, GD2beg text, GD2end text, GJbeg text, GJend text, Vbeg text, Vend text, D1beg text, D1end text, D2beg text, D2end text, Jbeg text, Jend text, Project text, Grouping text, SubGroup text, Species text, Sequence text, GermlineSequence text, CDR3DNA text, CDR3AA text, CDR3Length text, CDR3beg text, CDR3end text, Specificity text, Subspecificity text, ClonalPool text, ClonalRank text, VLocus text, JLocus text, DLocus text, DateEntered text, Comments text, Quality text, TotalMuts text, Mutations text, IDEvent text, CDR3MW, CDR3pI, Isotype, GCDR3beg, GCDR3end, Blank6, Blank7, Blank8, Blank9, Blank10, Blank11, Blank12, Blank13, Blank14, Blank15, Blank16, Blank17, Blank18, Blank19, Blank20, ID PRIMARY KEY NOT NULL)")
    # ID PRIMARY KEY,

    conn.commit()
    conn.close()



def CopyDatatoDB2(SQLSELECT, DBpathname, DB2path):
    # (dirname, filename) = os.path.split(DBpathname)
    conn = db.connect(DBpathname)
    cursor = conn.cursor()


    # ATTACH DATABASE "\mydir\data\beta.sqlite\" AS beta;
    # CREATE TABLE NewTable AS
    # SELECT * FROM beta.table3;
    # DETACH DATABASE beta;

    # INSERT INTO blog_posts
    # SELECT * FROM BlogProduction.dbo.blog_posts
    SQLStatement = 'ATTACH DATABASE "'+ DB2path + '" AS "DB2"'
    SQLStatement2 = 'INSERT INTO DB2.vgenesDB '+ SQLSELECT #SELECT * FROM vgenesDB WHERE ...' #'INSERT INTO + ' "' + vgenesDB.DB2 ' + '" '+ SQLSELECT #SELECT * FROM vgenesDB WHERE ...'
    SQLStatement3 = 'DETACH DATABASE DB2'
    try:
        cursor.execute(SQLStatement)
        cursor.execute(SQLStatement2)
        cursor.execute(SQLStatement3)
    except:

        print(SQLStatement + ', '+ SQLStatement2 + ', '+ SQLStatement3)


def UpdateMulti(SQLCommand, DBpathname):
    (dirname, filename) = os.path.split(DBpathname)

    os.chdir(dirname)
    # print(dirname)


    conn = db.connect(DBpathname)
    cursor = conn.cursor()
    # cursor.execute('drop table if exists vgenesdb')

    # SQLCommand = 'UPDATE vgenesDB SET ' + Field + ' = "' + Value + '" WHERE ID = ' + ID
    # if Field == 'SeqAlignment':
    #     SQLCommand = 'UPDATE vgenesdb SET Isotype = "' + Value + '" WHERE ID = ' + ID
    try:
        cursor.execute(SQLCommand)
    except:
        print(SQLCommand)
    # SQLCommand = 'UPDATE vgenesdb SET SeqName = "DeletionSeq2", V1 = "333" WHERE ID = 7'

    # ID PRIMARY KEY,

    conn.commit()
    conn.close
def UpdateField(ID, Value, Field, DBpathname):
    (dirname, filename) = os.path.split(DBpathname)

    os.chdir(dirname)
    # print(dirname)


    conn = db.connect(DBpathname)
    cursor = conn.cursor()
    # cursor.execute('drop table if exists vgenesdb')
    nID = str(ID)
    SQLCommand = 'UPDATE vgenesDB SET ' + Field + ' = "' + Value + '" WHERE ID = ' + nID
    # if Field == 'SeqAlignment':
    #     SQLCommand = 'UPDATE vgenesdb SET Isotype = "' + Value + '" WHERE ID = ' + ID
    try:
        cursor.execute(SQLCommand)
    except:
        print(SQLCommand)
    # SQLCommand = 'UPDATE vgenesdb SET SeqName = "DeletionSeq2", V1 = "333" WHERE ID = 7'

    # ID PRIMARY KEY,

    conn.commit()
    conn.close

def UpdateFieldbySeqName(ID, Value, Field, DBpathname):
    (dirname, filename) = os.path.split(DBpathname)

    os.chdir(dirname)
    # print(dirname)


    conn = db.connect(DBpathname)
    cursor = conn.cursor()
    # cursor.execute('drop table if exists vgenesdb')
    nID = str(ID)
    SQLCommand = 'UPDATE vgenesDB SET ' + Field + ' = "' + Value + '" WHERE SeqName = ' + '"' + nID + '"'
    # if Field == 'SeqAlignment':
    #     SQLCommand = 'UPDATE vgenesdb SET Isotype = "' + Value + '" WHERE ID = ' + ID
    try:
        cursor.execute(SQLCommand)
    except:
        print(SQLCommand)
    # SQLCommand = 'UPDATE vgenesdb SET SeqName = "DeletionSeq2", V1 = "333" WHERE ID = 7'

    # ID PRIMARY KEY,

    conn.commit()
    conn.close

def CreateAnalysisDB(FileName, DBpathname):
    # DBpathname = FileName   #os.path.join(os.path.expanduser('~'), 'Dropbox', 'VGenes', 'VDJGenes.db')
    currentFile = ProcessFASTA(FileName)

    (dirname, filename) = os.path.split(DBpathname)
    os.chdir(dirname)

    HumanV = ''
    HumanD = ''
    HumanJ = ''
    MouseV = ''
    MouseD = ''
    MouseJ = ''

    FieldValues = []
    FinalList = []
    Lines = []
    Lines = currentFile.split('\n')
    SeqNum = 0
    for FASTAline in Lines:
        FASTAline = FASTAline.replace('\n', '').replace('\r', '')
        if FASTAline != '' and FASTAline != None:
            itemNum = 0
            if FASTAline[0] == '>':
                # print(FASTAline)
                TitleLine  = FASTAline[1:]
                TitleLine = TitleLine.strip()
                Fields = TitleLine.split('|')
                  # have to remove last field as just from para mark
                for item in Fields:
                    if item == '': item = ' '
                    itemNum += 1
                    if itemNum == 2:
                        VlocusS = item
                        Vlocus = VlocusS[4:(len(VlocusS)-3)]
                        if VlocusS[:4] == 'IGHV':
                            Vlocus = 'VH' + Vlocus
                        elif VlocusS[:4] == 'IGHD':
                            Vlocus = 'DH' + Vlocus
                        elif VlocusS[:4] == 'IGHJ':
                            Vlocus = 'JH' + Vlocus
                        elif VlocusS[:4] == 'IGKV':
                            Vlocus = 'VK' + Vlocus
                        elif VlocusS[:4] == 'IGKJ':
                            Vlocus = 'JK' + Vlocus
                        elif VlocusS[:4] == 'IGLV':
                            Vlocus = 'VL' + Vlocus
                        elif VlocusS[:4] == 'IGLJ':
                            Vlocus = 'JL' + Vlocus
                        else:
                            Vlocus = VlocusS

                        # FieldValues.append(Vlocus)


                    elif itemNum == 3:
                        Strain = item
                        if Strain == "Homo sapiens":
                            # FieldValues.append('Human')
                            SeqName = 'Human' + Vlocus
                        elif Strain == 'Mus musculus':
                            # FieldValues.append('Mouse')
                            SeqName = 'Mouse' + Vlocus
                        else:
                            # FieldValues.append(item)
                            SeqName = Strain + Vlocus

                    # elif itemNum < 5:    # have to remove last field as just from para mark
                FieldValues.append(SeqName) #now put strain specifics



            else:
                ISequence  = FASTAline
                # print(Sequence)
                # FieldValues.append(ISequence)  #IMGT spaced
                # 'SeqName', 'Sequence', 'GermlineSequence', 'FR1To', 'CDR1To','FR2To', 'CDR2To','FR3To', 'CDR3end', 'Jend', 'IDEvent', 'Mutations', 'Vbeg', 'GVbeg', 'GJend'
                # FWR1 = ISequence[:78]
                # # FWR1 = FWR1.replace('.', '')
                # FR1To  = len(FWR1)
                #
                # CDR1 = ISequence[78:114]
                # # CDR1 = CDR1.replace('.', '')
                # CDR1To = len(CDR1)+FR1To
                #
                #
                # FWR2 = ISequence[114:165]
                # # FWR2 = FWR2.replace('.', '')
                # FWR2To = len(FWR2)+CDR1To
                #
                # CDR2 = ISequence[165:195]
                # # CDR2 = CDR2.replace('.', '')
                # CDR2To = len(CDR2)+FWR2To
                #
                FWR3 = ISequence[195:]
                # FWR3 = FWR3.replace('.', '')
                FWR3To = len(ISequence)-1

                # ISequence = ISequence.replace('.', '')
                CDR3end  = len(ISequence)-1

                FieldValues.append(ISequence)
                FieldValues.append(ISequence)

                # FieldValues.append(FR1To)
                # FieldValues.append(CDR1To)
                # FieldValues.append(FWR2To)
                # FieldValues.append(CDR2To)

                FieldValues.append(78)
                FieldValues.append(114)
                FieldValues.append(165)
                FieldValues.append(195)
                FieldValues.append(FWR3To)
                FieldValues.append(CDR3end)
                FieldValues.append(CDR3end)
                FieldValues.append('')
                FieldValues.append('')
                FieldValues.append(1)
                FieldValues.append(1)
                FieldValues.append(CDR3end)

                AAAlen = len(FieldValues)
                if AAAlen != 15:
                    print(FieldValues[1])
                # else:
                FinalList.append(tuple(FieldValues))
                SeqNum += 1
                FieldValues = []


# 'SeqName', 'Sequence', 'GermlineSequence', 'FR1To', 'CDR1To','FR2To', 'CDR2To','FR3To', 'CDR3end', 'Jend', 'IDEvent', 'Mutations', 'Vbeg', 'GVbeg', 'GJend'
        #
    # 'SeqName', 'Sequence', 'GermlineSequence', 'FR1To', 'CDR1To','FR2To', 'CDR2To','FR3To', 'CDR3end', 'Jend', 'IDEvent', 'Mutations', 'Vbeg', 'GVbeg', 'GJend'
    conn = db.connect(DBpathname)
    cursor = conn.cursor()
    cursor.execute('drop table if exists AnalysisDB')

    # >X60503|IGHV1-18*02|Homo sapiens|F|VH-REGION|142..417|276 nt|1|9 | | | |276+24=300|partial in 3'| |
    cursor.execute("create table AnalysisDB(SeqName text, Sequence text, GermlineSequence text, FR1To text, CDR1To text, FR2To text,  CDR2To text,  FR3To text,  CDR3end text,  Jend text, IDEvent text, Mutations text, Vbeg text, GVbeg text, GJend text)")

    if len(FinalList) > 0:
        cursor.executemany('''INSERT INTO AnalysisDB(SeqName, Sequence, GermlineSequence, FR1To, CDR1To, FR2To,  CDR2To,  FR3To,  CDR3end,  Jend, IDEvent, Mutations, Vbeg, GVbeg, GJend) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''', FinalList)


    conn.commit()
    conn.close()

def RunAnalysisSQL(DBpathname):
    # returns a dictionary with seqname as key and all other fileds  as a list as data
    # Note: always needs SeqName to be first field SQLed
    import os

    (dirname, filename) = os.path.split(DBpathname)

    os.chdir(dirname)


    conn = db.connect(DBpathname)
    #  then need to create a cursor that lets you traverse the database
    cursor = conn.cursor()


    SQLStatement = 'SELECT * FROM AnalysisDB'
    # else:
    cursor.execute(SQLStatement)

    DataIs = []
    KeyName  = ''
    Fields = []
    # rows = cursor.rowcount
    for row in cursor:
        # i = 0
        Fields.clear()
        for column in row:
            # if i == 0:
            #     KeyName = column
            # else:
            Fields.append(column)
            # i +=1

        DataIs.append(tuple(Fields))


    # conn.commit()  #  saves data into file
    conn.close()

    return DataIs

def CreateVDJDB(FASTAfilename):
    import os
    DBpathname = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'VDJGenes.db')
    currentFile = ProcessFASTA(FASTAfilename)

    (dirname, filename) = os.path.split(DBpathname)

    HumanV = ''
    HumanD = ''
    HumanJ = ''
    MouseV = ''
    MouseD = ''
    MouseJ = ''

    FieldValues = []
    FinalList = []
    Lines = []
    Lines = currentFile.split('\n')
    SeqNum = 0
    for FASTAline in Lines:
        FASTAline = FASTAline.replace('\n', '').replace('\r', '')
        if FASTAline != '' and FASTAline != None:
            itemNum = 0
            if FASTAline[0] == '>':
                # print(FASTAline)
                TitleLine  = FASTAline[1:]
                TitleLine = TitleLine.strip()
                Fields = TitleLine.split('|')
                  # have to remove last field as just from para mark
                for item in Fields:
                    if item == '': item = ' '
                    itemNum += 1
                    if itemNum == 2:
                        VlocusS = item
                        Vlocus = VlocusS[4:(len(VlocusS)-3)]
                        if VlocusS[:4] == 'IGHV':
                            Vlocus = 'VH' + Vlocus
                        elif VlocusS[:4] == 'IGHD':
                            Vlocus = 'DH' + Vlocus
                        elif VlocusS[:4] == 'IGHJ':
                            Vlocus = 'JH' + Vlocus
                        elif VlocusS[:4] == 'IGKV':
                            Vlocus = 'VK' + Vlocus
                        elif VlocusS[:4] == 'IGKJ':
                            Vlocus = 'JK' + Vlocus
                        elif VlocusS[:4] == 'IGLV':
                            Vlocus = 'VL' + Vlocus
                        elif VlocusS[:4] == 'IGLJ':
                            Vlocus = 'JL' + Vlocus

                        FieldValues.append(Vlocus)


                    elif itemNum == 3:
                        Strain = item
                        if Strain == "Homo sapiens":
                            FieldValues.append('Human')
                        else:
                            FieldValues.append('Mouse')
                    if itemNum < 16:    # have to remove last field as just from para mark
                        FieldValues.append(item) #now put strain specifics



            else:
                ISequence  = FASTAline
                # print(Sequence)
                FieldValues.append(ISequence)  #IMGT spaced
                Sequence = ISequence.replace('.', '')
                if Vlocus[0] == 'V':
                    if Strain == "Homo sapiens":
                        HumanV += '>' + VlocusS + '\n' + Sequence + '\n'
                    else:
                        MouseV += '>' + VlocusS + '\n' + Sequence + '\n'
                if Vlocus[0] == 'D':
                    if Strain == "Homo sapiens":
                        HumanD += '>' + VlocusS + '\n' + Sequence + '\n'
                    else:
                        MouseD += '>' + VlocusS + '\n' + Sequence + '\n'
                if Vlocus[0] == 'J':
                    if Strain == "Homo sapiens":
                        HumanJ += '>' + VlocusS + '\n' + Sequence + '\n'
                    else:
                        MouseJ += '>' + VlocusS + '\n' + Sequence + '\n'

                FieldValues.append(Sequence)  #Direct seq
                FieldValues.append(SeqNum)
                AAAlen = len(FieldValues)
                if AAAlen != 21:
                    print(FieldValues[1])
                # else:
                FinalList.append(tuple(FieldValues))
                SeqNum += 1
                FieldValues = []

    filename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'Database', 'HumanVGenes.nt')

    with open(filename, 'w') as currentfile:
        currentfile.write(HumanV)
    filename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'Database', 'HumanDGenes.nt')

    with open(filename, 'w') as currentfile:
        currentfile.write(HumanD)

    filename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'Database', 'HumanJGenes.nt')

    with open(filename, 'w') as currentfile:
        currentfile.write(HumanJ)

    filename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'Database', 'MouseVGenes.nt')

    with open(filename, 'w') as currentfile:
        currentfile.write(MouseV)

    filename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'Database', 'MouseDGenes.nt')

    with open(filename, 'w') as currentfile:
        currentfile.write(MouseD)

    filename = os.path.join(os.path.expanduser('~'), 'Applications', 'VGenes', 'Database', 'MouseJGenes.nt')

    with open(filename, 'w') as currentfile:
        currentfile.write(MouseJ)

    os.chdir(dirname)
    # print(dirname)




# 1. IMGT/LIGM-DB accession number(s)
# 2. gene and allele name
# 3. species
# 4. functionality
# 5. exon(s), region name(s), or extracted label(s)
# 6. start and end positions in the IMGT/LIGM-DB accession number(s)
# 7. number of nucleotides in the IMGT/LIGM-DB accession number(s)
# 8. codon start, or 'NR' (not relevant) for non coding labels and out-of-frame pseudogenes
# 9. +n: number of nucleotides (nt) added in 5' compared to the corresponding label extracted from IMGT/LIGM-DB
# 10. +n or -n: number of nucleotides (nt) added or removed in 3' compared to the corresponding label extracted from IMGT/LIGM-DB
# 11. +n, -n, and/or nS: number of added, deleted, and/or substituted nucleotides to correct sequencing errors, or 'not corrected' if non corrected sequencing errors
# 12. number of amino acids (AA): this field indicates that the sequence is in amino acids
# 13. number of characters in the sequence: nt (or AA)+IMGT gaps=total
# 14. partial (if it is)
# 15. reverse complementary (if it is)



    conn = db.connect(DBpathname)
    cursor = conn.cursor()
    cursor.execute('drop table if exists GermLineDB')

    # >X60503|IGHV1-18*02|Homo sapiens|F|VH-REGION|142..417|276 nt|1|9 | | | |276+24=300|partial in 3'| |
    cursor.execute("create table GermLineDB(IMGT_ID text, SeqName text, Allele text, Species text, Strain text, functionality text, GeneType text, IMGTDBnumber text, Length text, CodingStartNucleotide text, PlusN3IMGT text, PlusN5IMGT text, IMGTCorrections text, NumberAAifAA text, SeqNucsandIMGTgapsTotal text, PartialSeq text, RevComp text, IMGTSequence text, Sequence text, ID PRIMARY KEY NOT NULL)")

    if len(FinalList) > 0:
        cursor.executemany('''INSERT INTO GermLineDB(IMGT_ID, SeqName, Allele, Species, Strain, functionality, GeneType, IMGTDBnumber, Length, CodingStartNucleotide, PlusN3IMGT, PlusN5IMGT, IMGTCorrections, NumberAAifAA, SeqNucsandIMGTgapsTotal, PartialSeq, RevComp, IMGTSequence, Sequence, ID) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''', FinalList)


    conn.commit()
    conn.close()

def ProcessFASTA(FASTAfile):
    ErLog = ''
    ErlogFile = ''
    import os
    # FASTAfile = os.path.join(os.path.expanduser('~'), 'Dropbox', 'VGenes', 'Database', 'ALLIMGT.nt')
    if FASTAfile == '' or FASTAfile == None:
        return

    with open(FASTAfile, 'r') as currentFile:  #using with for this automatically closes the file even if you crash

        Newline = ''
        Titleline = ''
        ErLog = ''
        CleanSeq = ''
        for FASTAline in currentFile:

            FASTAline = FASTAline.replace('\n', '').replace('\r', '')
            if FASTAline == '':
                FASTAline = ' '

            if FASTAline[0] == '>':
                if len(Newline) > 1:
                    Titleline += '\n'
                    Newline += '\n'
                    CleanSeq += Titleline + Newline
                else:
                    if Newline != ' ' or  Newline != '':
                        if Titleline != '':
                            ErLog += Titleline + ': Sequence error\n'

                Titleline = FASTAline
                Newline = ''
            else:
                FASTAline = FASTAline.upper()
                for nuc in FASTAline:
                    if nuc == 'N' or nuc == 'A' or nuc == 'T' or nuc == 'G' or nuc == 'C' or nuc == '.':
                        Newline += nuc
        # need to write the last sequence into the FASTA file
        if len(Newline) > 1:
            Titleline += '\n'
            Newline += '\n'
            CleanSeq += Titleline + Newline
        else:
            if Newline != ' ' or  Newline != '':
                ErLog += Titleline + '\n'
    return CleanSeq

def OpenDB(DBFilename):
    # todo Need code to open DB and bind to form
    print(DBFilename)

# def updateData(self, DBFilename, IgBLASTAnalysis):
#     (dirname, filename) = os.path.split(DBFilename)
#
#     os.chdir(dirname)
#
#
#     conn = db.connect(DBFilename)
#     #  then need to create a cursor that lets you traverse the database
#     cursor = conn.cursor()
#     # todo must remember when adding fields to change all of this to fit
#     TopNum = 0
#     cursor.execute('''SELECT max(ID) FROM vgenesdb''')  #code to add unique primary keys
#     row = cursor.fetchone()
#     if row[0] == None:
#         TopNum = 0
#     else:
#         TopNum = int(row[0])
# #to get a unique ID
#     # need IgBLASTAnalysis to be list of lists for above operation because tupels are immuntable
#     # then convert to a list of tuples for the expandall to work
#     numberprocessed = 0
#     FinalBLASTed = []
#     for item in IgBLASTAnalysis:
#
#         # todo add code to check each seqname and verify not in the currently open db before adding
#         uId = item[0]
#         # queryis = 'SELECT SeqName FROM vgenesdb WHERE SeqName = ' + seqnamed
#         cursor.execute("SELECT SeqName FROM vgenesdb WHERE SeqName=:Id", {"Id": uId})
#
#         # cursor.execute(queryis)  #code to check if seqname is already in db
#         row = cursor.fetchone()
#
#         exists = True
#         if row == None:
#             exists = False
#         else:
#
#             query = 'A seqence named ' + row[0] + ' is already in this database.\n Add anyways?'
#             answer = questionMessage(self, query, 'YNC')
#             if answer == 'Yes':
#                 exists = False
#             elif answer == 'No':
#                 exists = True
#             elif answer == 'Cancel':
#                 break
#                 # conn.commit()  #  saves data into file
#                 # conn.close()
#                 # return numberprocessed
#
#
#         if exists == False:
#             TopNum += 1
#             item.append(TopNum)     #append a unique seqID
#
#             FinalBLASTed.append(tuple(item))
#             numberprocessed +=1
#
#     if len(FinalBLASTed) > 0:
#         cursor.executemany('''INSERT INTO vgenesDB(SeqName, SeqLen, GeneType, V1, V2, V3, D1, D2, D3, J1, J2, J3, StopCodon, ReadingFrame, productive, Strand, VSeqend, VDJunction, Dregion, DJJunction, begJ, VJunction, FR1From, FR1To, FR1length, FR1matches, FR1mis, FR1gaps, FR1PercentIdentity, CDR1From, CDR1To, CDR1length, CDR1matches, CDR1mis, CDR1gaps, CDR1PercentIdentity, FR2From, FR2To, FR2length, FR2matches, FR2mis, FR2gaps, FR2PercentIdentity, CDR2From, CDR2To, CDR2length, CDR2matches, CDR2mis, CDR2gaps, CDR2PercentIdentity, FR3From, FR3To, FR3length, FR3matches, FR3mis, FR3gaps, FR3PercentIdentity, TotMut, SeqAlignment, GVbeg, GVend, GD1beg, GD1end, GD2beg, GD2end, GJbeg, GJend, Vbeg, Vend, D1beg, D1end, D2beg, D2end, Jbeg, Jend, Project, Grouping, SubGroup, Species, Sequence, GermlineSequence, CDR3DNA, CDR3AA, CDR3Length, CDR3beg, CDR3end, Specificity, Subspecificity, ClonalPool, ClonalRank, VLocus, JLocus, DLocus, DateEntered, Comments, Quality, TotalMuts, Mutations, IDEvent, CDR3MW, CDR3pI, Isotype, GCDR3beg, GCDR3end, Blank6, Blank7, Blank8, Blank9, Blank10, Blank11, Blank12, Blank13, Blank14, Blank15, Blank16, Blank17, Blank18, Blank19, Blank20, ID) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''', FinalBLASTed)
#
#
#     # cursor.execute('''SELECT * FROM vgenesDB''')
#     # rows = cursor.rowcount
#     # for row in cursor:
#     #     for column in row:
#     #         print(column)
#
#     conn.commit()  #  saves data into file
#     conn.close()
#     # readData(DBpathname)
#
#     return numberprocessed
def enterData(self, DBpathname, IgBLASTAnalysis, answer3):

    (dirname, filename) = os.path.split(DBpathname)

    os.chdir(dirname)


    conn = db.connect(DBpathname)
    #  then need to create a cursor that lets you traverse the database
    cursor = conn.cursor()
    # todo must remember when adding fields to change all of this to fit
    TopNum = 0
    cursor.execute('''SELECT max(ID) FROM vgenesDB''')  #code to add unique primary keys
    row = cursor.fetchone()
    if row[0] == None:
        TopNum = 0
    else:
        TopNum = int(row[0])
#to get a unique ID
    # need IgBLASTAnalysis to be list of lists for above operation because tupels are immuntable
    # then convert to a list of tuples for the expandall to work
    numberprocessed = 0
    FinalBLASTed = []
    ErlogFile = '/Applications/IgBlast/database/ErLog.txt'
    ErLog = ''
    Recordlen = 0
    ToAllAnswer = 'none'
    if answer3 == "YesAll":
        ToAllAnswer = 'Yes'
    elif answer3 == 'NoAll':
        ToAllAnswer = 'No'

    for item in IgBLASTAnalysis:

        # todo add code to check each seqname and verify not in the currently open db before adding
        uId = item[0]
        # queryis = 'SELECT SeqName FROM vgenesdb WHERE SeqName = ' + seqnamed
        cursor.execute("SELECT SeqName, Project, Grouping, SubGroup FROM vgenesDB WHERE SeqName=:Id", {"Id": uId})

        # cursor.execute(queryis)  #code to check if seqname is already in db
        row = cursor.fetchone()

        exists = True
        if row == None:
            exists = False
        else:

            query = 'Duplicated sequences exist. Repair by appending project name (Y)? or consider each instance(N)?'
            answering = questionMessage(self, query, 'YN')
            if answering == 'Yes':
                    Project  =  item[75]
                    SeqName  = item[0]

                    NewerName  =  SeqName + Project
                    item[0] = NewerName


            if answering == 'No':


                if ToAllAnswer == 'none':

                    query = 'A sequence named ' + row[0] + ' is already in this database.\n Add anyways?'
                    answer = questionMessage(self, query, 'YNCA')
                    if answer == 'Yes':

                        query = 'Save duplicated sequence with a new name?'
                        answer2 = questionMessage(self, query, 'YN')
                        if answer2 == "Yes":
                            query = 'Existing name: ' + uId + '\nEnter a new name:'
                            DefaultTxt = uId + '_duplicate'
                            newName = setText(self, query, DefaultTxt)
                            item[0] = newName
                            exists = False


                    elif answer == 'No':

                        exists = True


                    elif answer == 'Cancel':
                        break
                    elif answer == 'YesToAll':
                        query = 'Save duplicated sequences with new names (appended with "_duplicated")?'
                        answer3 = questionMessage(self, query, 'YN')
                        if answer3 == 'Yes':
                            answer3 = 'YesAll'

                        exists = False

                        ToAllAnswer = 'Yes'

                    elif answer == 'NoToAll':
                        answer3 = 'NoAll'
                        exists = True
                        ToAllAnswer = 'No'


                else:
                    if ToAllAnswer == 'Yes' or 'YesAll':
                        exists = False
                    elif ToAllAnswer == 'No' or 'NoAll':
                        exists = True

                # conn.commit()  #  saves data into file
                # conn.close()
                # return numberprocessed


        if exists == False:
            Recordlen = len(item)
            if answer3 == 'No' or answer3  == 'NoAll':

                if Recordlen == 119:
                    TopNum += 1
                    item.append(TopNum)     #append a unique seqID

                    FinalBLASTed.append(tuple(item))
                    numberprocessed +=1

            elif answer3 == 'Yes' or answer3  == 'YesAll':

                if Recordlen == 119:
                    TopNum += 1
                    newName = item[0] + '_duplicate'
                    item[0] = newName

                    item.append(TopNum)     #append a unique seqID

                    FinalBLASTed.append(tuple(item))
                    numberprocessed +=1

            else:
                ErLog = uId + ' was problematic (line 399 VGenesSQL: record had '+ str(len(item)) + ' fields rather than 120)\n'
                with open(ErlogFile, 'a') as currentfile:
                    currentfile.write(ErLog)


    if len(FinalBLASTed) > 0:
        cursor.executemany('''INSERT INTO vgenesDB(SeqName, SeqLen, GeneType, V1, V2, V3, D1, D2, D3, J1, J2, J3, StopCodon, ReadingFrame, productive, Strand, VSeqend, VDJunction, Dregion, DJJunction, begJ, VJunction, FR1From, FR1To, FR1length, FR1matches, FR1mis, FR1gaps, FR1PercentIdentity, CDR1From, CDR1to, CDR1length, CDR1matches, CDR1mis, CDR1gaps, CDR1PercentIdentity, FR2From, FR2To, FR2length, FR2matches, FR2mis, FR2gaps, FR2PercentIdentity, CDR2From, CDR2to, CDR2length, CDR2matches, CDR2mis, CDR2gaps, CDR2PercentIdentity, FR3From, FR3To, FR3length, FR3matches, FR3mis, FR3gaps, FR3PercentIdentity, TotMut, SeqAlignment, GVbeg, GVend, GD1beg, GD1end, GD2beg, GD2end, GJbeg, GJend, Vbeg, Vend, D1beg, D1end, D2beg, D2end, Jbeg, Jend, Project, Grouping, SubGroup, Species, Sequence, GermlineSequence, CDR3DNA, CDR3AA, CDR3Length, CDR3beg, CDR3end, Specificity, Subspecificity, ClonalPool, ClonalRank, VLocus, JLocus, DLocus, DateEntered, Comments, Quality, TotalMuts, Mutations, IDEvent, CDR3MW, CDR3pI, Isotype, GCDR3beg, GCDR3end, Blank6, Blank7, Blank8, Blank9, Blank10, Blank11, Blank12, Blank13, Blank14, Blank15, Blank16, Blank17, Blank18, Blank19, Blank20, ID) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''', FinalBLASTed)


    # cursor.execute('''SELECT * FROM vgenesDB''')
    # rows = cursor.rowcount
    # for row in cursor:
    #     for column in row:
    #         print(column)

    conn.commit()  #  saves data into file
    conn.close()
    # readData(DBpathname)

    return numberprocessed, answer3

def ColName(DBpathname):
    import os

    (dirname, filename) = os.path.split(DBpathname)
    os.chdir(dirname)
    conn = db.connect(DBpathname)
    cursor = conn.cursor()

    cursor.execute("SELECT * FROM {}".format('vgenesDB'))
    col_name_list = [tuple[0] for tuple in cursor.description]

    return col_name_list


def RunSQL(DBpathname, SQLStatement):
    # returns a dictionary with seqname as key and all other fileds  as a list as data
    # Note: always needs SeqName to be first field SQLed
    import os

    (dirname, filename) = os.path.split(DBpathname)

    os.chdir(dirname)


    conn = db.connect(DBpathname)
    #  then need to create a cursor that lets you traverse the database
    cursor = conn.cursor()

    if SQLStatement == 'None':
        cursor.execute('''SELECT * FROM vgenesDB''')
    else:
        cursor.execute(SQLStatement)

    DataIs = []
    KeyName  = ''
    Fields = []
    # rows = cursor.rowcount
    for row in cursor:
        # i = 0
        Fields.clear()
        for column in row:
            # if i == 0:
            #     KeyName = column
            # else:
            Fields.append(column)
            # i +=1

        DataIs.append(tuple(Fields))


    # conn.commit()  #  saves data into file
    conn.close()

    return DataIs


def ImportVDB(pathname, DBFilename):
    (dirname, filename) = os.path.split(DBFilename)

    os.chdir(dirname)


    conn = db.connect(DBFilename)
    #  then need to create a cursor that lets you traverse the database
    cursor = conn.cursor()
    pathname+= '.vgenesDB'
    cursor.execute('''ATTACH pathname as VGDB2''')
    cursor.execute('''SELECT * FROM VGDB2''')
    DataIs = []
    for row in cursor:
        for column in row:
            DataIs.append(column)


    cursor.execute('''INSERT INTO DBpathname.vgenesDB SELECT * FROM VGDB2.vgenesDB''')

    conn.close()

def readData(DBpathname, SQLStatement):
    import os

    (dirname, filename) = os.path.split(DBpathname)

    os.chdir(dirname)


    conn = db.connect(DBpathname)
    #  then need to create a cursor that lets you traverse the database
    cursor = conn.cursor()

    if SQLStatement == 'None':
        cursor.execute('''SELECT * FROM vgenesDB''')
    else:
        cursor.execute(SQLStatement)

    DataIs = []
    # rows = cursor.rowcount
    for row in cursor:
        for column in row:
            DataIs.append(column)

    # conn.commit()  #  saves data into file
    conn.close()

    return DataIs

def FetchOneRecord(databasename):

    conn = db.connect(databasename)
    #  then need to create a cursor that lets you traverse the database
    cursor = conn.cursor()
    #  then need to create a cursor that lets you traverse the database

    # #can use sql to determine averages
    cursor.execute('select avg(temp) from temps')
    row = cursor.fetchone()

    conn.close()
    print('The average temp for the week was: %s' % row[0])

def deleterecords (DBFilename, SQLStatement):
    import os


    (dirname, filename) = os.path.split(DBFilename)

    os.chdir(dirname)


    conn = db.connect(DBFilename)
    #  then need to create a cursor that lets you traverse the database
    cursor = conn.cursor()
    SQLStatement2 = 'DELETE' + SQLStatement[18:]
    try:
        if SQLStatement2 != 'None':
            cursor.execute(SQLStatement2)
    except:
        SQLStatement2 = 'DELETE FROM vgenesDB WHERE Project = "Delete"'
        cursor.execute(SQLStatement2)



    # CommandSql = 'delete from vgenesDB where temp = ' + searchvalue
    # cursor.execute('delete from temps where temp = 40')

    # # need to requery and fetch data to display with changes
    # cursor.execute('select * from temps')
    # rows = cursor.fetchall()
    # for row in rows:
    #     print('%s %s' % (row[0], row[1]))
    conn.commit()
    conn.close()





def MakeSQLStatement(self, fields, SeqName):


    checkedProjects, checkedGroups, checkedSubGroups, checkedkids = self.getTreeChecked()

    SQLStatement = 'SELECT '

    if fields != 'All':
        fieldCount = len(fields)
        i = 1
        for field in fields:
            SQLStatement += field
            if i < fieldCount:
                SQLStatement += ', '
            else:
                SQLStatement += ' FROM vgenesDB'
            i += 1
    else:
        SQLStatement += '* FROM vgenesDB'  # 'SELECT * FROM vgenesDB WHERE ID = '

    firstmore = False

    if (len(checkedProjects) + len(checkedGroups) + len(checkedSubGroups) + len(
            checkedkids)) > 0:  # then something is seleected
        SQLStatement += ' WHERE '
        firstmore = True
    else:
        SQLStatement += ' WHERE SeqName = "'
        SQLStatement += SeqName + '"'


    i = 1
    if len(checkedProjects) > 0:
        if firstmore == True:
            firstmore = False
        else:
            SQLStatement += ', '
            firstmore = False
        project = self.ui.cboTreeOp1.currentText()
        fieldname = self.TransLateFieldtoReal(project, True)
        SQLStatement = SQLStatement + fieldname + ' = '
        for item in checkedProjects:
            SQLStatement += ('"' + item)
            if i < len(checkedProjects):
                SQLStatement += '" OR ' + fieldname + ' = '
            else:
                SQLStatement += '"'
            i += 1

            # if len(checkedGroups) > 0: SQLStatement += ', OR '

    i = 1
    if len(checkedGroups) > 0:
        if firstmore == True:
            firstmore = False
        else:
            if len(checkedProjects) > 0:
                SQLStatement += ' OR '
            else:
                SQLStatement += ', '

            firstmore = False
        group = self.ui.cboTreeOp2.currentText()
        fieldname = self.TransLateFieldtoReal(group, True)
        project = self.ui.cboTreeOp1.currentText()
        Projfieldname = self.TransLateFieldtoReal(project, True)

        # SQLStatement = SQLStatement + fieldname + ' = "'
        for item in checkedGroups:
            statement = fieldname + ' = "' + item[1] + '" AND ' + Projfieldname + ' = "' + item[0]
            SQLStatement += statement
            if i < len(checkedGroups):
                SQLStatement += '" OR '
            else:
                SQLStatement += '"'
            i += 1

    i = 1
    if len(checkedSubGroups) > 0:
        if firstmore == True:
            firstmore = False
        else:
            if len(checkedProjects) > 0 or len(checkedGroups) > 0:
                SQLStatement += ' OR '
            else:
                SQLStatement += ', '
            firstmore = False
        Subgroup = self.ui.cboTreeOp3.currentText()
        fieldname = self.TransLateFieldtoReal(Subgroup, True)
        group = self.ui.cboTreeOp2.currentText()
        Groupfieldname = self.TransLateFieldtoReal(group, True)
        project = self.ui.cboTreeOp1.currentText()
        Projfieldname = self.TransLateFieldtoReal(project, True)

        # SQLStatement = SQLStatement + fieldname + ' = "'
        for item in checkedSubGroups:
            statement = fieldname + ' = "' + item[2] + '" AND ' + Groupfieldname + ' = "' + item[
                1] + '" AND ' + Projfieldname + ' = "' + item[0]
            SQLStatement += statement
            if i < len(checkedSubGroups):
                SQLStatement += '" OR '
            else:
                SQLStatement += '"'
            i += 1

    i = 1
    if len(checkedkids) > 0:
        if firstmore == True:
            firstmore = False
        else:
            if len(checkedProjects) > 0 or len(checkedGroups) > 0 or len(checkedSubGroups) > 0:
                SQLStatement += ' OR '
            else:
                SQLStatement += ', '
            firstmore = False

        for item in checkedkids:
            SQLStatement += 'SeqName = "'
            SQLStatement += item
            if i < len(checkedkids):
                SQLStatement += '" OR '
            else:
                SQLStatement += '"'
            i += 1

    return SQLStatement

