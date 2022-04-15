__author__ = 'wilsonp'
import sqlite3 as db
import os
import re
# first need connect to a database
from VGenesDialogues import openFile, openFiles, newFile, questionMessage, setText

def creatnewDB(DBpathname):

    (dirname, filename) = os.path.split(DBpathname)

    os.chdir(dirname)
    # print(dirname)

    conn = db.connect(DBpathname)
    cursor = conn.cursor()
    cursor.execute('drop table if exists vgenesdb')
    cursor.execute("create table vgenesDB(SeqName text UNIQUE, SeqLen, GeneType text, V1 text, V2 text, V3 text, D1 text, D2 text, D3 text, J1 text, J2 text, J3 text, StopCodon text, ReadingFrame text, productive text, Strand text, VSeqend text, VDJunction text, Dregion text, DJJunction text, begJ text, VJunction text, FR1From text, FR1To text, FR1length text, FR1matches text, FR1mis text, FR1gaps text, FR1PercentIdentity text, CDR1From text, CDR1to text, CDR1length text, CDR1matches text, CDR1mis text, CDR1gaps text, CDR1PercentIdentity text, FR2From text, FR2To text, FR2length text, FR2matches text, FR2mis text, FR2gaps text, FR2PercentIdentity text, CDR2From text, CDR2to text, CDR2length text, CDR2matches text, CDR2mis text, CDR2gaps text, CDR2PercentIdentity text, FR3From text, FR3To text, FR3length text, FR3matches text, FR3mis text, FR3gaps text, FR3PercentIdentity text, TotMut text, SeqAlignment text, GVbeg text, GVend text, GD1beg text, GD1end text, GD2beg text, GD2end text, GJbeg text, GJend text, Vbeg text, Vend text, D1beg text, D1end text, D2beg text, D2end text, Jbeg text, Jend text, Project text, Grouping text, SubGroup text, Species text, Sequence text, GermlineSequence text, CDR3DNA text, CDR3AA text, CDR3Length text, CDR3beg text, CDR3end text, Specificity text, Subspecificity text, ClonalPool text, ClonalRank text, VLocus text, JLocus text, DLocus text, DateEntered text, Comments text, Quality text, TotalMuts text, Mutations text, IDEvent text, CDR3MW, CDR3pI, Isotype, GCDR3beg, GCDR3end, Blank6, Blank7, Blank8, Blank9, Blank10, Blank11, Blank12, Blank13, Blank14, Blank15, Blank16, Blank17, Blank18, Blank19, Blank20, ID PRIMARY KEY NOT NULL)")
    # ID PRIMARY KEY,
    cursor.execute('drop table if exists fieldsname')
    cursor.execute("CREATE TABLE IF NOT EXISTS fieldsname(ID int PRIMARY KEY NOT NULL, Field text, FieldNickName text, FieldType text, FieldComment text, display text, display_priority text)")

    FieldList = [
        [1, "SeqName", "Name", "Fixed", ""],
        [2, "SeqLen", "Length", "Fixed", ""],
        [3, "GeneType", "Type", "Fixed", ""],
        [4, "V1", "Vgene", "Fixed", ""],
        [5, "V2", "Vgene2ndchoice", "Fixed", ""],
        [6, "V3", "Vgene3rdchoice", "Fixed", ""],
        [7, "D1", "Dgene", "Fixed", ""],
        [8, "D2", "Dgene2ndchoice", "Fixed", ""],
        [9, "D3", "Dgene3rdchoice", "Fixed", ""],
        [10, "J1", "Jgene", "Fixed", ""],
        [11, "J2", "Jgene2ndchoice", "Fixed", ""],
        [12, "J3", "Jgene3rdchoice", "Fixed", ""],
        [13, "StopCodon", "Stopcodons?", "Fixed", ""],
        [14, "ReadingFrame", "Readingframe", "Fixed", ""],
        [15, "productive", "Productive?", "Fixed", ""],
        [16, "Strand", "Strand", "Fixed", ""],
        [17, "VSeqend", "EndofVgene", "Fixed", ""],
        [18, "VDJunction", "VtoDJunction", "Fixed", ""],
        [19, "Dregion", "Dregion", "Fixed", ""],
        [20, "DJJunction", "DtoJjunction", "Fixed", ""],
        [21, "begJ", "BeginningofJ", "Fixed", ""],
        [22, "VJunction", "VtoJjunction", "Fixed", ""],
        [23, "FR1From", "FWR1firstbase", "Fixed", ""],
        [24, "FR1To", "FWR1lastbase", "Fixed", ""],
        [25, "FR1length", "FWR1length", "Fixed", ""],
        [26, "FR1matches", "FWR1matches", "Fixed", ""],
        [27, "FR1mis", "FWR1mismatches", "Fixed", ""],
        [28, "FR1gaps", "FWR1gaps", "Fixed", ""],
        [29, "FR1PercentIdentity", "FWR1percentidentity", "Fixed", ""],
        [30, "CDR1From", "CDR1firstbase", "Fixed", ""],
        [31, "CDR1to", "CDR1lastbase", "Fixed", ""],
        [32, "CDR1length", "CDR1length", "Fixed", ""],
        [33, "CDR1matches", "CDR1matches", "Fixed", ""],
        [34, "CDR1mis", "CDR1mismatches", "Fixed", ""],
        [35, "CDR1gaps", "CDR1gaps", "Fixed", ""],
        [36, "CDR1PercentIdentity", "CDR1percentidentity", "Fixed", ""],
        [37, "FR2From", "FWR2firstbase", "Fixed", ""],
        [38, "FR2To", "FWR2lastbase", "Fixed", ""],
        [39, "FR2length", "FWR2length", "Fixed", ""],
        [40, "FR2matches", "FWR2matches", "Fixed", ""],
        [41, "FR2mis", "FWR2mismatches", "Fixed", ""],
        [42, "FR2gaps", "FWR2gaps", "Fixed", ""],
        [43, "FR2PercentIdentity", "FWR2percentidentity", "Fixed", ""],
        [44, "CDR2From", "CDR2firstbase", "Fixed", ""],
        [45, "CDR2to", "CDR2lastbase", "Fixed", ""],
        [46, "CDR2length", "CDR2length", "Fixed", ""],
        [47, "CDR2matches", "CDR2matches", "Fixed", ""],
        [48, "CDR2mis", "CDR2mismatches", "Fixed", ""],
        [49, "CDR2gaps", "CDR2gaps", "Fixed", ""],
        [50, "CDR2PercentIdentity", "CDR2percentidentity", "Fixed", ""],
        [51, "FR3From", "FWR3firstbase", "Fixed", ""],
        [52, "FR3To", "FWR3lastbase", "Fixed", ""],
        [53, "FR3length", "FWR3length", "Fixed", ""],
        [54, "FR3matches", "FWR3matches", "Fixed", ""],
        [55, "FR3mis", "FWR3mismatches", "Fixed", ""],
        [56, "FR3gaps", "FWR3gaps", "Fixed", ""],
        [57, "FR3PercentIdentity", "FWR3percentidentity", "Fixed", ""],
        [58, "TotMut", "IgBLASTmutationcount(V)", "Fixed", ""],
        [59, "SeqAlignment", "IgBLASTSequenceAlignment", "Fixed", ""],
        [60, "GVbeg", "GermlineVbegin", "Fixed", ""],
        [61, "GVend", "GermlineVend", "Fixed", ""],
        [62, "GD1beg", "GermlineD1begin", "Fixed", ""],
        [63, "GD1end", "GermlineD1end", "Fixed", ""],
        [64, "GD2beg", "GermlineD2begin", "Fixed", ""],
        [65, "GD2end", "GermlineD2end", "Fixed", ""],
        [66, "GJbeg", "GermlineJbegin", "Fixed", ""],
        [67, "GJend", "GermlineJend", "Fixed", ""],
        [68, "Vbeg", "Vbegin", "Fixed", ""],
        [69, "Vend", "Vend", "Fixed", ""],
        [70, "D1beg", "D1begin", "Fixed", ""],
        [71, "D1end", "D1end", "Fixed", ""],
        [72, "D2beg", "D2begin", "Fixed", ""],
        [73, "D2end", "D2end", "Fixed", ""],
        [74, "Jbeg", "Jbegin", "Fixed", ""],
        [75, "Jend", "Jend", "Fixed", ""],
        [76, "Project", "Project", "Fixed", ""],
        [77, "Grouping", "Grouping", "Fixed", ""],
        [78, "SubGroup", "Subgroup", "Fixed", ""],
        [79, "Species", "Species", "Fixed", ""],
        [80, "Sequence", "Sequence", "Fixed", ""],
        [81, "GermlineSequence", "Germlinesequence", "Fixed", ""],
        [82, "CDR3DNA", "CDR3DNA", "Fixed", ""],
        [83, "CDR3AA", "CDR3peptide", "Fixed", ""],
        [84, "CDR3Length", "CDR3length", "Fixed", ""],
        [85, "CDR3beg", "CDR3firstbase", "Fixed", ""],
        [86, "CDR3end", "CDR3lastbase", "Fixed", ""],
        [87, "Specificity", "Specificity", "Fixed", ""],
        [88, "Subspecificity", "Subspecificity", "Fixed", ""],
        [89, "ClonalPool", "ClonalPool", "Fixed", ""],
        [90, "ClonalRank", "ClonalRank", "Fixed", ""],
        [91, "VLocus", "Vlocus", "Fixed", ""],
        [92, "JLocus", "Jlocus", "Fixed", ""],
        [93, "DLocus", "Dlocus", "Fixed", ""],
        [94, "DateEntered", "Dateandtimeentered", "Fixed", ""],
        [95, "Comments", "Comments", "Fixed", ""],
        [96, "Quality", "Quality", "Fixed", ""],
        [97, "TotalMuts", "TotalMutations(VDJ)", "Fixed", ""],
        [98, "Mutations", "Mutationlist", "Fixed", ""],
        [99, "IDEvent", "Insertions&deletions", "Fixed", ""],
        [100, "CDR3MW", "CDR3molecularweight", "Fixed", ""],
        [101, "CDR3pI", "CDR3isoelectricpoint", "Fixed", ""],
        [102, "Isotype", "Isotype", "Fixed", ""],
        [103, "GCDR3beg", "GermlneCDR3begin", "Fixed", ""],
        [104, "GCDR3end", "GermlineCDR3end", "Fixed", ""],
        [105, "Blank6", "Autoreactivity", "Fixed", ""],
        [106, "Blank7", "ORF", "Fixed", ""],
        [107, "Blank8", "10xCluster", "Fixed", ""],
        [108, "Blank9", "Seuret_Cluster", "Fixed", ""],
        [109, "Blank10", "10xBarCode", "Fixed", ""],
        [110, "Blank11", "Population", "Fixed", ""],
        [111, "Blank12", "Label", "Fixed", ""],
        [112, "Blank13", "Status", "Fixed", ""],
        [113, "Blank14", "Blank14", "Fixed", ""],
        [114, "Blank15", "Blank15", "Fixed", ""],
        [115, "Blank16", "Blank16", "Fixed", ""],
        [116, "Blank17", "Blank17", "Fixed", ""],
        [117, "Blank18", "Blank18", "Fixed", ""],
        [118, "Blank19", "Blank19", "Fixed", ""],
        [119, "Blank20", "RawSeq", "Fixed", ""],
        [120, "ID", "ID", "Fixed", ""]
    ]

    yes_list = [1, 2, 3, 4, 7, 10, 13, 14, 15, 16, 76, 77, 78, 79, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94,
                95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115,
                116, 117, 118, 119, 120]

    for ele in FieldList:
        if ele[0] in yes_list:
            ele.append('yes')
        else:
            ele.append('no')

        if ele[0] == 1:
            ele.append(0)
        else:
            ele.append(9)

    cursor.executemany(
        '''INSERT INTO fieldsname(ID, Field, FieldNickName, FieldType, FieldComment, display, display_priority) VALUES(?,?,?,?,?,?,?)''',
        FieldList)

    conn.commit()
    conn.close()

def checkFieldTable(DBpathname):
    (dirname, filename) = os.path.split(DBpathname)

    os.chdir(dirname)

    conn = db.connect(DBpathname)
    cursor = conn.cursor()

    yes_list = [1, 2, 3, 4, 7, 10, 13, 14, 15, 16, 76, 77, 78, 79, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94,
                95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115,
                116, 117, 118, 119, 120]
    yes_str = [str(i) for i in yes_list]
    yes_str = ','.join(yes_str)

    try:
        cursor.execute("SELECT * FROM fieldsname WHERE ID = 1")
        HEADERStatement = 'PRAGMA table_info(fieldsname);'
        HeaderIn = RunSQL(DBpathname, HEADERStatement)
        header_list = [i[1] for i in HeaderIn]
        if 'display' in header_list:
            pass
        else:
            SQLStatement = "ALTER TABLE fieldsname ADD display text"
            cursor.execute(SQLStatement)
            conn.commit()
            SQLStatement = 'UPDATE fieldsname SET display = "no" WHERE 1'
            cursor.execute(SQLStatement)
            SQLStatement = 'UPDATE fieldsname SET display = "yes" WHERE ID in (' + yes_str + ')'
            cursor.execute(SQLStatement)
            SQLStatement = 'UPDATE fieldsname SET display = "yes" WHERE ID > 110'
            cursor.execute(SQLStatement)
            conn.commit()

        if 'display_priority' in header_list:
            pass
        else:
            SQLStatement = "ALTER TABLE fieldsname ADD display_priority text"
            cursor.execute(SQLStatement)
            conn.commit()
            SQLStatement = 'UPDATE fieldsname SET display_priority = 9 WHERE 1'
            cursor.execute(SQLStatement)
            SQLStatement = 'UPDATE fieldsname SET display_priority = 0 WHERE ID = 1'
            cursor.execute(SQLStatement)
            conn.commit()
    except:
        #cursor.execute('drop table if exists fieldsname')
        cursor.execute("CREATE TABLE IF NOT EXISTS fieldsname(ID int PRIMARY KEY NOT NULL, Field text, FieldNickName text, FieldType text, FieldComment text, display text, display_priority text)")

        HEADERStatement = 'PRAGMA table_info(vgenesDB);'
        HeaderIn = RunSQL(DBpathname, HEADERStatement)
        Fields = [i[1] for i in HeaderIn]

        keys = ['SeqName', 'SeqLen', 'GeneType', 'V1', 'V2', 'V3', 'D1', 'D2', 'D3', 'J1', 'J2', 'J3', 'StopCodon',
            'ReadingFrame', 'productive', 'Strand', 'VSeqend', 'VDJunction', 'Dregion', 'DJJunction', 'begJ',
            'VJunction', 'FR1From', 'FR1To', 'FR1length', 'FR1matches', 'FR1mis', 'FR1gaps', 'FR1PercentIdentity',
            'CDR1From', 'CDR1to', 'CDR1length', 'CDR1matches', 'CDR1mis', 'CDR1gaps', 'CDR1PercentIdentity', 'FR2From',
            'FR2To', 'FR2length', 'FR2matches', 'FR2mis', 'FR2gaps', 'FR2PercentIdentity', 'CDR2From', 'CDR2to',
            'CDR2length', 'CDR2matches', 'CDR2mis', 'CDR2gaps', 'CDR2PercentIdentity', 'FR3From', 'FR3To', 'FR3length',
            'FR3matches', 'FR3mis', 'FR3gaps', 'FR3PercentIdentity', 'TotMut', 'SeqAlignment', 'GVbeg', 'GVend',
            'GD1beg', 'GD1end', 'GD2beg', 'GD2end', 'GJbeg', 'GJend', 'Vbeg', 'Vend', 'D1beg', 'D1end', 'D2beg',
            'D2end', 'Jbeg', 'Jend', 'Project', 'Grouping', 'SubGroup', 'Species', 'Sequence', 'GermlineSequence',
            'CDR3DNA', 'CDR3AA', 'CDR3Length', 'CDR3beg', 'CDR3end', 'Specificity', 'Subspecificity', 'ClonalPool',
            'ClonalRank', 'VLocus', 'JLocus', 'DLocus', 'DateEntered', 'Comments', 'Quality', 'TotalMuts', 'Mutations',
            'IDEvent', 'CDR3MW', 'CDR3pI', 'Isotype', 'GCDR3beg', 'GCDR3end', 'Blank6', 'Blank7', 'Blank8', 'Blank9',
            'Blank10', 'Blank11', 'Blank12', 'Blank13', 'Blank14', 'Blank15', 'Blank16', 'Blank17', 'Blank18',
            'Blank19', 'Blank20', 'ID']
        values = ["Name", "Length", "Type", "V gene", "V gene 2nd choice", "V gene 3rd choice", "D gene",
             " D gene 2nd choice", "D gene 3rd choice", "J gene", " J gene 2nd choice", "J gene 3rd choice",
             "Stop codons?", "Reading frame", "Productive?", "Strand", "End of V gene", "V to D Junction",
             "D region", "D to J junction", "Beginning of J", "V to J junction", "FWR1 first base", "FWR1 last base",
             "FWR1 length", "FWR1 matches", "FWR1 mismatches", "FWR1 gaps", "FWR1 percent identity",
             "CDR1 first base", "CDR1 last base", "CDR1 length", "CDR1 matches", "CDR1 mismatches", "CDR1 gaps",
             "CDR1 percent identity", "FWR2 first base", "FWR2 last base", "FWR2 length", "FWR2 matches",
             "FWR2 mismatches", "FWR2 gaps", "FWR2 percent identity", "CDR2 first base", "CDR2 last base",
             "CDR2 length", "CDR2 matches", "CDR2 mismatches", "CDR2 gaps", "CDR2 percent identity",
             "FWR3 first base", "FWR3 last base", "FWR3 length", "FWR3 matches", "FWR3 mismatches", "FWR3 gaps",
             "FWR3 percent identity", "IgBLAST mutation count", "IgBLAST Sequence Alignment", "Germline V begin",
             "Germline V end", "Germline D1 begin", "Germline D1 end", "Germline D2 begin", "Germline D2 end",
             "Germline J begin", "Germline J end", " V begin", " V end", " D1 begin", " D1 end", " D2 begin",
             " D2 end", " J begin", " J end", "Project", "Grouping", "Subgroup", "Species", "Sequence",
             " Germline sequence", "CDR3 DNA", "CDR3 peptide", "CDR3 length", "CDR3 first base", "CDR3 last base",
             "Specificity", "Subspecificity", "Clonal Pool", "Clonal Rank", "V locus", "J locus", "D locus",
             "Date and time entered", "Comments", "Quality", "Total Mutations", "Mutation list",
             "Insertions & deletions", "CDR3 molecular weight", "CDR3 isoelectric point", "Isotype",
             "Germlne CDR3 begin", "Germline CDR3 end", "Autoreactivity", "Blank7", "10xCluster", "Seuret_Cluster",
             "10xBarCode", "Population",
             "Label", "Status", "Blank14", "Blank15", "Blank16", "Blank17", "Blank18", "Blank19", "Blank20", "ID"]

        Dict = dict(zip(keys, values))

        i = 1
        FieldList = []
        for field in Fields:
            if i < 106:
                if field in keys:
                    ele = [i, field, Dict[field], "Fixed", ""]
                else:
                    ele = [i, field, field, "Fixed", ""]
            else:
                if field in keys:
                    ele = [i, field, Dict[field], "Customized", ""]
                else:
                    ele = [i, field, field, "Customized", ""]

            if i in yes_list:
                ele.append('yes')
            else:
                ele.append('no')

            if i == 1:
                ele.append(0)
            else:
                ele.append(9)
            FieldList.append(ele)
            i += 1

        cursor.executemany('''INSERT INTO fieldsname(ID, Field, FieldNickName, FieldType, FieldComment, display, display_priority) VALUES(?,?,?,?,?,?,?)''',FieldList)
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

def UpdateFieldTable(ID, Value, Field, DBpathname):
    (dirname, filename) = os.path.split(DBpathname)

    os.chdir(dirname)
    # print(dirname)

    conn = db.connect(DBpathname)
    cursor = conn.cursor()
    # cursor.execute('drop table if exists vgenesdb')
    nID = str(ID)
    SQLCommand = 'UPDATE fieldsname SET ' + Field + ' = "' + Value + '" WHERE Field = "' + nID + '"'

    try:
        cursor.execute(SQLCommand)
    except:
        print(SQLCommand)

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

def DumpDB(DBpathname, tmp_path, recordNames):
    import time
    import re

    conn = db.connect(DBpathname)
    time_stamp = time.strftime("%Y-%m-%d-%H_%M_%S", time.localtime())
    file_path = os.path.join(tmp_path, time_stamp + '_dump.sql')

    if len(recordNames) == 0:
        with open(file_path, 'w') as f:
            for line in conn.iterdump():
                f.write('%s\n' % line)
    else:
        with open(file_path, 'w') as f:
            for line in conn.iterdump():
                if line.startswith('INSERT INTO "vgenesDB"'):
                    match = re.findall(r"VALUES\(\'([^\']+)", line)
                    if match[0] in recordNames:
                        f.write('%s\n' % line)
                else:
                    f.write('%s\n' % line)

    conn.close()
    return file_path

def ImportDB(DBpathname, SQLfile):
    try:
        conn = db.connect(DBpathname)
        cur = conn.cursor()
        SQL_f = open(SQLfile, 'r')
        SQL_str = SQL_f.read()
        cur.executescript(SQL_str)
        conn.close()
        return True
    except:
        return False

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
def enterData(self, DBpathname, IgBLASTAnalysis, answer3, ErlogFile):

    (dirname, filename) = os.path.split(DBpathname)
    try:
        os.chdir(dirname)
    except:
        pass


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
    #ErlogFile = os.path.join(temp_folder,'ErLog.txt')
    ErLog = ''
    Recordlen = 0
    ToAllAnswer = 'none'
    if answer3 == "YesAll":
        ToAllAnswer = 'Yes'
    elif answer3 == 'NoAll':
        ToAllAnswer = 'No'
    '''
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

                    FinalBLASTed.append(item)
                    numberprocessed +=1

            elif answer3 == 'Yes' or answer3  == 'YesAll':

                if Recordlen == 119:
                    TopNum += 1
                    newName = item[0] + '_duplicate'
                    item[0] = newName

                    item.append(TopNum)     #append a unique seqID

                    FinalBLASTed.append(item)
                    numberprocessed +=1

            else:
                ErLog = uId + ' was problematic (line 399 VGenesSQL: record had '+ str(len(item)) + ' fields rather than 120)\n'
                with open(ErlogFile, 'a') as currentfile:
                    currentfile.write(ErLog)
    '''
    for item in IgBLASTAnalysis:
        Recordlen = len(item)
        if Recordlen == 119:
            TopNum += 1
            item.append(TopNum)  # append a unique seqID

            FinalBLASTed.append(item)
            numberprocessed += 1

    with open(ErlogFile, 'a') as currentfile:
        currentfile.write("\n")
    dup_message = ''
    if len(FinalBLASTed) > 0:
        cursor.execute("SELECT SeqName FROM vgenesDB")
        rows = cursor.fetchall()
        names = [row[0] for row in rows]
        for item in FinalBLASTed:
            num = 1
            ori_name = item[0]
            cur_name = item[0]
            while cur_name in names:
                cur_name = ori_name + '_Duplicate' + str(num)
                num += 1
            item[0] = cur_name
            if num > 1:
                dup_message = ori_name + "    to    " + item[0] + "\n"
                with open(ErlogFile, 'a') as currentfile:
                    currentfile.write(dup_message)
            names.append(item[0])
        try:
            cursor.executemany('''INSERT INTO vgenesDB(SeqName, SeqLen, GeneType, V1, V2, V3, D1, D2, D3, J1, J2, J3, StopCodon, ReadingFrame, productive, Strand, VSeqend, VDJunction, Dregion, DJJunction, begJ, VJunction, FR1From, FR1To, FR1length, FR1matches, FR1mis, FR1gaps, FR1PercentIdentity, CDR1From, CDR1to, CDR1length, CDR1matches, CDR1mis, CDR1gaps, CDR1PercentIdentity, FR2From, FR2To, FR2length, FR2matches, FR2mis, FR2gaps, FR2PercentIdentity, CDR2From, CDR2to, CDR2length, CDR2matches, CDR2mis, CDR2gaps, CDR2PercentIdentity, FR3From, FR3To, FR3length, FR3matches, FR3mis, FR3gaps, FR3PercentIdentity, TotMut, SeqAlignment, GVbeg, GVend, GD1beg, GD1end, GD2beg, GD2end, GJbeg, GJend, Vbeg, Vend, D1beg, D1end, D2beg, D2end, Jbeg, Jend, Project, Grouping, SubGroup, Species, Sequence, GermlineSequence, CDR3DNA, CDR3AA, CDR3Length, CDR3beg, CDR3end, Specificity, Subspecificity, ClonalPool, ClonalRank, VLocus, JLocus, DLocus, DateEntered, Comments, Quality, TotalMuts, Mutations, IDEvent, CDR3MW, CDR3pI, Isotype, GCDR3beg, GCDR3end, Blank6, Blank7, Blank8, Blank9, Blank10, Blank11, Blank12, Blank13, Blank14, Blank15, Blank16, Blank17, Blank18, Blank19, Blank20, ID) VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)''', FinalBLASTed)
            conn.commit()  # saves data into file
            conn.close()
            if dup_message == "":
                return numberprocessed, answer3, ['n', dup_message]
            else:
                return numberprocessed, answer3, ['d', 'Please see log message for details!']
        except:
            dup_message = 'Importing error!'
            conn.close()
            return numberprocessed, answer3, ['e', dup_message]
    else:
        dup_message = 'No records!'
        conn.close()
        return numberprocessed, answer3, ['e', dup_message]



def ColName(DBpathname):
    import os

    (dirname, filename) = os.path.split(DBpathname)
    os.chdir(dirname)
    conn = db.connect(DBpathname)
    cursor = conn.cursor()

    cursor.execute("SELECT * FROM {}".format('vgenesDB'))
    col_name_list = [tuple[0] for tuple in cursor.description]

    return col_name_list

def RunUpdateSQL(DBpathname, SQLStatement):
    import os

    (dirname, filename) = os.path.split(DBpathname)
    os.chdir(dirname)
    conn = db.connect(DBpathname)
    cursor = conn.cursor()
    res = cursor.execute(SQLStatement)
    conn.commit()
    conn.close()
    return res.rowcount

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
        Fields.clear()
        for column in row:
            Fields.append(column)
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
    SQLStatement2 = 'DELETE' + SQLStatement
    try:
        if SQLStatement2 != 'None':
            cursor.execute(SQLStatement2)
    except:
        print(SQLStatement2)
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
    checkedProjects, checkedGroups, checkedSubGroups, checkedkids = self.getTreeCheckedChild()

    SQLStatement = 'SELECT '
    SQLStatement_all = ''

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
    SQLStatement_all = SQLStatement
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
        #fieldname = self.TransLateFieldtoReal(project, True)
        fieldname = re.sub(r'\(.+', '', project)
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
        #fieldname = self.TransLateFieldtoReal(group, True)
        fieldname = re.sub(r'\(.+', '', group)
        project = self.ui.cboTreeOp1.currentText()
        #Projfieldname = self.TransLateFieldtoReal(project, True)
        Projfieldname = re.sub(r'\(.+', '', project)

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
        #fieldname = self.TransLateFieldtoReal(Subgroup, True)
        fieldname = re.sub(r'\(.+', '', Subgroup)
        group = self.ui.cboTreeOp2.currentText()
        #Groupfieldname = self.TransLateFieldtoReal(group, True)
        Groupfieldname = re.sub(r'\(.+', '', group)
        project = self.ui.cboTreeOp1.currentText()
        #Projfieldname = self.TransLateFieldtoReal(project, True)
        Projfieldname = re.sub(r'\(.+', '', project)

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

        SQLStatement += 'SeqName IN ('
        for item in checkedkids:
            SQLStatement += '"' + item + '",'
            i += 1
        SQLStatement = SQLStatement.rstrip(',')
        SQLStatement += ')'
    else:
        question = 'You did not check any records, process all?'
        buttons = 'YN'
        answer = questionMessage(self, question, buttons)
        if answer == 'Yes':
            SQLStatement = SQLStatement_all

    return SQLStatement

def MakeSQLStatementNew(self, fields, SeqName):
    checkedkids = self.CheckedRecords

    SQLStatement = 'SELECT '
    SQLStatement_all = ''

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
    SQLStatement_all = SQLStatement

    i = 1
    if len(checkedkids) > 0:
        SQLStatement += ' WHERE SeqName IN ('
        for item in checkedkids:
            SQLStatement += '"' + item + '",'
            i += 1
        SQLStatement = SQLStatement.rstrip(',')
        SQLStatement += ')'
    else:
        SQLStatement = SQLStatement_all

    return SQLStatement

