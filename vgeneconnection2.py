__author__ = 'wilsonp'
import sqlite3 as db
# first need connect to a database
def createConnection():

    conn = db.connect('poop.db')
    #  then need to create a cursor that lets you traverse the database
    cursor = conn.cursor()
    #  then need to create a cursor that lets you traverse the database


    cursor.execute('drop table if exists vgenesdb')
    cursor.execute("create table vgenesdb(id integer, SeqName text, DNA text, Protein text)")  #, SeqLen integer, GeneType text, Isotype text, Vgene1 text, Vgene2 text, Vgene3 text, Dgene1 text, Dgene2 text, Dgene3 text, Jgene1 text, Jgene2 text, Jgene3 text, StopCodon text, ReadingFrame integer, Productive text, Strand text, VSeqEnd text, VDJunction text, Dregion text, DJJunction text, Jseqbeg text, VJjunction text, FR1from integer, FR1To integer, FR1length integer, FR1matched integer, FR1mis integer, FR1gaps integer, FR1PercentID real, FR2from integer, FR2To integer, FR2length integer, FR2matched integer, FR2mis integer, FR2gaps integer, FR2PercentID real, FR3from integer, FR3To integer, FR3length integer, FR3matched integer, FR3mis integer, FR3gaps integer, FR3PercentID real, CDR1from integer, CDR1To integer, CDR1length integer, CDR1matched integer, CDR1mis integer, CDr1gaps integer, CDR1PercentID real, CDR2from integer, CDR2To integer, CDR2length integer, CDR2matched integer, CDR2mis integer, CDR2gaps integer, CDR2PercentID real, IgBlastAlignMent text, TotMut integer, PercentIdentity real, CDR3DNA text, CDR3AA text, CDR3Length integer, CDR3pI real, CDR3MolWgt real, Specificity text, Species text, Project text, GroupID text, SubGroup text, ClonalPool integer, ClonalRank integer, GermlineSequence text, VLocus text, JLocus text, DLocus text, DateEntered text, Comments text, Quality text, InsDel1 text, InsDel2 text)")

    test = "insert into " + "vgenesdb values(101, 'FromString4', 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCCCTGAGACTCTCCTGTGAAGCCTCTGGATTCACGTTTGGCGGCAATGGCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGACTGGAGTGGGTCGCAGGTATCAGTGGTATTAGTGGGAATACATATTATTTAGGCTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATCCGAAGAGGACGTTATATCTACAAATGAATCGTCTGAGAGTCGAGGACACGGCCATTTATTACTGTGCGAAAGATCGTTTGATAGGAACAGATGGCGTCTCTTTCTTTGACCAATGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG', 'test')"


    cursor.execute("insert into vgenesdb values(101, 'new', 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCCCTGAGACTCTCCTGTGAAGCCTCTGGATTCACGTTTGGCGGCAATGGCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGACTGGAGTGGGTCGCAGGTATCAGTGGTATTAGTGGGAATACATATTATTTAGGCTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATCCGAAGAGGACGTTATATCTACAAATGAATCGTCTGAGAGTCGAGGACACGGCCATTTATTACTGTGCGAAAGATCGTTTGATAGGAACAGATGGCGTCTCTTTCTTTGACCAATGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG', 'test')")
    cursor.execute("insert into vgenesdb values(102, 'TestSeq2', 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCCCTGAGACTCTCCTGTGAAGCCTCTGGATTCACGTTTGGCGGCAATGGCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGACTGGAGTGGGTCGCAGGTATCAGTGGTATTAGTGGGAATACATATTATTTAGGCTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATCCGAAGAGGACGTTATATCTACAAATGAATCGTCTGAGAGTCGAGGACACGGCCATTTATTACTGTGCGAAAGATCGTTTGATAGGAACAGATGGCGTCTCTTTCTTTGACCAATGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG', 'test')")
    cursor.execute(test)


    conn.commit()  #  saves data into file

    conn.row_factory = db.Row
    cursor.execute('select * from vgenesdb')
    rows = cursor.fetchall()
    for row in rows:
        print('%s %s' % (row[0], row[1]))
    conn.close()

def FetchOneRecord(databasename):

    conn = db.connect(databasename)
    #  then need to create a cursor that lets you traverse the database
    cursor = conn.cursor()
    #  then need to create a cursor that lets you traverse the database

    # #can use sql to determine averages
    cursor.execute('select avg(temp) from temps')
    row = cursor.fetchone()


    print('The average temp for the week was: %s' % row[0])

def deleterecord (databasename, fieldname, searchvalue):
    conn = db.connect(databasename)
    #  then need to create a cursor that lets you traverse the database
    cursor = conn.cursor()
    #  then need to create a cursor that lets you traverse the database

    CommandSql = 'delete from ' + fieldname + ' where temp = ' + searchvalue
    cursor.execute('delete from temps where temp = 40')

    # # need to requery and fetch data to display with changes
    # cursor.execute('select * from temps')
    # rows = cursor.fetchall()
    # for row in rows:
    #     print('%s %s' % (row[0], row[1]))




createConnection()