__author__ = 'wilsonp'


from PyQt5.QtWidgets import QMessageBox
from PyQt5.QtSql import QSqlDatabase, QSqlQuery


def createConnection():
    db = QSqlDatabase.addDatabase('QSQLITE')
    db.setDatabaseName('poop.db') #to have only in memory make this (':memory:')


    if not db.open():
        QMessageBox.critical(None, "Cannot open database",
                "Unable to establish a database connection.\n"
                "This example needs SQLite support. Please read the Qt SQL "
                "driver documentation for information how to build it.\n\n"
                "Click Cancel to exit.",
                QMessageBox.Cancel)
        return False

    query = QSqlQuery()
    # query.exec_("create table vgenesdb(id int primary key, "
    #             "seqname varchar(20), sequence varchar(500))")

    # query.exec_("create table vgenesdb(id integer, seqname text, sequence text)")
    # query.exec_('drop table if exists vgenesdb')
    # query.exec_("create table vgenesdb(id integer, SeqName text, DNA text, Protein text)")  #, SeqLen integer, GeneType text, Isotype text, Vgene1 text, Vgene2 text, Vgene3 text, Dgene1 text, Dgene2 text, Dgene3 text, Jgene1 text, Jgene2 text, Jgene3 text, StopCodon text, ReadingFrame integer, Productive text, Strand text, VSeqEnd text, VDJunction text, Dregion text, DJJunction text, Jseqbeg text, VJjunction text, FR1from integer, FR1To integer, FR1length integer, FR1matched integer, FR1mis integer, FR1gaps integer, FR1PercentID real, FR2from integer, FR2To integer, FR2length integer, FR2matched integer, FR2mis integer, FR2gaps integer, FR2PercentID real, FR3from integer, FR3To integer, FR3length integer, FR3matched integer, FR3mis integer, FR3gaps integer, FR3PercentID real, CDR1from integer, CDR1To integer, CDR1length integer, CDR1matched integer, CDR1mis integer, CDr1gaps integer, CDR1PercentID real, CDR2from integer, CDR2To integer, CDR2length integer, CDR2matched integer, CDR2mis integer, CDR2gaps integer, CDR2PercentID real, IgBlastAlignMent text, TotMut integer, PercentIdentity real, CDR3DNA text, CDR3AA text, CDR3Length integer, CDR3pI real, CDR3MolWgt real, Specificity text, Species text, Project text, GroupID text, SubGroup text, ClonalPool integer, ClonalRank integer, GermlineSequence text, VLocus text, JLocus text, DLocus text, DateEntered text, Comments text, Quality text, InsDel1 text, InsDel2 text)")
    # query.exec_('delete from vgenesdb where SeqName = TestSeq')

    # query.exec_("insert into vgenesdb values(102, 'fromOther', 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCCCTGAGACTCTCCTGTGAAGCCTCTGGATTCACGTTTGGCGGCAATGGCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGACTGGAGTGGGTCGCAGGTATCAGTGGTATTAGTGGGAATACATATTATTTAGGCTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATCCGAAGAGGACGTTATATCTACAAATGAATCGTCTGAGAGTCGAGGACACGGCCATTTATTACTGTGCGAAAGATCGTTTGATAGGAACAGATGGCGTCTCTTTCTTTGACCAATGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG', 'test')")

    # query.exec_(record)



    return True