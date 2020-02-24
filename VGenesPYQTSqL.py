__author__ = 'wilsonp'

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QColor
from PyQt5.QtWidgets import QApplication, QTableView
from PyQt5.QtSql import QSqlQuery, QSqlQueryModel


from PyQt5.QtWidgets import QMessageBox
from PyQt5.QtSql import QSqlDatabase, QSqlQuery

global FieldList
# FieldList = []
FieldList = ['SeqName', 'LenSeq', 'GeneType', 'V1', 'V2', 'V3', 'D1', 'D2', 'D3', 'J1', 'J2', 'J3', 'StopCodon', 'ReadingFrame', 'productive', 'Strand', 'VSeqend', 'VDJunction', 'Dregion', 'DJJunction', 'begJ', 'VJunction', 'FR1From', 'FR1To', 'FR1length', 'FR1matches', 'FR1mis', 'FR1gaps', 'FR1PercentIdentity', 'CDR1From', 'CDR1To', 'CDR1length', 'CDR1matches', 'CDR1mis', 'CDR1gaps', 'CDR1PercentIdentity', 'FR2From', 'FR2To', 'FR2length', 'FR2matches', 'FR2mis', 'FR2gaps', 'FR2PercentIdentity', 'CDR2From', 'CDR2To', 'CDR2length', 'CDR2matches', 'CDR2mis', 'CDR2gaps', 'CDR2PercentIdentity', 'FR3From', 'FR3To', 'FR3length', 'FR3matches', 'FR3mis', 'FR3gaps', 'FR3PercentIdentity', 'TotMut', 'SeqAlignment', 'GVbeg', 'GVend', 'GD1beg', 'GD1end', 'GD2beg', 'GD2end', 'GJbeg', 'GJend', 'Vbeg', 'Vend', 'D1beg', 'D1end', 'D2beg', 'D2end', 'Jbeg', 'Jend', 'Project', 'Grouping', 'SubGroup', 'Species', 'Sequence', 'GermlineSequence', 'CDR3DNA', 'CDR3AA', 'CDR3Length', 'CDR3beg', 'CDR3end', 'Specificity', 'Subspecificity', 'ClonalPool', 'ClonalRank', 'VLocus', 'JLocus', 'DLocus', 'DateEntered', 'Comments', 'Quality', 'TotalMuts', 'Mutations',  'IDEvent', 'CDR3MW', 'CDR3pI', 'Blank3', 'Blank4', 'Blank5', 'Blank6', 'Blank7', 'Blank8', 'Blank9', 'Blank10', 'Blank11', 'Blank12', 'Blank13', 'Blank14', 'Blank15', 'Blank16', 'Blank17', 'Blank18', 'Blank19', 'Blank20', 'ID']

# SeqName, LenSeq, GeneType, V1, V2, V3, D1, D2, D3, J1, J2, J3, StopCodon, ReadingFrame, productive, Strand, VSeqend, VDJunction, Dregion, DJJunction, begJ, VJunction, FR1From, FR1To, FR1length, FR1matches, FR1mis, FR1gaps, FR1PercentIdentity, CDR1From, CDR1To, CDR1length, CDR1matches, CDR1mis, CDR1gaps, CDR1PercentIdentity, FR2From, FR2To, FR2length, FR2matches, FR2mis, FR2gaps, FR2PercentIdentity, CDR2From, CDR2To, CDR2length, CDR2matches, CDR2mis, CDR2gaps, CDR2PercentIdentity, FR3From, FR3To, FR3length, FR3matches, FR3mis, FR3gaps, FR3PercentIdentity, TotMut, SeqAlignment, GVbeg, GVend, GD1beg, GD1end, GD2beg, GD2end, GJbeg, GJend, Vbeg, Vend, D1beg, D1end, D2beg, D2end, Jbeg, Jend, Project, Grouping, SubGroup, Species, Sequence, GermlineSequence, CDR3DNA, CDR3AA, CDR3Length, CDR3pI, CDR3MolWgt, Specificity, Subspecificity, ClonalPool, ClonalRank, VLocus, JLocus, DLocus, DateEntered, Comments, Quality, InsDel1, InsDel2, ID

def createConnection(self, DBFilename):


    db = QSqlDatabase.addDatabase('QSQLITE')
    # db.setDatabaseName(':memory:')
    db.setDatabaseName(DBFilename)
    if not db.open():
        QMessageBox.critical(None, "Cannot open database",
                "Unable to establish a database connection.\n",
                QMessageBox.Cancel)
        return False


    return True

class EditableSqlModel(QSqlQueryModel):
    def flags(self, index):
        flags = super(EditableSqlModel, self).flags(index)

        if index.column() in (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118):
            flags |= Qt.ItemIsEditable

        return flags

    def setData(self, index, value, role):
        if index.column() not in (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118):
            return False

        primaryKeyIndex = self.index(index.row(), 119)
        id = self.data(primaryKeyIndex)

        self.clear()

        i = 0
        for item in FieldList:
            if index.column() == i:
                field  = item
                ok = self.setValue(id, value, field)
                self.refresh()
                return
            i += 1


        self.refresh()
        return

    def refresh(self):
        self.setQuery('select * from vgenesdb')
        i = 0
        for item in FieldList:
            self.setHeaderData(i, Qt.Horizontal, item)
            i += 1

    def setValue(self, Id, Newvalue, field):
        query = QSqlQuery()
        queryText = 'update vgenesdb set ' + field + ' = ? where ID = ?'
        query.prepare(queryText)
        query.addBindValue(Newvalue)
        query.addBindValue(Id)
        return query.exec_()



def initializeModel(model):
    model.setQuery('select * from vgenesdb')
    i = 0
    for item in FieldList:
        model.setHeaderData(i, Qt.Horizontal, item)

    model.setHeaderData(0, Qt.Horizontal, "SeqName")

# offset = 0
# views = []


    # view.setWindowTitle(title)
    # view.move(100 + offset, 100 + offset)
    # offset += 20
    # view.show()

# if __name__ == '__main__':
#
#     import sys
#
#     app = QApplication(sys.argv)


    # sys.exit(app.exec_())
