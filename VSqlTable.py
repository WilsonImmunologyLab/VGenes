__author__ = 'wilsonp'
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication, QTableView, QTreeView
from PyQt5.QtSql import QSqlTableModel

import vgenesconnection  #for some reason this only works with PyQt-SQL not python built in



def initializeModel(model):
    model.setTable('vgenesdb')

    model.setEditStrategy(QSqlTableModel.OnManualSubmit)
    model.select()

    # model.setHeaderData(0, Qt.Horizontal, "Index")
    # model.setHeaderData(1, Qt.Horizontal, "Name")
    # model.setHeaderData(2, Qt.Horizontal, "DNA")
    # model.setHeaderData(3, Qt.Horizontal, "Protein")


def createView(title, model):
    # view = QTreeView()
    view = QTableView()
    view.setModel(model)
    view.setWindowTitle(title)
    return view



if __name__ == '__main__':

    import sys

    app = QApplication(sys.argv)
    if not vgenesconnection.createConnection():
        sys.exit(1)

    # record = "insert into vgenesdb values(101, 'TestSeq', 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNTCCCTGAGACTCTCCTGTGAAGCCTCTGGATTCACGTTTGGCGGCAATGGCATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGACTGGAGTGGGTCGCAGGTATCAGTGGTATTAGTGGGAATACATATTATTTAGGCTCCGTGAAGGGCCGGTTCACCATCTCCAGAGACAATCCGAAGAGGACGTTATATCTACAAATGAATCGTCTGAGAGTCGAGGACACGGCCATTTATTACTGTGCGAAAGATCGTTTGATAGGAACAGATGGCGTCTCTTTCTTTGACCAATGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG', 'test')"
    #
    # vgenesconnection.AddRecord(record)

    model = QSqlTableModel()

    initializeModel(model)

    view1 = createView("Table Model (View 1)", model)
    # view2 = createView("Table Model (View 2)", model)


    view1.show()
    # view2.move(view1.x() + view1.width() + 20, view1.y())
    # view2.show()

    sys.exit(app.exec_())