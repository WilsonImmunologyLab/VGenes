# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'shmtabledialog.ui'
#
# Created by: PyQt5 UI code generator 5.13.0
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_SHMtableDialog(object):
    def setupUi(self, SHMtableDialog):
        SHMtableDialog.setObjectName("SHMtableDialog")
        SHMtableDialog.resize(819, 617)
        self.gridLayout = QtWidgets.QGridLayout(SHMtableDialog)
        self.gridLayout.setObjectName("gridLayout")
        self.tableWidget = QtWidgets.QTableWidget(SHMtableDialog)
        self.tableWidget.setObjectName("tableWidget")
        self.tableWidget.setColumnCount(0)
        self.tableWidget.setRowCount(0)
        self.gridLayout.addWidget(self.tableWidget, 0, 0, 1, 1)

        self.retranslateUi(SHMtableDialog)
        QtCore.QMetaObject.connectSlotsByName(SHMtableDialog)

    def retranslateUi(self, SHMtableDialog):
        _translate = QtCore.QCoreApplication.translate
        SHMtableDialog.setWindowTitle(_translate("SHMtableDialog", "Dialog"))
