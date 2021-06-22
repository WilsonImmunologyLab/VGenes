# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'hclctabledialog.ui'
#
# Created by: PyQt5 UI code generator 5.13.0
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_HCLCDialog(object):
    def setupUi(self, HCLCDialog):
        HCLCDialog.setObjectName("HCLCDialog")
        HCLCDialog.resize(808, 920)
        self.gridLayout = QtWidgets.QGridLayout(HCLCDialog)
        self.gridLayout.setObjectName("gridLayout")
        self.tableWidget = QtWidgets.QTableWidget(HCLCDialog)
        self.tableWidget.setObjectName("tableWidget")
        self.tableWidget.setColumnCount(0)
        self.tableWidget.setRowCount(0)
        self.gridLayout.addWidget(self.tableWidget, 1, 0, 1, 1)
        self.label = QtWidgets.QLabel(HCLCDialog)
        font = QtGui.QFont()
        font.setPointSize(16)
        self.label.setFont(font)
        self.label.setText("")
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)

        self.retranslateUi(HCLCDialog)
        QtCore.QMetaObject.connectSlotsByName(HCLCDialog)

    def retranslateUi(self, HCLCDialog):
        _translate = QtCore.QCoreApplication.translate
        HCLCDialog.setWindowTitle(_translate("HCLCDialog", "Dialog"))
