# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'table_dialog.ui'
#
# Created by: PyQt5 UI code generator 5.13.0
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_TableDialog(object):
    def setupUi(self, TableDialog):
        TableDialog.setObjectName("TableDialog")
        TableDialog.resize(781, 825)
        self.gridLayout_2 = QtWidgets.QGridLayout(TableDialog)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.checkBoxAll = QtWidgets.QCheckBox(TableDialog)
        self.checkBoxAll.setObjectName("checkBoxAll")
        self.gridLayout_2.addWidget(self.checkBoxAll, 1, 0, 1, 1)
        self.pushButtonSave = QtWidgets.QPushButton(TableDialog)
        self.pushButtonSave.setDefault(True)
        self.pushButtonSave.setObjectName("pushButtonSave")
        self.gridLayout_2.addWidget(self.pushButtonSave, 3, 0, 1, 1)
        self.frame = QtWidgets.QFrame(TableDialog)
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.gridLayout = QtWidgets.QGridLayout(self.frame)
        self.gridLayout.setObjectName("gridLayout")
        self.label = QtWidgets.QLabel(self.frame)
        font = QtGui.QFont()
        font.setPointSize(18)
        self.label.setFont(font)
        self.label.setAlignment(QtCore.Qt.AlignCenter)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.gridLayout_2.addWidget(self.frame, 0, 0, 1, 4)
        spacerItem = QtWidgets.QSpacerItem(465, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_2.addItem(spacerItem, 3, 2, 1, 1)
        self.tableWidget = QtWidgets.QTableWidget(TableDialog)
        self.tableWidget.setObjectName("tableWidget")
        self.tableWidget.setColumnCount(0)
        self.tableWidget.setRowCount(0)
        self.gridLayout_2.addWidget(self.tableWidget, 2, 0, 1, 4)
        self.pushButtonClose = QtWidgets.QPushButton(TableDialog)
        self.pushButtonClose.setObjectName("pushButtonClose")
        self.gridLayout_2.addWidget(self.pushButtonClose, 3, 3, 1, 1)
        self.checkBoxRow = QtWidgets.QCheckBox(TableDialog)
        self.checkBoxRow.setObjectName("checkBoxRow")
        self.gridLayout_2.addWidget(self.checkBoxRow, 1, 1, 1, 1)

        self.retranslateUi(TableDialog)
        QtCore.QMetaObject.connectSlotsByName(TableDialog)

    def retranslateUi(self, TableDialog):
        _translate = QtCore.QCoreApplication.translate
        TableDialog.setWindowTitle(_translate("TableDialog", "Dialog"))
        self.checkBoxAll.setText(_translate("TableDialog", "Check all"))
        self.pushButtonSave.setText(_translate("TableDialog", "Save"))
        self.label.setText(_translate("TableDialog", "Fields of data table"))
        self.pushButtonClose.setText(_translate("TableDialog", "Close"))
        self.checkBoxRow.setText(_translate("TableDialog", "Select entire row"))
