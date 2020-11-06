# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'patentdialog.ui'
#
# Created by: PyQt5 UI code generator 5.13.0
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_PatentDialog(object):
    def setupUi(self, PatentDialog):
        PatentDialog.setObjectName("PatentDialog")
        PatentDialog.resize(1094, 837)
        self.gridLayout_2 = QtWidgets.QGridLayout(PatentDialog)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.tableWidget = QtWidgets.QTableWidget(PatentDialog)
        self.tableWidget.setObjectName("tableWidget")
        self.tableWidget.setColumnCount(0)
        self.tableWidget.setRowCount(0)
        self.gridLayout_2.addWidget(self.tableWidget, 1, 0, 1, 3)
        spacerItem = QtWidgets.QSpacerItem(891, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_2.addItem(spacerItem, 2, 0, 1, 1)
        self.pushButtonCancel = QtWidgets.QPushButton(PatentDialog)
        self.pushButtonCancel.setObjectName("pushButtonCancel")
        self.gridLayout_2.addWidget(self.pushButtonCancel, 2, 1, 1, 1)
        self.pushButtonConfirm = QtWidgets.QPushButton(PatentDialog)
        self.pushButtonConfirm.setDefault(True)
        self.pushButtonConfirm.setObjectName("pushButtonConfirm")
        self.gridLayout_2.addWidget(self.pushButtonConfirm, 2, 2, 1, 1)
        self.frame = QtWidgets.QFrame(PatentDialog)
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.gridLayout = QtWidgets.QGridLayout(self.frame)
        self.gridLayout.setObjectName("gridLayout")
        self.label = QtWidgets.QLabel(self.frame)
        font = QtGui.QFont()
        font.setPointSize(16)
        self.label.setFont(font)
        self.label.setAlignment(QtCore.Qt.AlignCenter)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 0, 0, 1, 1)
        self.gridLayout_2.addWidget(self.frame, 0, 0, 1, 3)

        self.retranslateUi(PatentDialog)
        QtCore.QMetaObject.connectSlotsByName(PatentDialog)

    def retranslateUi(self, PatentDialog):
        _translate = QtCore.QCoreApplication.translate
        PatentDialog.setWindowTitle(_translate("PatentDialog", "Dialog"))
        self.pushButtonCancel.setText(_translate("PatentDialog", "Cancel"))
        self.pushButtonConfirm.setText(_translate("PatentDialog", "Confirm"))
        self.label.setText(_translate("PatentDialog", "Review your paired heavy chians and light chains"))
