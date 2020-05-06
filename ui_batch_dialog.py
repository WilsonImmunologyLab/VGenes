# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'batch_dialog.ui'
#
# Created by: PyQt5 UI code generator 5.13.0
#
# WARNING! All changes made in this file will be lost!


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_BatchDialog(object):
    def setupUi(self, BatchDialog):
        BatchDialog.setObjectName("BatchDialog")
        BatchDialog.resize(586, 672)
        self.gridLayout_2 = QtWidgets.QGridLayout(BatchDialog)
        self.gridLayout_2.setObjectName("gridLayout_2")
        self.label = QtWidgets.QLabel(BatchDialog)
        font = QtGui.QFont()
        font.setPointSize(15)
        self.label.setFont(font)
        self.label.setObjectName("label")
        self.gridLayout_2.addWidget(self.label, 0, 0, 1, 1)
        self.gridLayoutFig = QtWidgets.QGridLayout()
        self.gridLayoutFig.setObjectName("gridLayoutFig")
        self.gridLayout_2.addLayout(self.gridLayoutFig, 1, 0, 1, 4)
        self.radioButton = QtWidgets.QRadioButton(BatchDialog)
        self.radioButton.setObjectName("radioButton")
        self.gridLayout_2.addWidget(self.radioButton, 2, 0, 1, 2)
        self.scrollArea = QtWidgets.QScrollArea(BatchDialog)
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setObjectName("scrollArea")
        self.scrollAreaWidgetContents = QtWidgets.QWidget()
        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 560, 537))
        self.scrollAreaWidgetContents.setObjectName("scrollAreaWidgetContents")
        self.gridLayout_3 = QtWidgets.QGridLayout(self.scrollAreaWidgetContents)
        self.gridLayout_3.setObjectName("gridLayout_3")
        self.gridLayout = QtWidgets.QGridLayout()
        self.gridLayout.setObjectName("gridLayout")
        self.gridLayout_3.addLayout(self.gridLayout, 1, 0, 1, 1)
        spacerItem = QtWidgets.QSpacerItem(20, 498, QtWidgets.QSizePolicy.Minimum, QtWidgets.QSizePolicy.Expanding)
        self.gridLayout_3.addItem(spacerItem, 2, 0, 1, 1)
        self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        self.gridLayout_2.addWidget(self.scrollArea, 3, 0, 1, 4)
        spacerItem1 = QtWidgets.QSpacerItem(315, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout_2.addItem(spacerItem1, 4, 0, 1, 2)
        self.pushButtonCancel = QtWidgets.QPushButton(BatchDialog)
        self.pushButtonCancel.setMaximumSize(QtCore.QSize(120, 16777215))
        self.pushButtonCancel.setObjectName("pushButtonCancel")
        self.gridLayout_2.addWidget(self.pushButtonCancel, 4, 2, 1, 1)
        self.pushButtonOK = QtWidgets.QPushButton(BatchDialog)
        self.pushButtonOK.setMaximumSize(QtCore.QSize(120, 16777215))
        self.pushButtonOK.setDefault(True)
        self.pushButtonOK.setObjectName("pushButtonOK")
        self.gridLayout_2.addWidget(self.pushButtonOK, 4, 3, 1, 1)
        self.comboBox = QtWidgets.QComboBox(BatchDialog)
        self.comboBox.setObjectName("comboBox")
        self.gridLayout_2.addWidget(self.comboBox, 0, 1, 1, 3)

        self.retranslateUi(BatchDialog)
        QtCore.QMetaObject.connectSlotsByName(BatchDialog)

    def retranslateUi(self, BatchDialog):
        _translate = QtCore.QCoreApplication.translate
        BatchDialog.setWindowTitle(_translate("BatchDialog", "Dialog"))
        self.label.setText(_translate("BatchDialog", "Edit the values of :"))
        self.radioButton.setText(_translate("BatchDialog", "This column is numerical"))
        self.pushButtonCancel.setText(_translate("BatchDialog", "Cancel"))
        self.pushButtonOK.setText(_translate("BatchDialog", "Confirm"))
