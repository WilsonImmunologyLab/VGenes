# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'VGenesStartUpDialogue.ui'
#
# Created: Thu Sep  4 11:05:44 2014
#      by: PyQt5 UI code generator 5.3.1
#
# WARNING! All changes made in this file will be lost!

from PyQt5 import QtCore, QtGui, QtWidgets

class Ui_VGenesStartUpDialog(object):
    def setupUi(self, VGenesStartUpDialog):
        VGenesStartUpDialog.setObjectName("VGenesStartUpDialog")
        VGenesStartUpDialog.resize(480, 168)
        self.gridLayout = QtWidgets.QGridLayout(VGenesStartUpDialog)
        self.gridLayout.setObjectName("gridLayout")
        self.label_2 = QtWidgets.QLabel(VGenesStartUpDialog)
        self.label_2.setMinimumSize(QtCore.QSize(100, 0))
        self.label_2.setMaximumSize(QtCore.QSize(16777215, 16777215))
        self.label_2.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
        self.label_2.setObjectName("label_2")
        self.gridLayout.addWidget(self.label_2, 0, 0, 1, 3)
        self.label = QtWidgets.QLabel(VGenesStartUpDialog)
        self.label.setMinimumSize(QtCore.QSize(90, 0))
        self.label.setMaximumSize(QtCore.QSize(85, 16777215))
        self.label.setAlignment(QtCore.Qt.AlignRight|QtCore.Qt.AlignTrailing|QtCore.Qt.AlignVCenter)
        self.label.setObjectName("label")
        self.gridLayout.addWidget(self.label, 1, 0, 1, 1)
        self.cboRecent = QtWidgets.QComboBox(VGenesStartUpDialog)
        self.cboRecent.setMinimumSize(QtCore.QSize(300, 30))
        self.cboRecent.setMaxVisibleItems(15)
        self.cboRecent.setObjectName("cboRecent")
        self.cboRecent.addItem("")
        self.gridLayout.addWidget(self.cboRecent, 1, 1, 1, 3)
        spacerItem = QtWidgets.QSpacerItem(159, 17, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem, 2, 0, 1, 2)
        self.btnNew = QtWidgets.QPushButton(VGenesStartUpDialog)
        self.btnNew.setMinimumSize(QtCore.QSize(180, 45))
        self.btnNew.setMaximumSize(QtCore.QSize(180, 16777215))
        icon = QtGui.QIcon()
        icon.addPixmap(QtGui.QPixmap(":/PNG-Icons/blank_page.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.btnNew.setIcon(icon)
        self.btnNew.setIconSize(QtCore.QSize(24, 24))
        self.btnNew.setObjectName("btnNew")
        self.gridLayout.addWidget(self.btnNew, 2, 2, 1, 1)
        spacerItem1 = QtWidgets.QSpacerItem(95, 17, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem1, 2, 3, 1, 1)
        spacerItem2 = QtWidgets.QSpacerItem(200, 20, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem2, 3, 0, 1, 2)
        self.cboOpen = QtWidgets.QPushButton(VGenesStartUpDialog)
        self.cboOpen.setMinimumSize(QtCore.QSize(180, 45))
        self.cboOpen.setMaximumSize(QtCore.QSize(180, 16777215))
        icon1 = QtGui.QIcon()
        icon1.addPixmap(QtGui.QPixmap(":/PNG-Icons/folder.png"), QtGui.QIcon.Normal, QtGui.QIcon.Off)
        self.cboOpen.setIcon(icon1)
        self.cboOpen.setIconSize(QtCore.QSize(24, 24))
        self.cboOpen.setObjectName("cboOpen")
        self.gridLayout.addWidget(self.cboOpen, 3, 2, 1, 1)
        spacerItem3 = QtWidgets.QSpacerItem(95, 17, QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Minimum)
        self.gridLayout.addItem(spacerItem3, 3, 3, 1, 1)

        self.retranslateUi(VGenesStartUpDialog)
        QtCore.QMetaObject.connectSlotsByName(VGenesStartUpDialog)

    def retranslateUi(self, VGenesStartUpDialog):
        _translate = QtCore.QCoreApplication.translate
        VGenesStartUpDialog.setWindowTitle(_translate("VGenesStartUpDialog", "VGenes"))
        self.label_2.setText(_translate("VGenesStartUpDialog", "Choose a VGenes database file:"))
        self.label.setText(_translate("VGenesStartUpDialog", "Open recent:"))
        self.cboRecent.setToolTip(_translate("VGenesStartUpDialog", "Click file name to open"))
        self.cboRecent.setItemText(0, _translate("VGenesStartUpDialog", "Select a recently opened VGenes database file"))
        self.btnNew.setToolTip(_translate("VGenesStartUpDialog", "Start a new project"))
        self.btnNew.setText(_translate("VGenesStartUpDialog", "New"))
        self.cboOpen.setToolTip(_translate("VGenesStartUpDialog", "Find a VGenes database to open"))
        self.cboOpen.setText(_translate("VGenesStartUpDialog", "Open other"))

import VgenesResources_rc
