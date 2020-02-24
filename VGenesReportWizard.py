#!/usr/bin/env python


from PyQt5.QtWidgets import (QApplication, QGridLayout, QLabel, QLineEdit,
        QVBoxLayout, QWizard, QWizardPage)
from PyQt5 import QtCore
from PyQt5 import QtWidgets


def createIntroPage():
    page = QWizardPage()
    page.setTitle("Introduction")

    label = QLabel(
            "This wizard will help you create a custom report of your data. "
            )
    label.setWordWrap(True)

    layout = QVBoxLayout()
    layout.addWidget(label)
    page.setLayout(layout)

    return page


def createRegistrationPage():
    page = QWizardPage()
    page.setTitle("Report")
    page.setSubTitle("Please fill both fields.")

    nameLabel = QLabel("Name:")
    nameLineEdit = QLineEdit()

    emailLabel = QLabel("Email address:")
    emailLineEdit = QLineEdit()

    layout = QGridLayout()
    layout.addWidget(nameLabel, 0, 0)
    layout.addWidget(nameLineEdit, 0, 1)
    layout.addWidget(emailLabel, 1, 0)
    layout.addWidget(emailLineEdit, 1, 1)
    page.setLayout(layout)

    return page



def createReportDesignPage():
    page = QWizardPage()
    page.setTitle("Report Design")
    page.setSubTitle("Choose Fields to report.")

    listWidget = QtWidgets.QListWidget()
    listWidget.setGeometry(QtCore.QRect(130, 50, 131, 192))
    listWidget.setTabKeyNavigation(True)
    listWidget.setObjectName("FieldList")


    # nameLabel = QLabel("Name:")
    # nameLineEdit = QLineEdit()
    #
    # emailLabel = QLabel("Email address:")
    # emailLineEdit = QLineEdit()

    layout = QGridLayout()
    layout.addWidget(listWidget, 0, 0)
    # layout.addWidget(nameLineEdit, 0, 1)
    # layout.addWidget(emailLabel, 1, 0)
    # layout.addWidget(emailLineEdit, 1, 1)
    page.setLayout(layout)

    return page

def createConclusionPage():
    page = QWizardPage()
    page.setTitle("Conclusion")

    label = QLabel("You are now successfully registered. Have a nice day!")
    label.setWordWrap(True)

    layout = QVBoxLayout()
    layout.addWidget(label)
    page.setLayout(layout)

    return page


if __name__ == '__main__':

    import sys

    app = QApplication(sys.argv)

    wizard = QWizard()
    wizard.addPage(createIntroPage())
    wizard.addPage(createRegistrationPage())

    wizard.addPage(createReportDesignPage())
    wizard.addPage(createConclusionPage())

    wizard.setWindowTitle("VGenes Report Wizard")
    wizard.show()

    sys.exit(app.exec_())