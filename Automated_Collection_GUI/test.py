# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'variable_collection_mainwindow.ui'
#
# Created by: PyQt5 UI code generator 5.15.9
#
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets


class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(657, 494)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.title = QtWidgets.QLabel(self.centralwidget)
        self.title.setGeometry(QtCore.QRect(50, 0, 561, 111))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(36)
        self.title.setFont(font)
        self.title.setLayoutDirection(QtCore.Qt.LeftToRight)
        self.title.setIndent(-1)
        self.title.setObjectName("title")
        self.gene_name_box = QtWidgets.QLineEdit(self.centralwidget)
        self.gene_name_box.setGeometry(QtCore.QRect(50, 100, 271, 61))
        self.gene_name_box.setObjectName("gene_name_box")
        self.gene_label = QtWidgets.QLabel(self.centralwidget)
        self.gene_label.setGeometry(QtCore.QRect(100, 170, 161, 20))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(20)
        self.gene_label.setFont(font)
        self.gene_label.setObjectName("gene_label")
        self.mutations_label = QtWidgets.QLabel(self.centralwidget)
        self.mutations_label.setGeometry(QtCore.QRect(250, 290, 201, 31))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(20)
        self.mutations_label.setFont(font)
        self.mutations_label.setObjectName("mutations_label")
        self.submit_button = QtWidgets.QPushButton(self.centralwidget)
        self.submit_button.setGeometry(QtCore.QRect(290, 340, 100, 32))
        self.submit_button.setObjectName("submit_button")
        self.checkbox_scratch = QtWidgets.QCheckBox(self.centralwidget)
        self.checkbox_scratch.setGeometry(QtCore.QRect(60, 220, 101, 20))
        self.checkbox_scratch.setObjectName("checkbox_scratch")
        self.checkbox_update = QtWidgets.QCheckBox(self.centralwidget)
        self.checkbox_update.setGeometry(QtCore.QRect(60, 250, 71, 20))
        self.checkbox_update.setObjectName("checkbox_update")
        self.mutations_box = QtWidgets.QLineEdit(self.centralwidget)
        self.mutations_box.setGeometry(QtCore.QRect(200, 220, 271, 61))
        self.mutations_box.setObjectName("mutations_box")
        self.gene_button = QtWidgets.QPushButton(self.centralwidget)
        self.gene_button.setGeometry(QtCore.QRect(440, 110, 161, 41))
        self.gene_button.setObjectName("gene_button")
        self.gene_file_label = QtWidgets.QLabel(self.centralwidget)
        self.gene_file_label.setGeometry(QtCore.QRect(460, 160, 131, 20))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(18)
        self.gene_file_label.setFont(font)
        self.gene_file_label.setObjectName("gene_file_label")
        self.delimit_label = QtWidgets.QLabel(self.centralwidget)
        self.delimit_label.setGeometry(QtCore.QRect(250, 320, 201, 16))
        self.delimit_label.setObjectName("delimit_label")
        self.existing_button = QtWidgets.QPushButton(self.centralwidget)
        self.existing_button.setGeometry(QtCore.QRect(60, 270, 131, 32))
        font = QtGui.QFont()
        font.setFamily("Arial")
        font.setPointSize(13)
        self.existing_button.setFont(font)
        self.existing_button.setObjectName("existing_button")
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 657, 37))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "MainWindow"))
        self.title.setText(_translate("MainWindow", "Epitope Clinical Variable Collection"))
        self.gene_label.setText(_translate("MainWindow", "Input Gene Name"))
        self.mutations_label.setText(_translate("MainWindow", "Input Point Mutations"))
        self.submit_button.setText(_translate("MainWindow", "Submit"))
        self.checkbox_scratch.setText(_translate("MainWindow", "From Scratch"))
        self.checkbox_update.setText(_translate("MainWindow", "Update"))
        self.gene_button.setText(_translate("MainWindow", "Input Gene File"))
        self.gene_file_label.setText(_translate("MainWindow", "Input Gene File"))
        self.delimit_label.setText(_translate("MainWindow", "delimit w/commas and no spaces"))
        self.existing_button.setText(_translate("MainWindow", "Existing Data (.xlsx)"))


if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())
