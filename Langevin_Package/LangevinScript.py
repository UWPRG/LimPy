"""This is a script to create a GUI for inputs for Langevin Integrator."""
# import pdb
# import os
import sys
# import math
import csv

# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# import scipy as sp

# from mpl_toolkits.mplot3d import Axes3D
from PyQt4.QtCore import *
from PyQt4.QtGui import *

from LangevinFunctions import LIMD


class Form(QWidget):

    """Generate the GUI for receiving inputs for Langevin Integrator."""

    """You can skip to below because most of this is just GUI stuff.
    """

    def __init__(w):
        """Create the basic GUI with menu options to be selected."""
        QWidget.__init__(w)
        layout = QGridLayout()
        Label1 = QLabel('Steps')
        layout.addWidget(Label1, 1, 0)
        textbox1 = QLineEdit()
        layout.addWidget(textbox1, 1, 1)
        Label2 = QLabel('Time Step')
        layout.addWidget(Label2, 2, 0)
        textbox2 = QLineEdit()
        layout.addWidget(textbox2, 2, 1)
        Label4 = QLabel('T(K)')
        layout.addWidget(Label4, 3, 0)
        textbox4 = QLineEdit()
        layout.addWidget(textbox4, 3, 1)
        Label5 = QLabel('Mass')
        layout.addWidget(Label5, 4, 0)
        textbox5 = QLineEdit()
        layout.addWidget(textbox5, 4, 1)
        btn = QPushButton('Lock in Values')
        layout.addWidget(btn, 16, 0)
        cb = QCheckBox('Make Movie')
        layout.addWidget(cb, 14, 0)
        cb2 = QCheckBox('Read File for Inputs')
        layout.addWidget(cb2, 14, 2)
        Labelfn = QLabel('Data Filename:')
        layout.addWidget(Labelfn, 15, 0)
        textboxfn = QLineEdit()
        layout.addWidget(textboxfn, 15, 1)
        Labelrf = QLabel('Read Filename:')
        layout.addWidget(Labelrf, 15, 2)
        textboxrf = QLineEdit()
        layout.addWidget(textboxrf, 15, 3)
        combo = QComboBox()
        combo.addItems(["Choose a Potential Dimension", "1-D Potential",
                        "2-D Potential"])
        combo2 = QComboBox()
        combo2.addItems(["Choose a Method", "MD", "Metadynamics",
                        "Well-Tempered Metadynamics", "Infrequent WT MetaD"])
        global count
        global cc
        cc = 0
        count = 0

        def onActivated(text):
            global Label3
            global textbox3
            global Label31
            global textbox31
            global Label32
            global textbox32
            global Label33
            global textbox33
            global Labely
            global textboxy
            global Labely1
            global textboxy1
            global Labely2
            global textboxy2
            global Labely3
            global textboxy3
            global count
            # Dropdown menu for 1-D or 2-D Potential
            if (text == "1-D Potential"):
                if count == 1:
                    layout.removeWidget(Labely)
                    layout.removeWidget(textboxy)
                    Labely.deleteLater()
                    textboxy.deleteLater()
                    layout.removeWidget(Labely1)
                    layout.removeWidget(textboxy1)
                    Labely1.deleteLater()
                    textboxy1.deleteLater()
                    layout.removeWidget(Labely2)
                    layout.removeWidget(textboxy2)
                    Labely2.deleteLater()
                    textboxy2.deleteLater()
                    layout.removeWidget(Labely3)
                    layout.removeWidget(textboxy3)
                    Labely3.deleteLater()
                    textboxy3.deleteLater()
                    count = 0
                Label3 = QLabel('X0')
                layout.addWidget(Label3, 5, 0)
                textbox3 = QLineEdit()
                layout.addWidget(textbox3, 5, 1)
                Label31 = QLabel('Xmin')
                layout.addWidget(Label31, 6, 0)
                textbox31 = QLineEdit()
                layout.addWidget(textbox31, 6, 1)
                Label32 = QLabel('Xmax')
                layout.addWidget(Label32, 7, 0)
                textbox32 = QLineEdit()
                layout.addWidget(textbox32, 7, 1)
                Label33 = QLabel('X increment')
                layout.addWidget(Label33, 8, 0)
                textbox33 = QLineEdit()
                layout.addWidget(textbox33, 8, 1)
                w.setLayout(layout)
            if (text == "2-D Potential"):
                Label3 = QLabel('X0')
                layout.addWidget(Label3, 5, 0)
                textbox3 = QLineEdit()
                layout.addWidget(textbox3, 5, 1)
                Label31 = QLabel('Xmin')
                layout.addWidget(Label31, 6, 0)
                textbox31 = QLineEdit()
                layout.addWidget(textbox31, 6, 1)
                Label32 = QLabel('Xmax')
                layout.addWidget(Label32, 7, 0)
                textbox32 = QLineEdit()
                layout.addWidget(textbox32, 7, 1)
                Label33 = QLabel('X increment')
                layout.addWidget(Label33, 8, 0)
                textbox33 = QLineEdit()
                layout.addWidget(textbox33, 8, 1)
                Labely = QLabel('Y0')
                layout.addWidget(Labely, 9, 0)
                textboxy = QLineEdit()
                layout.addWidget(textboxy, 9, 1)
                Labely1 = QLabel('Ymin')
                layout.addWidget(Labely1, 10, 0)
                textboxy1 = QLineEdit()
                layout.addWidget(textboxy1, 10, 1)
                Labely2 = QLabel('Ymax')
                layout.addWidget(Labely2, 11, 0)
                textboxy2 = QLineEdit()
                layout.addWidget(textboxy2, 11, 1)
                Labely3 = QLabel('Y increment')
                layout.addWidget(Labely3, 12, 0)
                textboxy3 = QLineEdit()
                layout.addWidget(textboxy3, 12, 1)
                count = 1

        def onActive(text):
            """Create the dropdown menu options for the method."""
            global Label6
            global textbox6
            global Label7
            global textbox7
            global Label8
            global textbox8
            global Label9
            global textbox9
            global Label10
            global textbox10
            global cc
            # Dropdown menu for method selection
            if (text == "MD"):
                if cc == 1:
                    layout.removeWidget(Label6)
                    layout.removeWidget(textbox6)
                    Label6.deleteLater()
                    textbox6.deleteLater()
                    layout.removeWidget(Label7)
                    layout.removeWidget(textbox7)
                    Label7.deleteLater()
                    textbox7.deleteLater()
                    layout.removeWidget(Label8)
                    layout.removeWidget(textbox8)
                    Label8.deleteLater()
                    textbox8.deleteLater()
                if cc == 2:
                    layout.removeWidget(Label6)
                    layout.removeWidget(textbox6)
                    Label6.deleteLater()
                    textbox6.deleteLater()
                    layout.removeWidget(Label7)
                    layout.removeWidget(textbox7)
                    Label7.deleteLater()
                    textbox7.deleteLater()
                    layout.removeWidget(Label8)
                    layout.removeWidget(textbox8)
                    Label8.deleteLater()
                    textbox8.deleteLater()
                    layout.removeWidget(Label9)
                    layout.removeWidget(textbox9)
                    Label9.deleteLater()
                    textbox9.deleteLater()
                if cc == 3:
                    layout.removeWidget(Label6)
                    layout.removeWidget(textbox6)
                    Label6.deleteLater()
                    textbox6.deleteLater()
                    layout.removeWidget(Label7)
                    layout.removeWidget(textbox7)
                    Label7.deleteLater()
                    textbox7.deleteLater()
                    layout.removeWidget(Label8)
                    layout.removeWidget(textbox8)
                    Label8.deleteLater()
                    textbox8.deleteLater()
                    layout.removeWidget(Label9)
                    layout.removeWidget(textbox9)
                    Label9.deleteLater()
                    textbox9.deleteLater()
                    layout.removeWidget(Label10)
                    layout.removeWidget(textbox10)
                    Label10.deleteLater()
                    textbox10.deleteLater()
                cc = 0
            if (text == "Metadynamics"):
                if cc == 2:
                    layout.removeWidget(Label9)
                    layout.removeWidget(textbox9)
                    Label9.deleteLater()
                    textbox9.deleteLater()
                if cc == 3:
                    layout.removeWidget(Label9)
                    layout.removeWidget(textbox9)
                    Label9.deleteLater()
                    textbox9.deleteLater()
                    layout.removeWidget(Label10)
                    layout.removeWidget(textbox10)
                    Label10.deleteLater()
                    textbox10.deleteLater()
                Label6 = QLabel('Gaussian Height')
                layout.addWidget(Label6, 1, 3)
                textbox6 = QLineEdit()
                layout.addWidget(textbox6, 1, 2)
                Label7 = QLabel('Gaussian Width')
                layout.addWidget(Label7, 2, 3)
                textbox7 = QLineEdit()
                layout.addWidget(textbox7, 2, 2)
                Label8 = QLabel('Deposition Frequency')
                layout.addWidget(Label8, 3, 3)
                textbox8 = QLineEdit()
                layout.addWidget(textbox8, 3, 2)
                cc = 1
            if (text == "Well-Tempered Metadynamics"):
                if cc == 3:
                    layout.removeWidget(Label10)
                    layout.removeWidget(textbox10)
                    Label10.deleteLater()
                    textbox10.deleteLater()
                Label6 = QLabel('Gaussian Height')
                layout.addWidget(Label6, 1, 3)
                textbox6 = QLineEdit()
                layout.addWidget(textbox6, 1, 2)
                Label7 = QLabel('Gaussian Width')
                layout.addWidget(Label7, 2, 3)
                textbox7 = QLineEdit()
                layout.addWidget(textbox7, 2, 2)
                Label8 = QLabel('Deposition Frequency')
                layout.addWidget(Label8, 3, 3)
                textbox8 = QLineEdit()
                layout.addWidget(textbox8, 3, 2)
                Label9 = QLabel('Well Temperature')
                layout.addWidget(Label9, 4, 3)
                textbox9 = QLineEdit()
                layout.addWidget(textbox9, 4, 2)
                cc = 2
            if (text == "Infrequent WT MetaD"):
                Label6 = QLabel('Gaussian Height')
                layout.addWidget(Label6, 1, 3)
                textbox6 = QLineEdit()
                layout.addWidget(textbox6, 1, 2)
                Label7 = QLabel('Gaussian Width')
                layout.addWidget(Label7, 2, 3)
                textbox7 = QLineEdit()
                layout.addWidget(textbox7, 2, 2)
                Label8 = QLabel('Deposition Frequency')
                layout.addWidget(Label8, 3, 3)
                textbox8 = QLineEdit()
                layout.addWidget(textbox8, 3, 2)
                Label9 = QLabel('Well Temperature')
                layout.addWidget(Label9, 4, 3)
                textbox9 = QLineEdit()
                layout.addWidget(textbox9, 4, 2)
                Label10 = QLabel('Number of Events')
                layout.addWidget(Label10, 5, 3)
                textbox10 = QLineEdit()
                layout.addWidget(textbox10, 5, 2)
                cc = 3
        combo.activated[str].connect(onActivated)
        layout.addWidget(combo, 0, 1)
        combo2.activated[str].connect(onActive)
        layout.addWidget(combo2, 0, 2)
        w.setLayout(layout)
        w.connect(btn, SIGNAL("clicked()"), w, SLOT("close()"))

        def on_click():
            """Assign variables once button is clicked to lock in inputs."""
            global inps
            global movieflag
            global potdim
            global sm
            global mdps
            movieflag = 0
            # Read in inputs from csv file
            if cb2.isChecked():
                inputsfile = str(textboxrf.text())
                inputs = pd.read_csv(inputsfile+'.csv')
                for jp in range(inputs.shape[0]):
                    print ('Running System ' + str(jp+1) + ' of ' +
                           str(inputs.shape[0]))
                    potdim = str(inputs['Dimension'][jp])
                    sm = str(inputs['Method'][jp])
                    filetitle = str(inputs['Data Filename'][jp])
                    inps = np.zeros(12)
                    inps[0] = float(inputs['Steps'][jp])
                    inps[1] = float(inputs['Step size'][jp])
                    inps[2] = float(inputs['X0'][jp])
                    inps[3] = float(inputs['Temperature'][jp])
                    inps[4] = float(inputs['Mass'][jp])
                    inps[5] = float(inputs['Xmin'][jp])
                    inps[6] = float(inputs['Xmax'][jp])
                    inps[7] = float(inputs['Xincrement'][jp])
                    inps[8] = float(inputs['Y0'][jp])
                    inps[9] = float(inputs['Ymin'][jp])
                    inps[10] = float(inputs['Ymax'][jp])
                    inps[11] = float(inputs['Yincrement'][jp])

                    mdps = np.zeros(4)
                    mdps[0] = float(inputs['Gaussian Height'][jp])
                    mdps[1] = float(inputs['Gaussian Width'][jp])
                    mdps[2] = float(inputs['Deposition Frequency'][jp])
                    mdps[3] = float(inputs['Well Temperature'][jp])
                    # Make Movie if checked
                    if cb.isChecked():
                        movieflag = 1
                    # Infrequent Metadynamics Loop
                    if (sm == "Infrequent WT MetaD"):
                        loops = int(inputs['Number of Events'][jp])
                        INFdata = np.zeros([loops, 2])
                        for its in range(0, loops):
                            result = LIMD(inps, mdps, potdim, sm, movieflag)
                            INFdata[its, 0] = result[0]
                            INFdata[its, 1] = result[1]
                            print 'Trial ' + str(its+1) + 'of ' + str(loops)
                        np.savetxt(filetitle + '.csv', INFdata, delimiter=',')
                        with open(filetitle + 'info.csv', "wb") as f:
                                writer = csv.writer(f)
                                writer.writerow([result[2]])
                    # Get converged FES
                    else:
                            result = LIMD(inps, mdps, potdim, sm, movieflag)
                            with open(filetitle+'info.csv', "wb") as f:
                                writer = csv.writer(f)
                                writer.writerow(['RMSD', 'RMSDKLD',
                                                'RMSDAligned'])
                                writer.writerow([result[0], result[1],
                                                result[2]])
                                writer.writerow([result[3]])
                    print "Jobs Completed"
            # Read in inputs from the GUI
            else:
                sm = str(combo2.currentText())
                potdim = str(combo.currentText())
                filetitle = str(textboxfn.text())
                if (potdim == '1-D Potential'):
                    inps = np.array([float(textbox1.text()),
                                    float(textbox2.text()),
                                    float(textbox3.text()),
                                    float(textbox4.text()),
                                    float(textbox5.text()),
                                    float(textbox31.text()),
                                    float(textbox32.text()),
                                    float(textbox33.text())])
                elif (potdim == '2-D Potential'):
                    inps = np.array([float(textbox1.text()),
                                    float(textbox2.text()),
                                    float(textbox3.text()),
                                    float(textbox4.text()),
                                    float(textbox5.text()),
                                    float(textbox31.text()),
                                    float(textbox32.text()),
                                    float(textbox33.text()),
                                    float(textboxy.text()),
                                    float(textboxy1.text()),
                                    float(textboxy2.text()),
                                    float(textboxy3.text())])
                if (sm == 'Metadynamics'):
                    mdps = np.array([float(textbox6.text()),
                                    float(textbox7.text()),
                                    float(textbox8.text())])
                elif (sm == "MD"):
                        mdps = np.array([0])
                elif (sm == 'Well-Tempered Metadynamics'):
                    mdps = np.array([float(textbox6.text()),
                                    float(textbox7.text()),
                                    float(textbox8.text()),
                                    float(textbox9.text())])
                elif (sm == "Infrequent WT MetaD"):
                    mdps = np.array([float(textbox6.text()),
                                    float(textbox7.text()),
                                    float(textbox8.text()),
                                    float(textbox9.text())])
                    loops = int(textbox10.text())
                if cb.isChecked():
                    movieflag = 1
                # Infrequent Metadynamics Loop
                if (sm == "Infrequent WT MetaD"):
                    INFdata = np.zeros([loops, 2])
                    for its in range(0, loops):
                        result = LIMD(inps, mdps, potdim, sm, movieflag)
                        INFdata[its, 0] = result[0]
                        INFdata[its, 1] = result[1]
                        print (str(its+1) + ' Trial(s) completed of ' +
                               str(loops))
                    np.savetxt(filetitle+'.csv', INFdata, delimiter=',')
                    with open('info'+filetitle+'.csv', "wb") as f:
                        writer = csv.writer(f)
                        writer.writerow([result[2]])
                # Get converged FES
                else:
                        result = LIMD(inps, mdps, potdim, sm, movieflag)
                        with open('info'+filetitle+'.csv', "wb") as f:
                            writer = csv.writer(f)
                            writer.writerow(['RMSD', 'RMSDKLD', 'RMSDAligned'])
                            writer.writerow([result[0], result[1], result[2]])
                            writer.writerow([result[3]])

        btn.clicked.connect(on_click)

app = QApplication(sys.argv)
form = Form()
form.show()
app.exec_()
