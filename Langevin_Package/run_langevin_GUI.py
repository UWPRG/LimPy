"""This is a script to create a GUI for inputs for Langevin Integrator."""
# import pdb
import os
import sys
# import math
import csv

# import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import markdown
# import scipy as sp

# from mpl_toolkits.mplot3d import Axes3D
from PyQt4.QtCore import *
from PyQt4.QtGui import *

import langevin_functions as lf
import potential_functions as pf
from simulate1D import simulate_1Dsystem
from simulate2D import simulate_2Dsystem


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
        textbox1.setText('1000000')
        layout.addWidget(textbox1, 1, 1)
        Label2 = QLabel('Time Step')
        layout.addWidget(Label2, 2, 0)
        textbox2 = QLineEdit()
        textbox2.setText('0.01')
        layout.addWidget(textbox2, 2, 1)
        Label4 = QLabel('T(K)')
        layout.addWidget(Label4, 3, 0)
        textbox4 = QLineEdit()
        textbox4.setText('300')
        layout.addWidget(textbox4, 3, 1)
        Label5 = QLabel('Mass')
        layout.addWidget(Label5, 4, 0)
        textbox5 = QLineEdit()
        textbox5.setText('1')
        layout.addWidget(textbox5, 4, 1)
        btn = QPushButton('Lock in Values')
        layout.addWidget(btn, 16, 0)
        helpbtn = QPushButton('HELP')
        layout.addWidget(helpbtn, 16, 4)
        cbsp = QCheckBox('Show Plot')
        layout.addWidget(cbsp, 1, 4)
        Labelguf = QLabel('Graph Update Frequency (steps):')
        layout.addWidget(Labelguf, 3, 4)
        textboxguf = QLineEdit()
        textboxguf.setText('1000')
        layout.addWidget(textboxguf, 4, 4)
        cb = QCheckBox('Make Movie')
        layout.addWidget(cb, 2, 4)
        cb2 = QCheckBox('Read File for Inputs')
        layout.addWidget(cb2, 14, 3)
        Labelfn = QLabel('Data Filename:')
        layout.addWidget(Labelfn, 15, 0)
        textboxfn = QLineEdit()
        layout.addWidget(textboxfn, 15, 1)
        Labelrf = QLabel('Read Filename:')
        layout.addWidget(Labelrf, 15, 3)
        textboxrf = QLineEdit()
        layout.addWidget(textboxrf, 15, 4)
        Labelgm = QLabel('Gamma:')
        layout.addWidget(Labelgm, 13, 0)
        textboxgm = QLineEdit()
        textboxgm.setText('5')
        layout.addWidget(textboxgm, 13, 1)
        potential_options = pf.get_potential_dict()
        po = ["Choose a Potential Function"] + potential_options.keys()
        combo = QComboBox()
        combo.addItems(po)
        combo2 = QComboBox()
        combo2.addItems(["Choose a Method", "MD", "Metadynamics",
                        "Well-Tempered Metadynamics", "Infrequent WT MetaD"])
        combo3 = QComboBox()
        combo3.addItems(["Energy Units", "Kcal/mol", "KJ/mol",
                        "Joule/mol", "cal/mol"])
        global count
        global cc
        cc = 0
        count = 0

        def showhelp():
            msg = QMessageBox()
            msg.setText('Steps: Number of time steps in simulation' + 2*'\n' +
                        'Time step: Step size in units of time' + 2*'\n' +
                        'T(K): Temperature in units of Kelvin' + 2*'\n' +
                        'Mass: Mass' + 2*'\n' +
                        'Gamma: Friction Factor' + 2*'\n' +
                        'X0: Starting X value' + 2*'\n' +
                        'Xmin: Minimum Value on X axis' + 2*'\n' +
                        'Xmax: Maximum Value on X Axis' + 2*'\n' +
                        'Xinc: X axis grid spacing' + 2*'\n' +
                        'Y0: Starting Y value' + 2*'\n' +
                        'Ymin: Minimum Y value on X axis' + 2*'\n' +
                        'Ymax: Maximum Y value on X Axis' + 2*'\n' +
                        'Yinc: Y axis grid spacing' + '\n'
                        )
            retval = msg.exec_()

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
            if (text == "cosine_potential"):
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

                Label3 = QLabel('X0')
                layout.addWidget(Label3, 5, 0)
                textbox3 = QLineEdit()
                textbox3.setText('3.14')
                layout.addWidget(textbox3, 5, 1)
                Label31 = QLabel('Xmin')
                layout.addWidget(Label31, 6, 0)
                textbox31 = QLineEdit()
                textbox31.setText('-0.31')
                layout.addWidget(textbox31, 6, 1)
                Label32 = QLabel('Xmax')
                layout.addWidget(Label32, 7, 0)
                textbox32 = QLineEdit()
                textbox32.setText('6.5')
                layout.addWidget(textbox32, 7, 1)
                Label33 = QLabel('X increment')
                layout.addWidget(Label33, 8, 0)
                textbox33 = QLineEdit()
                textbox33.setText('0.01')
                layout.addWidget(textbox33, 8, 1)
                count = 0
            if text == "two_gaussian_potential":
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

                Label3 = QLabel('X0')
                layout.addWidget(Label3, 5, 0)
                textbox3 = QLineEdit()
                textbox3.setText('2.67')
                layout.addWidget(textbox3, 5, 1)
                Label31 = QLabel('Xmin')
                layout.addWidget(Label31, 6, 0)
                textbox31 = QLineEdit()
                textbox31.setText('-4.0')
                layout.addWidget(textbox31, 6, 1)
                Label32 = QLabel('Xmax')
                layout.addWidget(Label32, 7, 0)
                textbox32 = QLineEdit()
                textbox32.setText('4.0')
                layout.addWidget(textbox32, 7, 1)
                Label33 = QLabel('X increment')
                layout.addWidget(Label33, 8, 0)
                textbox33 = QLineEdit()
                textbox33.setText('0.01')
                layout.addWidget(textbox33, 8, 1)
                count = 0

            if (text == "pv_2D_potential"):

                Label3 = QLabel('X0')
                layout.addWidget(Label3, 5, 0)
                textbox3 = QLineEdit()
                textbox3.setText('1.5')
                layout.addWidget(textbox3, 5, 1)
                Label31 = QLabel('Xmin')
                layout.addWidget(Label31, 6, 0)
                textbox31 = QLineEdit()
                textbox31.setText('0')
                layout.addWidget(textbox31, 6, 1)
                Label32 = QLabel('Xmax')
                layout.addWidget(Label32, 7, 0)
                textbox32 = QLineEdit()
                textbox32.setText('3.0')
                layout.addWidget(textbox32, 7, 1)
                Label33 = QLabel('X increment')
                layout.addWidget(Label33, 8, 0)
                textbox33 = QLineEdit()
                textbox33.setText('0.01')
                layout.addWidget(textbox33, 8, 1)
                Labely = QLabel('Y0')
                layout.addWidget(Labely, 9, 0)
                textboxy = QLineEdit()
                textboxy.setText('0.6')
                layout.addWidget(textboxy, 9, 1)
                Labely1 = QLabel('Ymin')
                layout.addWidget(Labely1, 10, 0)
                textboxy1 = QLineEdit()
                textboxy1.setText('-2.0')
                layout.addWidget(textboxy1, 10, 1)
                Labely2 = QLabel('Ymax')
                layout.addWidget(Labely2, 11, 0)
                textboxy2 = QLineEdit()
                textboxy2.setText('2.0')
                layout.addWidget(textboxy2, 11, 1)
                Labely3 = QLabel('Y increment')
                layout.addWidget(Labely3, 12, 0)
                textboxy3 = QLineEdit()
                textboxy3.setText('0.01')
                layout.addWidget(textboxy3, 12, 1)
                count = 1

            if text == "mueller_brown_potential":
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
        layout.addWidget(combo3, 0, 4)
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
                inputsfile = (os.getcwd() + '/inputfiles/' +
                              str(textboxrf.text()))
                inputdata = lf.get_parameters(inputsfile)
                inps = inputdata[0]
                mdps = inputdata[1]
                dimension = inputdata[2]
                method = inputdata[3]
                potfunc = inputdata[4]
                filetitle = inputdata[5]
                makeplot = inputdata[6]
                plot_freq = int(inputdata[7])
                make_movie = inputdata[8]
                trials = int(mdps[-1])
                if method == "Infrequent WT MetaD":
                    for its in range(0, trials):
                        if dimension == '1-D Potential':
                            trial = simulate_1Dsystem(inps, mdps,
                                                      dimension, method,
                                                      potfunc, filetitle,
                                                      makeplot, plot_freq,
                                                      make_movie)
                        else:
                            trial = simulate_2Dsystem(inps, mdps,
                                                      dimension, method,
                                                      potfunc, filetitle,
                                                      makeplot, plot_freq,
                                                      make_movie)
                        if its == 0:
                            timedata = np.array([trial[0], trial[1]])
                        else:
                            timedata = np.append(timedata,
                                                 np.array([trial[0],
                                                          trial[1]]))
                    collect = np.asarray(timedata)
                    collect = np.reshape(collect, (num_iter*size, 2))
                    np.savetxt(filetitle+'_Allevents.csv', collect,
                               delimiter=',')
                    with open(filetitle + 'info.csv', "wb") as f:
                            writer = csv.writer(f)
                            writer.writerow([trial[2]])
                    # Get converged FES
                else:

                    if dimension == '1-D Potential':
                        trial = simulate_1Dsystem(inps, mdps,
                                                  dimension, method,
                                                  potfunc, filetitle,
                                                  makeplot, plot_freq,
                                                  make_movie)
                        colvar = pd.DataFrame({'CV': trial[0], 'E': trial[1]})
                        colvar.reset_index('CV')
                    else:
                        trial = simulate_2Dsystem(inps, mdps,
                                                  dimension, method,
                                                  potfunc, filetitle,
                                                  makeplot, plot_freq,
                                                  make_movie)
                        colvar = pd.DataFrame({'CV1': trial[0][:, 0],
                                               'CV2': trial[0][:, 1],
                                               'E': trial[1]})
                        colvar.reset_index('CV1')
                    colvar.index.name = 'Step'
                    if method == "MD":
                        colvar.to_csv(filetitle+'_MD.csv')
                    else:
                        colvar.to_csv(filetitle+'_COLVAR.csv')
                    with open(filetitle+'info.csv', "wb") as f:
                        writer = csv.writer(f)
                        writer.writerow(['RMSD', 'RMSDKLD',
                                        'RMSDAligned'])
                        writer.writerow([trial[2]])
                        writer.writerow([trial[3]])
                    with open(filetitle+'info.csv', "wb") as f:
                        writer = csv.writer(f)
                        writer.writerow([trial[0]])
                        writer.writerow([trial[1]])

            else:
                if cbsp.isChecked():
                    makeplot = 'True'
                else:
                    makeplot = 'False'
                method = str(combo2.currentText())
                potential = str(combo.currentText())
                units = str(combo3.currentText())
                if units == "Kcal/mol" or units == "Energy Units":
                    kb = 0.00198
                elif units == "KJ/mol":
                    kb = 0.008314
                elif units == "Joule/mol":
                    kb = 8.314
                elif units == "cal/mol":
                    kb = 1.98
                if (potential == "cosine_potential" or
                   potential == "two_gaussian_potential"):
                    dimension = '1-D Potential'

                else:
                    dimension = '2-D Potential'
                filetitle = str(textboxfn.text())
                if (dimension == '1-D Potential'):
                    inps = np.array([float(textbox1.text()),
                                    float(textbox2.text()),
                                    float(textbox3.text()),
                                    float(textbox4.text()),
                                    float(textbox5.text()),
                                    float(textbox31.text()),
                                    float(textbox32.text()),
                                    float(textbox33.text()), 0, 0, 0, 0,
                                    float(textboxgm.text()),
                                    float(kb)])
                elif (dimension == '2-D Potential'):
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
                                    float(textboxy3.text()),
                                    float(textboxgm.text()),
                                    float(kb)])
                if (method == 'Metadynamics'):
                    mdps = np.array([float(textbox6.text()),
                                    float(textbox7.text()),
                                    float(textbox8.text())])
                elif (method == "MD"):
                        mdps = np.array([0])
                elif (method == 'Well-Tempered Metadynamics'):
                    mdps = np.array([float(textbox6.text()),
                                    float(textbox7.text()),
                                    float(textbox8.text()),
                                    float(textbox9.text())])
                elif (method == "Infrequent WT MetaD"):
                    mdps = np.array([float(textbox6.text()),
                                    float(textbox7.text()),
                                    float(textbox8.text()),
                                    float(textbox9.text())])
                    loops = int(textbox10.text())
                plot_freq = float(textboxguf.text())
                if cb.isChecked():
                    make_movie = 'True'
                    makeplot = 'True'
                else:
                    make_movie = 'False'
                # Infrequent Metadynamics Loop
                if (method == "Infrequent WT MetaD"):
                    for its in range(0, loops):
                        if dimension == '1-D Potential':
                            trial = simulate_1Dsystem(inps, mdps,
                                                      dimension, method,
                                                      potential, filetitle,
                                                      makeplot, plot_freq,
                                                      make_movie)
                        else:
                            trial = simulate_2Dsystem(inps, mdps,
                                                      dimension, method,
                                                      potential, filetitle,
                                                      makeplot, plot_freq,
                                                      make_movie)
                        if its == 0:
                            timedata = np.array([trial[0], trial[1]])
                        else:
                            timedata = np.append(timedata,
                                                 np.array([trial[0],
                                                          trial[1]]))
                    collect = np.asarray(timedata)
                    collect = np.reshape(collect, (loops, 2))
                    np.savetxt(filetitle+'_Allevents.csv', collect,
                               delimiter=',')
                    with open(filetitle + 'info.csv', "wb") as f:
                        writer = csv.writer(f)
                        writer.writerow([trial[2]])
                    # Get converged FES
                else:
                    if dimension == '1-D Potential':
                        trial = simulate_1Dsystem(inps, mdps,
                                                  dimension, method,
                                                  potential, filetitle,
                                                  makeplot, plot_freq,
                                                  make_movie)
                        colvar = pd.DataFrame({'CV': trial[0], 'E': trial[1]})
                        colvar.reset_index('CV')
                    else:
                        trial = simulate_2Dsystem(inps, mdps,
                                                  dimension, method,
                                                  potential, filetitle,
                                                  makeplot, plot_freq,
                                                  make_movie)
                        colvar = pd.DataFrame({'CV1': trial[0][:, 0],
                                               'CV2': trial[0][:, 1],
                                               'E': trial[1]})
                        colvar.reset_index('CV1')
                    colvar.index.name = 'Step'
                    if method == "MD":
                        colvar.to_csv(filetitle+'_MD.csv')
                    else:
                        colvar.to_csv(filetitle+'_COLVAR.csv')
                    with open(filetitle+'info.csv', "wb") as f:
                        writer = csv.writer(f)
                        writer.writerow(['RMSD', 'RMSDKLD',
                                        'RMSDAligned'])
                        writer.writerow([trial[2]])
                        writer.writerow([trial[3]])
                    with open(filetitle+'info.csv', "wb") as f:
                        writer = csv.writer(f)
                        writer.writerow([trial[0]])
                        writer.writerow([trial[1]])

        btn.clicked.connect(on_click)
        helpbtn.clicked.connect(showhelp)
app = QApplication(sys.argv)
form = Form()
form.show()
app.exec_()
