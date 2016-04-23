"""This is a script to create a GUI for inputs for Langevin Integrator."""

import langevin_functions as lf
import potential_class
import boundarycondition
from simulate1D import simulate_1Dsystem
from simulate2D import simulate_2Dsystem
from statistical_functions import perform_ks_analysis, sampling
import sys
import os
import numpy as np
import pandas as pd
import pdb
import csv

from PyQt4.QtCore import *
from PyQt4.QtGui import *


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
        layout.addWidget(btn, 20, 0)
        helpbtn = QPushButton('HELP')
        layout.addWidget(helpbtn, 16, 4)
        cbsp = QCheckBox('Show Plot')
        layout.addWidget(cbsp, 1, 4)
        # cb2 = QCheckBox('Read File for Inputs')
        # layout.addWidget(cb2, 14, 3)
        Labelfn = QLabel('Data Filename:')
        layout.addWidget(Labelfn, 19, 0)
        textboxfn = QLineEdit()
        layout.addWidget(textboxfn, 19, 1)
        # Labelrf = QLabel('Read Filename:')
        # layout.addWidget(Labelrf, 15, 3)
        # textboxrf = QLineEdit()
        # layout.addWidget(textboxrf, 15, 4)
        Labelgm = QLabel('Gamma:')
        layout.addWidget(Labelgm, 5, 0)
        textboxgm = QLineEdit()
        textboxgm.setText('5')
        layout.addWidget(textboxgm, 5, 1)
        potential_options = potential_class.get_potential_dict()
        po = ["Potential Function"] + potential_options.keys()
        combo = QComboBox()
        combo.addItems(po)
        combo2 = QComboBox()
        combo2.addItems(["Method", "MD", "Metadynamics",
                        "Well-Tempered Metadynamics", "Infrequent WT MetaD"])
        combo3 = QComboBox()
        combo3.addItems(["Energy Units", "Kcal/mol", "KJ/mol",
                        "Joule/mol", "cal/mol"])
        combo4 = QComboBox()
        combo4.addItems(["Boundary Conditions", "Periodic", "Quartic" ,
                         "No BC"])
        global count
        global cc
        global countpp
        global bccount
        bccount = 0
        cc = 0
        count = 0
        countpp = 0

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
            global Labelpp1
            global Labelpp2
            global textboxpp1
            global textboxpp2
            global Labelpp3
            global Labelpp4
            global textboxpp3
            global textboxpp4
            global countpp
            global pot
            global paramlength
            # Dropdown menu for 1-D or 2-D Potential

            pot_dict= potential_class.get_potential_dict()
            pot = pot_dict[str(text)]
            pot = pot()
            paramlength = len(pot.parameters)
            pdim = pot.dimension
            psets = potential_class.get_GUI_presets_dict()
            try:
                presetvals = psets[str(text)]

            except:
                presetvals = np.array([0,0,0,0,0,0,0,0]).astype(str)
            if (pdim == "1-D Potential"):
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
                if paramlength == 4:
                    Labelpp1 = QLabel('Potential Parameter 1:')
                    layout.addWidget(Labelpp1, 6, 0)
                    textboxpp1 = QLineEdit()
                    layout.addWidget(textboxpp1, 6, 1)
                    Labelpp2 = QLabel('Potential Parameter 2:')
                    layout.addWidget(Labelpp2, 7, 0)
                    textboxpp2 = QLineEdit()
                    layout.addWidget(textboxpp2, 7, 1)
                    Labelpp3 = QLabel('Potential Parameter 3:')
                    layout.addWidget(Labelpp3, 8, 0)
                    textboxpp3 = QLineEdit()
                    layout.addWidget(textboxpp3, 8, 1)
                    Labelpp4 = QLabel('Potential Parameter 4:')
                    layout.addWidget(Labelpp4, 9, 0)
                    textboxpp4 = QLineEdit()
                    layout.addWidget(textboxpp4, 9, 1)
                    countpp == 1
                if paramlength == 2:
                    Labelpp1 = QLabel('Potential Parameter 1:')
                    layout.addWidget(Labelpp1, 6, 0)
                    textboxpp1 = QLineEdit()
                    layout.addWidget(textboxpp1, 6, 1)
                    Labelpp2 = QLabel('Potential Parameter 2:')
                    layout.addWidget(Labelpp2, 7, 0)
                    textboxpp2 = QLineEdit()
                    layout.addWidget(textboxpp2, 7, 1)
                    if countpp == 1:
                        layout.removeWidget(Labelpp3)
                        layout.removeWidget(textboxpp3)
                        Labelpp3.deleteLater()
                        textboxpp3.deleteLater()
                        layout.removeWidget(Labelpp4)
                        layout.removeWidget(textboxpp4)
                        Labelpp4.deleteLater()
                        textboxpp4.deleteLater()
                    countpp = 0
                Label3 = QLabel('X0')
                layout.addWidget(Label3, 10, 0)
                textbox3 = QLineEdit()
                textbox3.setText(presetvals[0])
                layout.addWidget(textbox3, 10, 1)
                Label31 = QLabel('Xmin')
                layout.addWidget(Label31, 11, 0)
                textbox31 = QLineEdit()
                textbox31.setText(presetvals[1])
                layout.addWidget(textbox31, 11, 1)
                Label32 = QLabel('Xmax')
                layout.addWidget(Label32, 12, 0)
                textbox32 = QLineEdit()
                textbox32.setText(presetvals[2])
                layout.addWidget(textbox32, 12, 1)
                Label33 = QLabel('X increment')
                layout.addWidget(Label33, 13, 0)
                textbox33 = QLineEdit()
                textbox33.setText(presetvals[3])
                layout.addWidget(textbox33, 13, 1)
                count = 0

            if (pdim == "2-D Potential"):

                if paramlength == 4:
                    Labelpp1 = QLabel('Potential Parameter 1:')
                    layout.addWidget(Labelpp1, 6, 0)
                    textboxpp1 = QLineEdit()
                    textboxpp1.setText('1')
                    layout.addWidget(textboxpp1, 6, 1)
                    Labelpp2 = QLabel('Potential Parameter 2:')
                    layout.addWidget(Labelpp2, 7, 0)
                    textboxpp2 = QLineEdit()
                    textboxpp2.setText('0')
                    layout.addWidget(textboxpp2, 7, 1)
                    Labelpp3 = QLabel('Potential Parameter 3:')
                    layout.addWidget(Labelpp3, 8, 0)
                    textboxpp3 = QLineEdit()
                    textboxpp3.setText('1')
                    layout.addWidget(textboxpp3, 8, 1)
                    Labelpp4 = QLabel('Potential Parameter 4:')
                    layout.addWidget(Labelpp4, 9, 0)
                    textboxpp4 = QLineEdit()
                    textboxpp4.setText('0')
                    layout.addWidget(textboxpp4, 9, 1)
                    countpp == 1
                if paramlength == 2:
                    Labelpp1 = QLabel('Potential Parameter 1:')
                    layout.addWidget(Labelpp1, 6, 0)
                    textboxpp1 = QLineEdit()
                    textboxpp1.setText('0')
                    layout.addWidget(textboxpp1, 6, 1)
                    Labelpp2 = QLabel('Potential Parameter 2:')
                    layout.addWidget(Labelpp2, 7, 0)
                    textboxpp2 = QLineEdit()
                    textboxpp2.setText('0')
                    layout.addWidget(textboxpp2, 7, 1)
                    if countpp == 1:
                        layout.removeWidget(Labelpp3)
                        layout.removeWidget(textboxpp3)
                        Labelpp3.deleteLater()
                        textboxpp3.deleteLater()
                        layout.removeWidget(Labelpp4)
                        layout.removeWidget(textboxpp4)
                        Labelpp4.deleteLater()
                        textboxpp4.deleteLater()
                    countpp = 0
                Label3 = QLabel('X0')
                layout.addWidget(Label3, 10, 0)
                textbox3 = QLineEdit()
                textbox3.setText(presetvals[0])
                layout.addWidget(textbox3, 10, 1)
                Label31 = QLabel('Xmin')
                layout.addWidget(Label31, 11, 0)
                textbox31 = QLineEdit()
                textbox31.setText(presetvals[1])
                layout.addWidget(textbox31, 11, 1)
                Label32 = QLabel('Xmax')
                layout.addWidget(Label32, 12, 0)
                textbox32 = QLineEdit()
                textbox32.setText(presetvals[2])
                layout.addWidget(textbox32, 12, 1)
                Label33 = QLabel('X increment')
                layout.addWidget(Label33, 13, 0)
                textbox33 = QLineEdit()
                textbox33.setText(presetvals[3])
                layout.addWidget(textbox33, 13, 1)
                Labely = QLabel('Y0')
                layout.addWidget(Labely, 14, 0)
                textboxy = QLineEdit()
                textboxy.setText(presetvals[4])
                layout.addWidget(textboxy, 14, 1)
                Labely1 = QLabel('Ymin')
                layout.addWidget(Labely1, 15, 0)
                textboxy1 = QLineEdit()
                textboxy1.setText(presetvals[5])
                layout.addWidget(textboxy1, 15, 1)
                Labely2 = QLabel('Ymax')
                layout.addWidget(Labely2, 16, 0)
                textboxy2 = QLineEdit()
                textboxy2.setText(presetvals[6])
                layout.addWidget(textboxy2, 16, 1)
                Labely3 = QLabel('Y increment')
                layout.addWidget(Labely3, 17, 0)
                textboxy3 = QLineEdit()
                textboxy3.setText(presetvals[7])
                layout.addWidget(textboxy3, 17, 1)
                count = 1

        def on_BCActive(text):
            global pot
            global bccount
            global cbxbclow
            global cbxbchigh
            global cbybclow
            global cbybchigh
            global bc_actx
            global bc_acty
            if text == 'Quartic':
                if pot.dimension == '1-D Potential':
                    if bccount == 2:
                        layout.removeWidget(cbybclow)
                        layout.removeWidget(cbybchigh)
                        cbybclow.deleteLater()
                        cbybchigh.deleteLater()

                    elif bccount == 4:
                        layout.removeWidget(bc_actx)
                        layout.removeWidget(bc_acty)
                        bc_actx.deleteLater()
                        bc_acty.deleteLater()
                    elif bccount == 3:
                        layout.removeWidget(bc_actx)
                        bc_actx.deleteLater()

                    cbxbclow = QCheckBox('Enable Lower X BC')
                    layout.addWidget(cbxbclow , 9, 3)
                    cbxbchigh = QCheckBox('Enable Higher X BC')
                    layout.addWidget(cbxbchigh, 10, 3)
                    bccount = 1
                else:
                    if bccount == 1:
                        layout.removeWidget(cbxbclow)
                        layout.removeWidget(cbxbchigh)
                        cbxbclow.deleteLater()
                        cbxbchigh.deleteLater()
                    elif bccount == 4:
                        layout.removeWidget(bc_actx)
                        layout.removeWidget(bc_acty)
                        bc_actx.deleteLater()
                        bc_acty.deleteLater()
                    elif bccount == 3:
                        layout.removeWidget(bc_actx)
                        bc_actx.deleteLater()
                    cbxbclow = QCheckBox('Enable Lower X BC')
                    layout.addWidget(cbxbclow , 9, 3)
                    cbxbchigh = QCheckBox('Enable Higher X BC')
                    layout.addWidget(cbxbchigh, 10, 3)
                    cbybclow = QCheckBox('Enable Lower Y BC')
                    layout.addWidget(cbybclow , 11, 3)
                    cbybchigh = QCheckBox('Enable High Y BC')
                    layout.addWidget(cbybchigh, 12, 3)
                    bccount = 2

            elif text == 'Periodic':
                if bccount == 1:
                    layout.removeWidget(cbxbclow)
                    layout.removeWidget(cbxbchigh)
                    cbxbclow.deleteLater()
                    cbxbchigh.deleteLater()
                elif bccount == 2:
                    layout.removeWidget(cbxbclow)
                    layout.removeWidget(cbxbchigh)
                    cbxbclow.deleteLater()
                    cbxbchigh.deleteLater()
                    layout.removeWidget(cbybclow)
                    layout.removeWidget(cbybchigh)
                    cbybclow.deleteLater()
                    cbybchigh.deleteLater()
                bc_actx = QCheckBox('Enable X BC')
                layout.addWidget(bc_actx, 9, 3)
                bccount = 3
                if pot.dimension == '2-D Potential':
                    bc_acty = QCheckBox('Enable Y BC')
                    layout.addWidget(bc_acty, 10, 3)
                    bccount = 4
            else:
                if bccount == 1:
                    layout.removeWidget(cbxbclow)
                    layout.removeWidget(cbxbchigh)
                    cbxbclow.deleteLater()
                    cbxbchigh.deleteLater()
                elif bccount == 2:
                    layout.removeWidget(cbxbclow)
                    layout.removeWidget(cbxbchigh)
                    cbxbclow.deleteLater()
                    cbxbchigh.deleteLater()
                    layout.removeWidget(cbybclow)
                    layout.removeWidget(cbybchigh)
                    cbybclow.deleteLater()
                    cbybchigh.deleteLater()
                elif bccount == 3:
                    layout.removeWidget(bc_actx)
                    bc_actx.deleteLater()
                elif bccount == 4:
                    layout.removeWidget(bc_actx)
                    layout.removeWidget(bc_acty)
                    bc_actx.deleteLater()
                    bc_acty.deleteLater()
                bccount = 0


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
            global LabelRE1
            global textboxRE1
            global LabelRE2
            global textboxRE2
            global pot
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
                    layout.removeWidget(LabelRE1)
                    layout.removeWidget(textboxRE1)
                    LabelRE1.deleteLater()
                    textboxRE1.deleteLater()
                    layout.removeWidget(LabelRE2)
                    layout.removeWidget(textboxRE2)
                    LabelRE2.deleteLater()
                    textboxRE2.deleteLater()
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
                    layout.removeWidget(LabelRE1)
                    layout.removeWidget(textboxRE1)
                    LabelRE1.deleteLater()
                    textboxRE1.deleteLater()
                    layout.removeWidget(LabelRE2)
                    layout.removeWidget(textboxRE2)
                    LabelRE2.deleteLater()
                    textboxRE2.deleteLater()
                Label6 = QLabel('Gaussian Height')
                layout.addWidget(Label6, 1, 2)
                textbox6 = QLineEdit()
                layout.addWidget(textbox6, 1, 3)
                Label7 = QLabel('Gaussian Width')
                layout.addWidget(Label7, 2, 2)
                textbox7 = QLineEdit()
                layout.addWidget(textbox7, 2, 3)
                Label8 = QLabel('Deposition Frequency')
                layout.addWidget(Label8, 3, 2)
                textbox8 = QLineEdit()
                layout.addWidget(textbox8, 3, 3)
                cc = 1
            if (text == "Well-Tempered Metadynamics"):
                if cc == 4:
                    layout.removeWidget(Label10)
                    layout.removeWidget(textbox10)
                    Label10.deleteLater()
                    textbox10.deleteLater()
                    layout.removeWidget(LabelRE1)
                    layout.removeWidget(textboxRE1)
                    LabelRE1.deleteLater()
                    textboxRE1.deleteLater()
                    layout.removeWidget(LabelRE2)
                    layout.removeWidget(textboxRE2)
                    LabelRE2.deleteLater()
                    textboxRE2.deleteLater()
                    layout.removeWidget(LabelRE3)
                    layout.removeWidget(textboxRE3)
                    LabelRE3.deleteLater()
                    textboxRE3.deleteLater()
                    layout.removeWidget(LabelRE4)
                    layout.removeWidget(textboxRE4)
                    LabelRE4.deleteLater()
                    textboxRE4.deleteLater()
                if cc == 3:
                    layout.removeWidget(Label10)
                    layout.removeWidget(textbox10)
                    Label10.deleteLater()
                    textbox10.deleteLater()
                    layout.removeWidget(LabelRE1)
                    layout.removeWidget(textboxRE1)
                    LabelRE1.deleteLater()
                    textboxRE1.deleteLater()
                    layout.removeWidget(LabelRE2)
                    layout.removeWidget(textboxRE2)
                    LabelRE2.deleteLater()
                    textboxRE2.deleteLater()
                Label6 = QLabel('Gaussian Height')
                layout.addWidget(Label6, 1, 2)
                textbox6 = QLineEdit()
                layout.addWidget(textbox6, 1, 3)
                Label7 = QLabel('Gaussian Width')
                layout.addWidget(Label7, 2, 2)
                textbox7 = QLineEdit()
                layout.addWidget(textbox7, 2, 3)
                Label8 = QLabel('Deposition Frequency')
                layout.addWidget(Label8, 3, 2)
                textbox8 = QLineEdit()
                layout.addWidget(textbox8, 3, 3)
                Label9 = QLabel('Bias Factor')
                layout.addWidget(Label9, 4, 2)
                textbox9 = QLineEdit()
                layout.addWidget(textbox9, 4, 3)
                cc = 2
            if (text == "Infrequent WT MetaD"):
                Label6 = QLabel('Gaussian Height')
                layout.addWidget(Label6, 1, 2)
                textbox6 = QLineEdit()
                layout.addWidget(textbox6, 1, 3)
                Label7 = QLabel('Gaussian Width')
                layout.addWidget(Label7, 2, 2)
                textbox7 = QLineEdit()
                layout.addWidget(textbox7, 2, 3)
                Label8 = QLabel('Deposition Frequency')
                layout.addWidget(Label8, 3, 2)
                textbox8 = QLineEdit()
                layout.addWidget(textbox8, 3, 3)
                Label9 = QLabel('Bias Factor')
                layout.addWidget(Label9, 4, 2)
                textbox9 = QLineEdit()
                layout.addWidget(textbox9, 4, 3)
                Label10 = QLabel('Number of Events')
                layout.addWidget(Label10, 5, 2)
                textbox10 = QLineEdit()
                layout.addWidget(textbox10, 5, 3)
                LabelRE1 = QLabel('X Lower Rare Event')
                layout.addWidget(LabelRE1, 6, 2)
                textboxRE1 = QLineEdit()
                layout.addWidget(textboxRE1, 6, 3)
                LabelRE2 = QLabel('X Upper Rare Event')
                layout.addWidget(LabelRE2, 7, 2)
                textboxRE2 = QLineEdit()
                layout.addWidget(textboxRE2, 7, 3)
                if cc == 4 and pot.dimension != '2-D Potential':
                    LabelRE3.deleteLater()
                    textboxRE3.deleteLater()
                    layout.removeWidget(LabelRE4)
                    layout.removeWidget(textboxRE4)
                    LabelRE4.deleteLater()
                    textboxRE4.deleteLater()
                if pot.dimension == '2-D Potential':
                    LabelRE3 = QLabel('Y Lower Rare Event')
                    layout.addWidget(LabelRE3, 8, 2)
                    textboxRE3 = QLineEdit()
                    layout.addWidget(textboxRE3, 8, 3)
                    LabelRE4 = QLabel('Y Upper Rare Event')
                    layout.addWidget(LabelRE4, 9, 2)
                    textboxRE4 = QLineEdit()
                    layout.addWidget(textboxRE4, 9, 3)
                    cc = 4
                cc = 3

        def on_checked():
            global plotit
            global cb
            global Labelguf
            global textboxguf
            global Labelmine
            global textboxmine
            global Labelmaxe
            global textboxmaxe
            if cbsp.isChecked():
                plotit = 0
                cb = QCheckBox('Make Movie')
                layout.addWidget(cb, 2, 4)
                Labelguf = QLabel('Graph Update Frequency (steps):')
                layout.addWidget(Labelguf, 3, 4)
                textboxguf = QLineEdit()
                textboxguf.setText('1000')
                layout.addWidget(textboxguf, 3, 5)
                Labelmine = QLabel('Min Energy Value to Plot:')
                layout.addWidget(Labelmine, 4, 4)
                textboxmine = QLineEdit()
                textboxmine.setText('0')
                layout.addWidget(textboxmine, 4, 5)
                Labelmaxe = QLabel('Min Energy Value to Plot:')
                layout.addWidget(Labelmaxe, 5, 4)
                textboxmaxe = QLineEdit()
                textboxmaxe.setText('100')
                layout.addWidget(textboxmaxe, 5, 5)
            elif plotit == 0:
                layout.removeWidget(cb)
                layout.removeWidget(Labelguf)
                layout.removeWidget(textboxguf)
                layout.removeWidget(Labelmine)
                layout.removeWidget(textboxmine)
                layout.removeWidget(Labelmaxe)
                layout.removeWidget(textboxmaxe)
                cb.deleteLater()
                Labelguf.deleteLater()
                textboxguf.deleteLater()
                Labelmine.deleteLater()
                textboxmine.deleteLater()
                Labelmaxe.deleteLater()
                textboxmaxe.deleteLater()
                plotit = 1

        combo.activated[str].connect(onActivated)
        layout.addWidget(combo, 0, 1)
        combo2.activated[str].connect(onActive)
        combo4.activated[str].connect(on_BCActive)
        cbsp.stateChanged.connect(on_checked)
        layout.addWidget(combo2, 0, 3)
        layout.addWidget(combo3, 0, 4)
        layout.addWidget(combo4,8,3)
        w.setLayout(layout)
        w.resize(1000,500)
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
            # if cb2.isChecked():
            #     inputsfile = (os.getcwd() + '/inputfiles/' +
            #                   str(textboxrf.text()))
            #     inputdata = lf.get_parameters(inputsfile)
            #     inps = inputdata[0]
            #     mdps = inputdata[1]
            #     dimension = inputdata[2]
            #     method = inputdata[3]
            #     potfunc = inputdata[4]
            #     filetitle = inputdata[5]
            #     makeplot = inputdata[6]
            #     plot_freq = int(inputdata[7])
            #     make_movie = inputdata[8]
            #     trials = int(mdps[-1])
            #     if method == "Infrequent WT MetaD":
            #         for its in range(0, trials):
            #             if dimension == '1-D Potential':
            #                 trial = simulate_1Dsystem(inps, mdps,
            #                                           dimension, method,
            #                                           potfunc, filetitle,
            #                                           makeplot, plot_freq,
            #                                           make_movie)
            #             else:
            #                 trial = simulate_2Dsystem(inps, mdps,
            #                                           dimension, method,
            #                                           potfunc, filetitle,
            #                                           makeplot, plot_freq,
            #                                           make_movie)
            #             if its == 0:
            #                 timedata = np.array([trial[0], trial[1],trial[3]])
            #             else:
            #                 timedata = np.append(timedata,
            #                                      np.array([trial[0],
            #                                               trial[1],trial[3]]))
            #         collect = np.asarray(timedata)
            #         collect = np.reshape(collect, (num_iter*size, 2))
            #         np.savetxt(filetitle+'_Allevents.csv', collect,
            #                    delimiter=',')
            #         with open(filetitle + 'info.csv', "wb") as f:
            #                 writer = csv.writer(f)
            #                 writer.writerow([trial[2]])
            #         # Get converged FES
            #     else:
            #
            #         if dimension == '1-D Potential':
            #             trial = simulate_1Dsystem(inps, mdps,
            #                                       dimension, method,
            #                                       potfunc, filetitle,
            #                                       makeplot, plot_freq,
            #                                       make_movie)
            #             colvar = pd.DataFrame({'CV': trial[0], 'E': trial[1]})
            #             colvar.reset_index('CV')
            #         else:
            #             trial = simulate_2Dsystem(inps, mdps,
            #                                       dimension, method,
            #                                       potfunc, filetitle,
            #                                       makeplot, plot_freq,
            #                                       make_movie)
            #             colvar = pd.DataFrame({'CV1': trial[0][:, 0],
            #                                    'CV2': trial[0][:, 1],
            #                                    'E': trial[1]})
            #             colvar.reset_index('CV1')
            #         colvar.index.name = 'Step'
            #         if method == "MD":
            #             colvar.to_csv(filetitle+'_MD.csv')
            #         else:
            #             colvar.to_csv(filetitle+'_COLVAR.csv')
            #         with open(filetitle+'info.csv', "wb") as f:
            #             writer = csv.writer(f)
            #             writer.writerow(['RMSD', 'RMSDKLD',
            #                             'RMSDAligned'])
            #             writer.writerow([trial[2]])
            #             writer.writerow([trial[3]])
            #         with open(filetitle+'info.csv', "wb") as f:
            #             writer = csv.writer(f)
            #             writer.writerow([trial[0]])
            #             writer.writerow([trial[1]])



            if cbsp.isChecked():
                makeplot = 'True'
                plot_freq = float(textboxguf.text())
                ebound = np.array([float(textboxmine.text()),
                                   float(textboxmaxe.text())])
                if cb.isChecked():
                    make_movie = 'True'
                else:
                    make_movie = 'False'
            else:
                makeplot = 'False'
                make_movie = 'False'
                plot_freq = 0.0
                ebound = np.array([0,1])
            method = str(combo2.currentText())
            potential = str(combo.currentText())
            potential_dictionary = potential_class.get_potential_dict()
            potential_type = potential_dictionary[potential]
            potential = potential_type()
            new_parameters = np.zeros(len(potential.parameters))
            if len(new_parameters) == 2:
                try:
                    new_parameters[0] = float(textboxpp1.text())
                except:
                    new_parameters[0] = potential.parameters[0]
                try:
                    new_parameters[1] = float(textboxpp2.text())
                except:
                    new_parameters[1] = potential.parameters[1]
            if len(new_parameters) == 4:
                try:
                    new_parameters[0] = float(textboxpp1.text())
                except:
                    new_parameters[0] = potential.parameters[0]
                try:
                    new_parameters[1] = float(textboxpp2.text())
                except:
                    new_parameters[1] = potential.parameters[1]
                try:
                    new_parameters[2] = float(textboxpp1.text())
                except:
                    new_parameters[2] = potential.parameters[2]
                try:
                    new_parameters[3] = float(textboxpp2.text())
                except:
                    new_parameters[3] = potential.parameters[3]
            if len(potential.rare_event) == 2:
                try:
                    potential.set_lower_rare_event(float(textboxRE1.text()))
                except:
                    potential.set_lower_rare_event(potential.rare_event[0])
                try:
                    potential.set_upper_rare_event(float(textboxRE2.text()))
                except:
                    potential.set_upper_rare_event(potential.rare_event[1])
            elif len(potential.rare_event) == 4:
                try:
                    potential.set_lower_rare_event(float(textboxRE1.text()),
                                                   float(textboxRE2.text()))
                except:
                    potential.set_lower_rare_event(potential.rare_event[0][0],
                                                   potential.rare_event[0][1])
                try:
                    potential.set_upper_rare_event(float(textboxRE3.text()),
                                                   float(textboxRE4.text()))
                except:
                    potential.set_upper_rare_event(potential.rare_event[1][0],
                                                   potential.rare_event[1][1])

            potential.set_parameters(new_parameters)
            units = str(combo3.currentText())
            if units == "Kcal/mol" or units == "Energy Units":
                kb = 0.00198
            elif units == "KJ/mol":
                kb = 0.008314
            elif units == "Joule/mol":
                kb = 8.314
            elif units == "cal/mol":
                kb = 1.98
            dimension = potential.dimension
            filetitle = str(textboxfn.text())
            mdps = np.zeros(5)
            bc_dict = boundarycondition.get_boundary_condition_dict()

            if (str(combo4.currentText()) == 'Periodic' and
                  bc_actx.isChecked()):
                 xbc = bc_dict[str(combo4.currentText())]
                 xbc = xbc()
                 xbc.set_new_upper_location(float(textbox32.text()))
                 xbc.set_new_lower_location(float(textbox31.text()))
            elif str(combo4.currentText()) == 'Quartic':
                 xbc = bc_dict[str(combo4.currentText())]
                 xbc = xbc()
                 if cbxbchigh.isChecked():
                    xbc.set_new_upper_location(float(textbox32.text()))
                 if cbxbclow.isChecked():
                    xbc.set_new_lower_location(float(textbox31.text()))
            else:
                xbc = bc_dict['Empty']
                xbc = xbc()
            bcs = [xbc]
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
                if (str(combo4.currentText()) == 'Periodic' and
                      bc_acty.isChecked()):
                     ybc = bc_dict[str(combo4.currentText())]
                     ybc = ybc()
                     ybc.set_new_upper_location(float(textboxy2.text()))
                     ybc.set_new_lower_location(float(textboxy1.text()))
                elif str(combo4.currentText()) == 'Quartic':
                    if cbxbchigh.isChecked():
                        ybc.set_new_upper_location(float(textboxy2.text()))
                    if cbxbclow.isChecked():
                        ybc.set_new_lower_location(float(textboxy1.text()))
                else:
                    ybc = bc_dict['Empty']
                    ybc = ybc()
                bcs.append(ybc)
            if (method == 'Metadynamics'):
                mdps = np.array([float(textbox6.text()),
                                float(textbox7.text()),
                                float(textbox8.text())])
            elif method == 'MD':
                mdps = np.array([0.0,0.01,inps[0]*2,1.0,1.0])
            elif (method == 'Well-Tempered Metadynamics'):
                biasfactor = float(textbox9.text())
                DT = (biasfactor-1)*float(textbox4.text())
                mdps = np.array([float(textbox6.text()),
                                float(textbox7.text()),
                                float(textbox8.text()),
                                DT,
                                1.0])
            elif (method == "Infrequent WT MetaD"):
                biasfactor = float(textbox9.text())
                DT = (biasfactor-1)*float(textbox4.text())
                mdps = np.array([float(textbox6.text()),
                                float(textbox7.text()),
                                float(textbox8.text()),
                                DT,
                                float(textbox10.text())])
                try:
                    potential.set_lower_rare_event(float(textboxRE1.text()))
                except:
                    potential.set_lower_rare_event(float("inf")*-1)
                try:
                    potential.set_upper_rare_event(float(textboxRE2.text()))
                except:
                    potential.set_upper_rare_event(float("inf"))
            # Infrequent Metadynamics Loop
            trials = mdps[-1]
            checkprogress = 0
            potfunc = potential
            while checkprogress < trials:
                if potfunc.dimension == '1-D Potential':
                    trial = simulate_1Dsystem(inps, mdps, method, potfunc, bcs,
                                              filetitle, makeplot, plot_freq,
                                              make_movie, ebound)
                else:
                    trial = simulate_2Dsystem(inps, mdps, method, potfunc, bcs,
                                              filetitle, makeplot, plot_freq,
                                              make_movie, ebound)

                if method == 'Infrequent WT MetaD':
                    if checkprogress == 0:
                        timedata = pd.DataFrame({'Time': [trial[0]],
                                                 'Teff': [trial[1]],
                                                 'Event': [trial[3]]})
                        checkprogress = len(timedata)
                    else:
                        newdata = pd.DataFrame({'Time': [trial[0]],
                                                'Teff': [trial[1]],
                                                'Event': [trial[3]]})
                        timedata = timedata.append(newdata, ignore_index=True)
                        checkprogress = len(timedata)
                else:

                    if potfunc.dimension == '1-D Potential':
                        colvar = pd.DataFrame({'CV': trial[1][0],
                                               'E': trial[1][1]})
                        colvar.reset_index('CV')
                    else:
                        colvar = pd.DataFrame({'CV1': trial[0][:, 0],
                                               'CV2': trial[0][:, 1],
                                               'E': trial[1]})
                        colvar.reset_index('CV1')
                    colvar.index.name = 'Step'
                    colvar.to_csv(filetitle+'_COLVAR.csv')
                    with open(filetitle + '_info.csv', "ab") as f:
                            writer = csv.writer(f)
                            writer.writerow(['RMSD', 'RMSDkld',
                                             'RMSD alignerr'])
                            writer.writerow([trial[2]])
                            writer.writerow([trial[3]])
                    break
            if os.path.isfile(filetitle + '_info.csv') is False:
                with open(filetitle + '_info.csv', "ab") as f:
                            writer = csv.writer(f)
                            writer.writerow([trial[-2]])

            if method == 'Infrequent WT MetaD':
                timedata.to_csv(filetitle + '_Allevents.csv', delimiter=',')
                ks_results = perform_ks_analysis(timedata)
                if len(timedata[timedata['Event'] == 'A']) > 0:
                    ks_resultsA = perform_ks_analysis(timedata[timedata['Event'] == 'A'])
                if len(timedata[timedata['Event'] == 'B']) > 0:
                    ks_resultsB = perform_ks_analysis(timedata[timedata['Event'] == 'B'])
                # if os.path.isfile('bootstrapped.csv') is False:
                #     with open('bootstrapped.csv', "ab") as f:
                #             writer = csv.writer(f)
                #             writer = writer.writerow(['Means', 'Pvals', 'Rejected'])
                # while monitor <= 1000:
                #     (means, pvals, reject) = sampling(filetitle + '_Allevents.csv', 1000,
                #                                       round(len(timedata)/2))
                #     with open('bootstrapped.csv', "ab") as f:
                #             writer = csv.writer(f)
                #             writer.writerow([means, pvals, reject])
                #     checkprogress = pd.read_csv('bootstrapped.csv')
                #     checkaccept = checkprogress[checkprogress['Rejected'] == 'No']
                #     if ((len(checkprogress) - len(checkaccept))/len(checkprogress) >
                #        0.90 and len(checkprogress) > 100):
                #         break
                #     monitor = len(checkaccept)
                #
                # finisheddata = pd.read_csv('bootstrapped.csv')
                # validdata = finisheddata[finisheddata['Rejected'] == 'No']
                # rejectedtrials = (len(finisheddata) - len(validdata))
                if os.path.isfile(filetitle + '_statistics.csv') is False:
                    with open(filetitle + '_statistics.csv', "ab") as f:
                            writer = csv.writer(f)
                            # writer.writerow(['Bootstrapped: Mean Escape Time',
                            #                  'Mean p-value', '# of Trials Rejected'])
                            # writer.writerow([validdata['Means'].mean(),
                            #                  validdata['Pvals'].mean(),
                            #                  rejectedtrials])
                            writer.writerow(['All events'])
                            writer.writerow([ks_results])
                            if len(timedata[timedata['Event'] == 'A']) > 0:
                                writer.writerow(['A events'])
                                writer.writerow([ks_resultsA])
                            if len(timedata[timedata['Event'] == 'B']) > 0:
                                writer.writerow(['B events'])
                                writer.writerow([ks_resultsB])


        btn.clicked.connect(on_click)
        helpbtn.clicked.connect(showhelp)
app = QApplication(sys.argv)
form = Form()
form.show()
app.exec_()
