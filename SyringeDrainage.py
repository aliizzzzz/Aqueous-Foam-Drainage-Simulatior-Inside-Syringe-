#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Created on Mon Nov 15 12:01:53 2021

@author: alizz

Mathematical simulation of sclerosing foam drainage kinetics inside a syringe during injection
'''
# https://towardsdatascience.com/how-to-connect-objects-with-each-other-in-different-situations-with-pythonic-ways-d3aaf4c89553
from IPython import get_ipython
import sys
from time import sleep
import math
import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
import time


class Foam(object):
    ''' Class (composite) Foam taking user input for foam volume. Calculates quasistatic foam
    drainage kinetic during foam injection inside a syringe (uses an instance of the Syringe class).
    Intermediate calculations: liquid length, liquid volume and liquid height.
    Auxilary conditions: maximum liquid height < syringe diameter. '''

    def __init__(self):
        ''' Class constructor. Creates an instance of the Syringe() class, while setting initial
        values of foam volume (VF), syringe length (L) '''
        print('\r')
        print('Quasi Static Drainage Simulator for Sclerosing Foams\n\n')
        self.syrObj = Syringe()
        self.inputValidator = dataValidation()
        self.setFoamHalfTime()
        self.setGasFraction()
        self.setFoamVolume()
        self.syrObj.setFlowrate()
        self.setAccuracy()
        self.setSampleFreq()
        self.VL = [x*self.GFraction for x in self.VF]
        self.L = self.syrObj.V*1000/self.syrObj.getSyringeCrossArea()
        self.setTimeParams()
        self.figs = []

    def setFoamHalfTime(self):
        ''' Invokes the validatedInput function to obtain faom type from user. Sets FHT based on
        user input '''
        self.foamType = self.inputValidator.validatedInput(
            'Enter foam type (D/d) for DSS, (T/t) for Tessari:\n>> ', str, range_=['t', 'd'])
        if self.foamType in ['t']:
            self.foamName = 'Tessari'
            self.FHT = 90.0
        if self.foamType in ['d']:
            self.foamName = 'DSS'
            self.FHT = 160

    def setGasFraction(self):
        ''' Takes user input for foam L:G ratio '''
        self.LG = self.inputValidator.validatedInput('Enter ratio of gas (3,4,5):\n>> ', int,
                                                     range_=[3, 4, 5])
        self.GFraction = 1/(1+self.LG)

    def setFoamVolume(self):
        ''' Takes user input for foam volume in mL '''
        prompt = f'Enter foam volumes (maximum {self.syrObj.V}) in mL (maximum 4 integers, '
        prompt += 'separated with spaces):\n>> '
        self.VF = [float(x) for x in
                   self.inputValidator.validatedInput(prompt, list, max_=self.syrObj.V+1e-10)]

    def setTimeParams(self):
        ''' Sets the simulation time parameters - simulation time (T), and
        sampling frequency (Tfrequency in Hz) - for every given flowrate '''
        self.T = {}
        self.tdata = {}
        for flow in self.syrObj.Q:
            self.T[flow] = {}
            self.tdata[flow] = {}
            for vol in self.VF:
                self.T[flow][vol] = vol/flow
                self.tdata[flow][vol] = np.linspace(0, vol/flow, (int(vol/flow)*self.Tfrequency)+1)

    def setAccuracy(self):
        ''' Sets the number of decimal points for cross-referencing modelled area and liquid area.
        Higher values  yield smoother plots, but may take more time.'''
        self.AccuracyDP = self.inputValidator.validatedInput(
            'Enter no. of decimal points (maximum 5) for central angle estimation:\n>> ', int,
            max_=5+1e-10)

    def setSampleFreq(self):
        ''' Sets the sampling frequency. This defines how many time points are created per second
        of flow time '''
        self.Tfrequency = self.inputValidator.validatedInput(
            'Enter sampling frequency in Hz (maximum 50):\n>> ', int, max_=50+1e-10)

    def printMessage(self):
        ''' Prints parameters of constructed simulation '''
        print('\n\nSetting up...')
        sleep(0.25)
        print(f'\nSimulation Paramters\n{50*"="}')
        print(f'\n{self.syrObj.V} mL syringe class (centric tip) constructed')
        print(f'ID = {self.syrObj.D:.2f} mm')
        print(f'\n1:{self.LG} {self.foamName} foam, half time: {self.FHT} seconds')
        print('Foam volumes: %s mL' % self.VF)
        print('Injection flowrates: %s mL/min\n' % [x*60 for x in self.syrObj.Q])
        sleep(1)

    def getFoamParams(self):
        ''' Prints a dictionary containing foam parameters '''
        foamDict = {'Formulation: ': self.foamName, 'Half-time (s)': self.FHT, 'L:G': '1:'+self.LG}
        return foamDict

    def getLiquidHeight(self, tdata, Q, T, foamLiquidContent, foamVol):
        ''' Returns liquid height in mm2
        Uses an iterative method to solve for the central angle (theta)
        A list of simulated theta values (0,2π) are used to calculated simulated liquid area values
        Simulated liquid area values are rounded to 2 d.p., and are checked against the actual
        liquid areas (2 d.p.) Number of iterations may need to be adjusted
        to find all solutions for theta'''
        iterations = 10000
        print(f'\n\nSolving {Q*60} mL/min case (foam volume = {foamVol} mL):\n{50*"-"}')
        print(f'Using {iterations} iterations to solve for central angle...\n')
        while iterations < 500000:
            areas = [round(x, self.AccuracyDP) for x in self.getLiquidArea(tdata, Q,
                                                                           foamLiquidContent)]
            simulatedThetas = np.linspace(0, 2 * math.pi, iterations)
            thetas = []
            simulatedAreas = []

            for theta in simulatedThetas:
                simArea = round(self.modelArea(theta, self.syrObj.D/2), self.AccuracyDP)
                simulatedAreas.append(simArea)

            for area in areas:
                if area not in simulatedAreas:
                    iterations += 5000
                    sys.stdout.write('\rFailed! Trying %d iterations...  ' % (iterations))
                    sys.stdout.flush()
                    break
                thetas.append(simulatedThetas[simulatedAreas.index(area)])
            else:
                print(f'\n\nSolution found using {iterations} iterations!\n')
                sleep(0.5)
                iterations = 500000
                continue

        liquidHeight = [self.calculateHeight(self.syrObj.D, x) for x in thetas]

        # to simulate liquid exiting the syringe:
        liquidHeight = [self.syrObj.D/2 if x > self.syrObj.D/2 else x for x in liquidHeight]
        tdata = np.linspace(0, T, (int(T)*self.Tfrequency)+1)[0:len(liquidHeight)]
        return (liquidHeight, tdata)

    def getLiquidArea(self, tdata, Q, foamLiquidContent):
        ''' Returns liquid corss-section area in mm2 '''
        area = [x/y for x, y in zip([w * 1000 for w in self.getLiquidVolume(tdata,
                                                                            foamLiquidContent)],
                self.getLiquidLength(tdata, Q))
                if x/y < self.syrObj.getSyringeCrossArea() and x/y > 0]  # mm2
        return area

    def getLiquidLength(self, tdata, Q):
        ''' Returns length of liquid inside the syringe in mm '''
        length = self.L
        lengths = []
        dt = tdata[1] - tdata[0]
        for value in tdata:
            dl = dt*(-Q*1000)/self.syrObj.getSyringeCrossArea()
            length += dl
            lengths.append(length)
        return lengths

    def getLiquidVolume(self, tdata, foamLiquidContent):
        ''' Calculates liquid volume based on foam drainage kinetics. Returns liquid volume in mL '''
        drainedLiquidVol = 0
        volumes = []
        dt = tdata[1] - tdata[0]
        for value in tdata:
            dv = dt*foamLiquidContent/(2*self.FHT)
            foamLiquidContent -= dv
            drainedLiquidVol += dv
            volumes.append(drainedLiquidVol)
        return volumes

    def modelArea(self, theta, radius):
        ''' Models cross section area of syringe occupied by the accumulating drained liquid '''
        area = ((radius ** 2) / 2) * (theta-math.sin(theta))
        return area

    def calculateHeight(self, syringeDiameter, angle):
        ''' Calculates height of liquid in a horizontal syrnge given
        syringeDiameter and central angle '''
        height = (syringeDiameter/2)*(1-math.cos(angle/2))
        return height

    def newFigure(self, Q):
        ''' Creates a new figure window '''
        fig = plt.figure(figsize=[12, 9])
        fig.suptitle(f'1:{self.LG} {self.foamName} Foam\nInjection flowrate: {Q*60} mL/min',
                     weight='bold', fontsize=30)
        self.figs.append(fig)
        ax = plt.gca()
        return fig, ax

    def plotLiquidHeights(self, tdata, hdata, foamVol):
        ''' Plots liquid height results'''
        plt.minorticks_on()
        plt.plot(tdata, hdata, label=self.getDataLabels(foamVol))
        plt.tick_params(axis='both', labelsize=18)
        plt.tick_params(axis='both', which='major', length=7, width=1.7)
        plt.tick_params(axis='both', which='minor', length=4, width=1)
        plt.xlabel('Injection Time (s)', fontsize=20, labelpad=20)
        plt.ylabel('Liquid Height (mm)', fontsize=20, labelpad=15)
        plt.legend(ncol=1, prop={'family': 'monospace', 'size': 18}, loc='upper left',
                   handletextpad=0.8)
        ax = plt.gca()
        ax.grid(True, linestyle=':', linewidth=1.5, which='minor')
        ax.grid(True, linestyle='-', linewidth=1.7, which='major')

    def getDataLabels(self, foamVol):
        ''' returns a list of legend labels for all simulations'''
        # label = r'$\mathregular{V_{Foam}}$: %.1f mL' % (foamVol)
        label = r'$\mathregular{V_{Foam}}$: %s mL' % (str(int(foamVol)).zfill(2))
        return label

    def saveFigs(self):
        print(f'Accuracy: {self.AccuracyDP} Decimal points')
        print(f'Sampling frequency: {self.Tfrequency: .1f} Hz')
        print(f'Simulation time: {time.time()-self.start: .2f} s\n')
        userInput = self.inputValidator.validatedInput(
            'Save plots in current directory (y/n)?\n>> ', str, range_=['', 'y', 'n'])
        if userInput in ['', 'y']:
            if os.path.exists(os.getcwd()+'/Figs/'+f'{self.syrObj.V:.0f}mLSyringe'):
                shutil.rmtree(os.getcwd()+'/Figs/'+f'{self.syrObj.V:.0f}mLSyringe')
            os.mkdir(os.getcwd()+'/Figs/'+f'{self.syrObj.V:.0f}mLSyringe')
            for fig in self.figs:
                title = fig._suptitle.get_text()
                for s in ['/', ':', ' ', '\n']:
                    title = title.replace(s, '-')
                fig.savefig(os.getcwd()+'/Figs/'+f'{self.syrObj.V:.0f}mLSyringe/'+title+'.png',
                            dpi=450, transparent=False)
        else:
            pass

    def simulateInjection(self):
        ''' Simulates drainage kinetics inside a syringe of known diameter and volume.
        Calls getLiquidHeight() which initiates the method resolution of simulation.
        User inputs: list of flowrates, list of foam volumes
        Plots a figure for each flowrate '''
        self.start = time.time()
        self.printMessage()
        self.H = []
        for i, flow in enumerate(self.syrObj.Q):
            fig, ax = self.newFigure(self.syrObj.Q[i])
            self.figs.append(fig)
            tempTs = []
            for j, vol in enumerate(self.VF):
                heightAndTimeData = self.getLiquidHeight(self.tdata[flow][vol], self.syrObj.Q[i],
                                                         self.T[flow][vol],
                                                         self.VL[j], self.VF[j])
                self.H.append(heightAndTimeData[0])
                tempT = heightAndTimeData[1]
                tempTs.append(tempT)
                self.plotLiquidHeights(tempT, self.H[-1], self.VF[j])
            ax.set_ylim(bottom=0)
            ax.set_xlim(left=0)


class Syringe():
    ''' Class Syringe (component), constructs a syringe object with given diameter,
    volume and flowrate '''

    def __init__(self):
        ''' Class constructor, sets default values of volume (V, mL) and
        diameter (D, mm) for a 10 mL leur-lock NormJect syringe (12 mL), and takes
        user input for injection rate in mL/min'''

        self.syringeDiameters = {'20': {'V': 20.0, 'd': 20.10}, '10': {'V': 10.0, 'd': 15.96},
                                 '5': {'V': 5.0, 'd': 12.46}, '3': {'V': 3.0, 'd': 9.83}}
        self.inputValidator = dataValidation()
        self.setSyringe()

    def setSyringe(self):
        ''' Initiates an instance of the Syringe class() '''
        # print('Default syringe: 10 (12) mL NotmJect. Proceed (y/n)?')
        userInput = self.inputValidator.validatedInput(
            'Default syringe: 10 (12) mL NotmJect. Proceed (y/n)?\n>> ', str, range_=['', 'y', 'n'])
        if userInput in ['', 'y']:
            self.V = self.syringeDiameters['10']['V']
            self.D = self.syringeDiameters['10']['d']
        else:
            self.V = self.setSyringeVolume()
            self.D = self.syringeDiameters[str(int(self.V))]['d']
        sleep(0.5)

    def setSyringeVolume(self):
        ''' Method for user to modify syringe volume in mL'''

        self.V = float(self.inputValidator.validatedInput(
            'Enter syringe volume in mL (20, 10, 5, 3 mL):\n>> ', int, range_=[3, 5, 10, 20]))
        return self.V

    def setFlowrate(self):
        ''' Method sets foam injection rate, takes input in mL/min, stores Q in mL/s.
        Instantaneous injetion is assumed to occur at 100 mL/min (maximum allowed value)'''
        self.Q = [float(x)/60 for x in self.inputValidator.validatedInput(
            'Enter injection flowrates in mL/min (maximum 6 integers, separated with spaces):\n>> ',
            list, max_=101)]

    def getSyringeCrossArea(self):
        ''' Method returns corss-section area of syringe in mm2 '''
        return (math.pi*(self.D/2)**2)


class dataValidation():
    ''' Data validation class '''

    def validatedInput(self, prompt, type_, max_=None, range_=None):
        ''' Input validation function '''
        while True:
            try:
                inputVar = input(prompt).lower()
                if type_ == int:
                    inputVar = int(inputVar)
                    if max_ is not None:
                        if inputVar >= max_:
                            raise self.MaxError
                    if range_ is not None:
                        if inputVar not in range_:
                            raise ValueError
                    if inputVar < 0:
                        raise ValueError
                    else:
                        break
                if type_ == str:
                    if inputVar not in range_:
                        raise ValueError
                    else:
                        break
                if type_ == list:
                    if not inputVar:
                        raise ValueError
                    if len(inputVar.split()) > 6:
                        raise ValueError
                    for value in inputVar.split():
                        if max_ is not None:
                            if int(value) >= max_:
                                raise self.MaxError
                        if int(value) < 0:
                            raise ValueError
                    break
            except self.MaxError:
                print(f'\nInput must be less than {round(max_,0)}...')
            except ValueError:
                print('\nInvalid input, try again...')
        if type_ == list:
            return np.sort([float(x) for x in inputVar.split()])
        else:
            return inputVar

    class MaxError(Exception):
        ''' Custom exception class '''
        pass


# %% Main
if __name__ == '__main__':
    get_ipython().magic('clear')  # Clear console
    get_ipython().magic('reset -sf')  # Reset variable space
    plt.style.use('ggplot')
    foamObj = Foam()
    foamObj.simulateInjection()
    foamObj.saveFigs()
