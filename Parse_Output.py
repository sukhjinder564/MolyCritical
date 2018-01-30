##############################################################################################################################################################################
#                                                            This script parses the output file for the HEU                                                                  #
#                                                          critical experiment for Mo cross section validation.                                                              #
##############################################################################################################################################################################
#                                                                                                                                                                            #
# parse() flags:                                                                                                                                                             #
# 1 - use for single output file                                                                                                                                             #
# 2 - use in optimization script to return max sensitivity in 10 eV - 100 keV                                                                                                #
# 3 - use in optimization script to return keff                                                                                                                              #
#                                                                                                                                                                            #
##############################################################################################################################################################################

import os
import argparse
import subprocess
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import ndimage

# ==========================================================================> Main Parse Function <========================================================================= #       
    
def Parse(fname, flag):

    # ========================================================================> Global Variables <========================================================================== #       

    global isotopes                                                                                                    # Global variable - list of molybdenum isotope masses
    global abundances                                                                                                  # Global variable - list of Mo isotope abundances
                                                                                                                       
    isotopes = ['92', '94', '95', '96', '97', '98', '100']
    abundances = [0.1453, 0.0915, 0.1584, 0.1667, 0.0960, 0.2439, 0.0982]

    # =================================================================> Extract Multiplication Factor <==================================================================== #       

    def get_k(line):                                                                                                   # Extracts mutliplication factor from output file
        return [float(line.split()[2]), float(line.split()[3])]

    # ==================================================================> Extract Kinetics Parameters <===================================================================== #       

    def get_kinetics(i, lines):                                                                                        # Extracts kinetics parameters from output file
        return [float(lines[i].split()[2]), float(lines[i].split()[3])], [float(lines[i+1].split()[1]), float(lines[i+1].split()[2])], [float(lines[i+2].split()[1]), float(lines[i+2].split()[2])]
        
    # ===================================================================> Extract Fission Tally Data <===================================================================== #       

    def get_tally(i, lines):                                                                                           # Extracts fission tally data from output file
        tally = []
        
        for j in range(i+3, 1000000):
            if lines[j] == '\n':
                break
            else:
                tally.append(float(lines[j].split()[1]))
                
        return tally
        
    # ===============================================================> Extract Isotopic Sensitivity Data <================================================================== #       
        
    def get_sensblocks(sensstart, sensend, lines):                                                                     # Extracts individual isotope sensitivity data from
        sens, senserror, blocks = [], [], []                                                                           # output file
        
        senslines = [x for x in lines[sensstart:sensend] if x != ' \n']
        headeridx = [senslines.index(y) for y in senslines if '.80c' in y]
        
        for i in range(0, len(headeridx)-1):
            blocks.append(senslines[headeridx[i]+2:headeridx[i+1]])
            
        for j in range(0, len(blocks)):
            tempsens, temperror = [], []
            for line in blocks[j]:
                tempsens.append(float(line.split()[2]))
                temperror.append(float(line.split()[3]))
            sens.append(tempsens)
            senserror.append(temperror)
            
        return sens, senserror
                
    # ===================================================> Calculate Fraction of Flux in Intermediate Energy Region <======================================================= #       
        
    def calc_fluxfrac():                                                                                               # Calculates the fraction of the flux in the
        try:                                                                                                           # intermediate energy region
            data = pd.read_csv('meshtal_24', header=None, delim_whitespace=True, dtype=str)
            data.columns = ['Bin', 'x', 'y', 'z', 'Tally', 'Error']
            
            totalfluxdata = data[data.Bin == 'Total'].Tally.astype(float)
            thermalfluxdata = data[data.Bin == '1.000E-05'].Tally.astype(float)
            interfluxdata = data[data.Bin == '1.000E-01'].Tally.astype(float)
            fastfluxdata = data[data.Bin == '2.000E+01'].Tally.astype(float)
            
            tff = sum(thermalfluxdata)/sum(totalfluxdata)
            iff = sum(interfluxdata)/sum(totalfluxdata)
            fff = sum(fastfluxdata)/sum(totalfluxdata)
            
            return tff, iff, fff        
        except:
            raise Exception('meshtal file not found!')
            
    # =================================================> Calculate Fraction of Fissions in Intermediate Energy Region <===================================================== #       
            
    def calc_tallyfrac(tally):                                                                                          # Calculates the fraction of fissions occuring in the 
        tff = tally[0]/tally[-1]
        iff = tally[1]/tally[-1]
        fff = tally[2]/tally[-1]
        
        return tff, iff, fff
        
    # ===============================================> Calculate Natural XS Sensitivity in Intermediate Energy Region <===================================================== #       
        
    def calc_naturalsens(sens, error):                                                                                 # Calculates the system sensitivity to natural isotopic
        absoluteerror = np.zeros_like(error)                                                                           # molybdenum cross sections
        
        for i in range(0, len(error)):
            for j in range(0, len(error[i])):
                absoluteerror[i][j] = float(sens[i][j]) * float(error[i][j])
                
        naturalsens, temperror = np.zeros_like(sens[0]), np.zeros_like(absoluteerror[0])
        
        for k in range(0, len(sens)):
            naturalsens += abundances[k]*np.array(sens[k])
            
        for l in range(0, len(abundances)):
            for m in range(0, len(absoluteerror)):
                temperror += (abundances[l]**2)*(np.array(absoluteerror[m])**2)
                
        naturalerror = temperror**0.5
        
        return naturalsens, naturalerror
        
    # ===========================================> Calculate Maximum Isotopic XS Sensitivity in Intermediate Energy Region <================================================ #       
    
    def calc_sensmax(sens):                                                                                            # Calculates the maximum intermediate energy sensitivity
        maxsens = max(abs(sens))
        return maxsens
        
    # =======================================================================> Plot Sensitivity <=========================================================================== #       

    def plot_sensitivity(sens, error):                                                                                 # Plot sensitivity as a function of energy 
        bins = np.logspace(-5, -1, 202)
        naturalsens, naturalerror = [], []
        
        for i in range(0, len(sens)):
            naturalsens.append(calc_naturalsens(sens[i], error[i])[0])
            naturalerror.append(calc_naturalsens(sens[i], error[i])[1])
                    
        fig, ax = plt.subplots()
        
        for j in range(0, len(naturalsens)):
            if j == 0:
                label = 'Capture'
            elif j == 1:
                label = 'Elastic Scattering'
                
            ax.step(bins[:-1], naturalsens[j], where='post', label=label)
        
        ax.set_title('Natural Isotopic Molybdenum Cross Section Sensitivity', y=1.025)
        ax.set_ylabel('Sensitivity')
        ax.set_xlabel('Energy [MeV]')
        ax.set_xscale('log')
        ax.legend()
        ax.grid()
            
        return ax
    
    # ======================================================================> Split meshtal File <========================================================================== #

    def split_meshtal():
        with open('meshtal') as ofile:
            lines = ofile.readlines()
        
        startidx = []
        
        for i in range(0, len(lines)):
            if 'Mesh Tally Number' in lines[i]:
                startidx.append(i)
                
        for flag in ['24', '34', '44']:
            with open('meshtal_' + flag, 'w') as wfile:
                if flag == '24':
                    writelines = lines[startidx[0]+11:startidx[1]-1]
                elif flag == '34':
                    writelines = lines[startidx[1]+11:startidx[2]-1]
                elif flag == '44':
                    writelines = lines[startidx[2]+11:-1]
                  
                for j in range(0, len(writelines)):
                    wfile.write(writelines[j])
    
    def remove_meshtal():
        for flag in ['24', '34', '44']:
            try:
                os.remove(os.getcwd() + '\\' + 'meshtal_' + flag)
            except OSError:
                pass
    
    # =======================================================================> Plot Fmesh Tally <=========================================================================== #
            
    def plot_fmesh(flag):
    
        if flag == '24':
            axes = []
            binnums = [3,2,1]
            
            for bin in binnums:
                with open('commands.txt', 'w') as wfile:
                    wfile.write('FMESH 24\nBASIS 0 1 0 0 0 1\nEBIN {0}\nFILE\nEND\nEND'.format(bin))
                
                with open('Fmeshplot.bat', 'w') as wfile:
                    wfile.write('FOR %%a IN (*.inp.r) do (\n')
                    wfile.write('   start /b /low /wait C:\MCNP\MCNP_CODE\\bin\mcnp611.exe z RUNTPE=%%a COM=commands.txt )\n')
                    wfile.write('gswin64c -o FMESH{0}.png -sDEVICE=pngalpha plotm.ps'.format(bin))
        
                p1 = subprocess.Popen('env.bat', cwd=os.getcwd())
                stdout, stderr = p1.communicate()
        
                p2 = subprocess.Popen("Fmeshplot.bat", cwd=os.getcwd())
                stdout, stderr = p2.communicate()
        
                img = ndimage.imread('FMESH{0}.png'.format(bin))
                rotated = ndimage.rotate(img, 270)
            
                fig, ax = plt.subplots()
                
                if bin == 1:
                    title = 'Thermal (0 - 10 eV) Flux Shape'
                elif bin == 2:
                    title = 'Intermediate (10 eV - 100 keV) Flux Shape'
                else:
                    title = 'Fast (100 keV - 20 MeV) Flux Shape'
                  
                ax.imshow(rotated)
                ax.set_title(title)
                ax.axis('off')
                
                delfiles = ['comout', 'outp', 'FMESH{0}.png'.format(bin), 'plotm.ps', 'commands.txt', 'Fmeshplot.bat']
                for file in delfiles:
                    try:
                        os.remove(os.getcwd() + '\\' + file)
                    except OSError:
                        pass
                        
                axes.append(ax)
            
            return axes
            
        elif flag == '34' or flag == '44':
                    
            data = pd.read_csv('meshtal_' + flag, header=None, delim_whitespace=True)
            data.columns = ['Energy', 'R', 'Z', 'Theta', 'Tally', 'Error']
            data = data[data.Energy != 'Total']
            
            thermal = data[data.Energy == '1.000E-05']
            inter = data[data.Energy == '1.000E-01']
            fast = data[data.Energy == '2.000E+01']
            
            fig, ax = plt.subplots()
            
            if flag == '34':
                ax.plot(thermal.R, thermal.Tally, label='Thermal', color = 'r')
                ax.plot(inter.R, inter.Tally, label='Intermediate', color = 'g')
                ax.plot(fast.R, fast.Tally, label='Fast', color = 'b')
                ax.set_title('Radial Grouped Flux Profile')
                ax.set_ylabel('Flux')
                ax.set_xlabel('Radial Position [cm]')
            elif flag == '44':
                ax.plot(thermal.Z, thermal.Tally, label='Thermal', color='r')
                ax.plot(inter.Z, inter.Tally, label='Intermediate', color='g')
                ax.plot(fast.Z, fast.Tally, label='Fast', color='b')
                ax.set_title('Axial Grouped Flux Profile')
                ax.set_ylabel('Flux')
                ax.set_xlabel('Axial Position Relative to Bottom [cm]')
                
            ax.grid()
            ax.legend()
                
            return ax
            
    # =========================================================> Calculate Stationary and Mobile Unit Weights <============================================================= #
    # Stationary unit weight - weight of stationary unit cells + reflector
    # Mobile unit weight - weight of mobile unit cells
    
    def calc_weights(lines):
        Inventory = {'solid':5, '2.5':7, '6':7}
        UsedPlates = {'solid':0, '2.5':0, '6':0}
        UpperPlateTypes = ['solid', '2.5', '6']
        LowerPlateTypes = ['2.5', '6']
        FuelThick = 0.29972
        FuelDen = 18.8525
        MolyDen = 10.22
    
        for i in range(0, len(lines)):
            if 'Moderator Thickness' in lines[i]:
                ModThick = float(lines[i].split()[5])
            elif 'Moderator Density' in lines[i]:
                ModDen = -1*float(lines[i].split()[5])
            elif 'Molybdenum Thickness' in lines[i]:
                MolyThick = float(lines[i].split()[5])
            elif 'Reflector Density' in lines[i]:
                RefDen = -1*float(lines[i].split()[5])
            elif 'Number of Moderator' in lines[i]:
                ModNum = int(lines[i].split()[10])
            elif 'Number of Fuel' in lines[i]:
                FuelNum = int(lines[i].split()[10])
            elif 'Number of Molybdenum' in lines[i]:
                MolyNum = int(lines[i].split()[10])
            elif 'Number of Top' in lines[i]:
                TopUnits = int(lines[i].split()[7])
            elif 'Number of Bottom' in lines[i]:
                BotUnits = int(lines[i].split()[7])
            elif 'Corner Reflector Bottom' in lines[i]:
                CRBot = float(lines[i].split()[6])
            elif 'Corner Reflector Top' in lines[i]:
                CRTop = float(lines[i].split()[6])
        
        def calc_plate_weight(loc, type, hole, den, thick):
            OuterR = 26.67
            
            if loc == 'bottom' and type != 'fuel':
                InnerR = (2.48/2)*2.54
            elif loc == 'top' and type != 'fuel':
                InnerR = 0
            else:
                if hole == 'solid':
                    InnerR = 0
                elif hole == '2.5':
                    InnerR = (2.5/2) * 2.54
                elif hole == '6':
                    InnerR = (6/2) * 2.54
                
            V = np.pi*thick*(OuterR**2 - InnerR**2)
            m = den * V
            return m
            
        def calc_ref_weight(den, CRBot, CRTop):
        
            TopRefV = ((22*22*5.68) - np.pi*((0.25/2)**2)*5.68)*(2.54**3)
            BotRefV = np.pi*(10.5**2 - 1.24**2)*5.68  *(2.54**3)
            SideRefV = ((34.76*34.76*48.78) - (22*22*48.78)) * (2.54**3)
            CRefV = ((22*22*((CRTop-CRBot)/2.54))-(np.pi*10.5*2*((CRTop-CRBot)/2.54))) * (2.54**3)
            
            TotalV = TopRefV + BotRefV + SideRefV + CRefV

            m = den*TotalV
            return m
        
        TopFuelWt = 0
        BotFuelWt = 0
        
        for type in UpperPlateTypes:
            if Inventory[type] < UsedPlates[type]:
                TopFuelWt += calc_plate_weight('top', 'fuel', type, FuelDen, FuelThick) * FuelNum
            else:
                pass
                
        for type in LowerPlateTypes:
            if Inventory[type] < UsedPlates[type]:
                BotFuelWt += calc_plate_weight('bottom', 'fuel', type, FuelDen, FuelThick) * FuelNum
            else:
                pass
        
        TopModWt = calc_plate_weight('top', 'mod', 'solid', ModDen, ModThick) * ModNum
        TopMolyWt  = calc_plate_weight('top', 'moly', 'solid', MolyDen, MolyThick) * MolyNum
        TopUCellWt = TopFuelWt + TopModWt + TopMolyWt
        
        BotModWt = calc_plate_weight('bottom', 'mod', 'solid', ModDen, ModThick) * ModNum
        BotMolyWt = calc_plate_weight('bottom', 'moly', 'solid', MolyDen, MolyThick) * MolyNum
        BotUCellWt = BotFuelWt + BotModWt + BotMolyWt
        
        RefWt = calc_ref_weight(RefDen, CRBot, CRTop)
        
        StatWeight = (TopUnits*TopUCellWt + RefWt) * (1/453.592)
        MobWeight = (BotUnits*BotUCellWt) * (1/453.592)
        
        return StatWeight, MobWeight

    # =====================================================================> Check Tally Statistics <======================================================================= #       
        
    def check_tallies(line):                                                                                           # Checks if the tallies passed all statistical checks
        checklist = line.split()[1:]
        
        if 'no' in checklist:
            return False
        else:
            return True
        
    # ======================================================================> Write Run Summary <=========================================================================== #       
        
    def write_summary(data):                                                                                           # Writes the run summary
    
        def align(line):
            while len(line) < 78:
                line += ' ' 
            line += '#'
            return line        
        
        print('\n')
        print('###############################################################################')
        print('#                                                                             #')
        print('# ==============================> RUN SUMMARY <============================== #')
        print('#                                                                             #')
        print('###############################################################################')
        print('#                                                                             #')
        print('# This summary provides detailed information concerning the run, namely the   #')
        print('# system eigenvalues and kinetics paramters as well as the fraction of the    #')
        print('# flux and fission rate occuring in the intermediate (10 eV - 100 keV) range  #')
        print('#                                                                             #')
        print('###############################################################################')
        print('#                                                                             #')
        print(align('# Multiplication Factor:                {0} +/- {1}'.format(data['keff'][0], data['keff'][1])))
        print(align('# Mean Generation Time:                 {0} +/- {1} [ns]'.format(data['gt'][0], data['gt'][1])))
        print(align('# Rossi-Alpha Eigenvalue:               {0} +/- {1} [/ns]'.format(data['ra'][0], data['ra'][1])))
        print(align('# Effective Delayed Neutron Fraction:   {0} +/- {1}'.format(data['beff'][0], data['beff'][1])))
        print('#                                                                             #')
        print(align('# Fraction of flux between 0 and 10 eV:            {0}'.format(round(data['thermalfluxfrac'], 5))))
        print(align('# Fraction of flux between 10 eV and 100 keV:      {0}'.format(round(data['interfluxfrac'], 5))))
        print(align('# Fraction of flux between 100 keV and 20 MeV:     {0}'.format(round(data['fastfluxfrac'], 5))))
        print('#                                                                             #')
        print(align('# Fraction of fissions between 0 and 10 eV:        {0}'.format(round(data['thermalfissionfrac'], 5))))
        print(align('# Fraction of fissions between 10 eV and 100 keV:  {0}'.format(round(data['interfissionfrac'], 5))))
        print(align('# Fraction of fissions between 100 keV and 20 MeV: {0}'.format(round(data['fastfissionfrac'], 5))))
        print('#                                                                             #')
        print(align('# Intermediate Molybdenum Capture Rate:            {0}'.format(round(data['mocaprate'], 5))))
        print('#                                                                             #')
        print(align('# Maximum Mo capture sensitivity:                  {0}'.format(round(data['capturesensmax'], 5))))
        print(align('# Maximum Mo elastic scattering sensitivity:       {0}'.format(round(data['escattersensmax'], 5))))
        print('#                                                                             #')
        print(align('# Weight of stationary portion:                    {0}'.format(round(data['statweight'], 0))))
        print(align('# Weight of mobile portion:                        {0}'.format(round(data['mobweight'],0))))
        print('#                                                                             #')
        print(align('# Source entropy converged?                        {0}'.format(data['converged'])))
        print(align('# Tallies met statistical checks?                  {0}'.format(data['checklist'])))
        print('#                                                                             #')
        print('###############################################################################')
        
    # =========================================================================> Parse Output <============================================================================= #       
        
    def parse_output(lines, flag):                                                                                     # Parses the output file and returns necessary data,
        font = {'size':18}                                                                                             # if applicable
        mpl.rc('font', **font)                                                                                         # Flag = 1: Write summary, plot tallies and
        mpl.rcParams['figure.figsize'] = 9, 7                                                                          # sensitivity
                                                                                                                       # Flag = 2: Returns optimization design variables
        data = {}                                                                                                      
        checklist, sensidx, sens, senserror = [], [], [], []                                            
        tallynum = 1
        converged = False

        for i in range(0, len(lines)):

            if 'final result' in lines[i]:
                data['keff'] = get_k(lines[i])
                
            if 'gen. time' in lines[i]:
                data['gt'], data['ra'], data['beff'] = get_kinetics(i, lines)

            if 'cell union total' in lines[i]:
                if tallynum == 1:
                    ft = get_tally(i, lines)
                    tallynum += 1
                else:
                    mot = get_tally(i, lines)
                           
            if 'passed?' in lines[i]:
                checklist.append(check_tallies(lines[i]))
                
            if 'Source entropy convergence check passed' in lines[i]:
                converged = True
                                
            if 'sensitivity profile' in lines[i]:
                sensidx.append(i)
            if '1average' in lines[i]:
                sensidx.append(i)
        
        if False in checklist:
            data['checklist'] = 'No'
        else:
            data['checklist'] = 'Yes'
            
        if converged == True:
            data['converged'] = 'Yes'
        else:
            data['converged'] = 'No'
        
        for j in range(0, len(sensidx)-1):
            sens.append(get_sensblocks(sensidx[j], sensidx[j+1], lines)[0])
            senserror.append(get_sensblocks(sensidx[j], sensidx[j+1], lines)[1])
        
        data['thermalfissionfrac'], data['interfissionfrac'], data['fastfissionfrac'] = calc_tallyfrac(ft)
        data['capturesensmax'] = calc_sensmax(calc_naturalsens(sens[0], senserror[0])[0])
        data['escattersensmax'] = calc_sensmax(calc_naturalsens(sens[1], senserror[1])[0])                
        data['statweight'], data['mobweight'] = calc_weights(lines)
        
        if flag == 1:
            split_meshtal()
            data['thermalfluxfrac'], data['interfluxfrac'], data['fastfluxfrac'] = calc_fluxfrac()
            data['mocaprate'] = mot[1]
            ax456 = plot_fmesh('24')
            ax4, ax5, ax6 = ax456[0], ax456[1], ax456[2]
            ax3 = plot_fmesh('44')
            ax2 = plot_fmesh('34')
            ax1 = plot_sensitivity(sens, senserror)
            write_summary(data)            
            remove_meshtal()
            plt.show()
        elif flag == 2:
            passed = True
            if data['statweight'] > 20000:
                passed = False
            if data['mobweight'] > 2000:
                passed = False
            if (data['keff'][0] - data['keff'][1]) < 1 or (data['keff'][0] + data['keff'][1]) > (1 + 0.8*data['beff']):
                passed = False
            if data['interfissionfrac'] != max([data['thermalfissionfrac'], data['interfissionfrac'], data['fastfissionfrac']]):
                passed = False
            return data['capturesensmax'], data['escattersensmax'], passed        
        else:
            raise Exception('Incorrect flag specification')
            
    # ===============================================================================> Run <================================================================================ #       
            
    with open(fname) as ofile:
        lines = ofile.readlines()
        
    output = parse_output(lines, flag)
    if output == None:
        pass
    else:
        return output
    
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = 'Parse MCNP Output File')
    parser.add_argument('fname', help = 'Output File Name', metavar = '<fname>')
    args = parser.parse_args()
    
    Parse(args.fname, 1)