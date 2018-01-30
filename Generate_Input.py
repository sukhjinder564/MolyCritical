##############################################################################################################################################################################
#                                                         This script configures the base input file for the HEU                                                             #
#                                                          critical experiment for Mo cross section validation.                                                              #
##############################################################################################################################################################################

import os.path
import argparse
import numpy as np

def Generate(TU, BU, MT, MOT, flag):

    Moderators = {'Graphite':[5, 1.7029], 'Be':[6, 1.848], 'Teflon':[7, 2.25], 'Lucite':[8, 1.19]}                         # Options for moderator materials with densities
    Reflectors = {'Copper':[9, 8.96], 'Nickel':[10, 8.902]}                                                                # Options for reflector materials with densities
    PlateTypes = {'solid':5, '2.5':7, '6':7}                                                                               # Plate types and inventory number
    PlateUniverses = {'solid':7, '2.5':8, '6':9}                                                                           # Plate type universe numbers
    UpperTypes = ['solid', '2.5', '6']                                                                                     # Plate types that can be used in the upper unit
    LowerTypes = ['2.5', '6']                                                                                              # Plate types that can be used in the lower unit

    # ====================================================================> Script Configuration Parameters <=================================================================== #       

    Mod = 'Graphite'                                                                                                       # Moderator material
    Ref = 'Copper'                                                                                                         # Reflector material

    CellOrder = ['mo', 'm', 'f', 'm']                                                                                      # Order of Unit Cell (top down)

    PerCycle = 10000                                                                                                       # Number of particles per kcode cycle
    SkipCycle = 25                                                                                                         # Number of kcode cycles to skip
    TotalCycle = 150                                                                                                       # Total number of kcode cycles

    # =================================================================> Command Line Configurations Parameters <=============================================================== #

    TopUnits = TU                                                                                                          # Number of upper unit unit cells
    BottomUnits = BU                                                                                                       # Number of lower unit unit cells

    ModThick = MT                                                                                                          # Thickness of moderator plates
    MolyThick = MOT                                                                                                        # Thickness of molybdenum plates
    
    # =============================================================================> Constraints <============================================================================== #       

    assert TopUnits + BottomUnits <= 5 + 7 + 7                                                                             # Asserts that the plate inventory is met
    assert BottomUnits <= 14 - (TopUnits - 5)                                                                              # Asserts that the plate inventory is met
    
    # ==============================================================================> Set Values <============================================================================== #       

    FuelThick = 0.118 * 2.54                                                                                               # Thickness of fuel plates

    roundval = 5                                                                                                           # Decimal place to round printed variables to

    CRThick = 5.68 * 2.54                                                                                                  # Corner reflector thickness
    SDThick = 0.104 * 2.54                                                                                                 # Steel diaphragm thickness
    PlatenThick = 1 * 2.54                                                                                                 # Platen thickness
    PAPThick = 1.5 * 2.54                                                                                                  # PAP thickness
    TopRefThick = 17.48                                                                                                    # Top reflector thickness
    BottomRefThick = 5.68 * 2.54                                                                                           # Bottom reflector thickness
    TubeHeight = (25 + 0.25) * 2.54                                                                                        # Alignment tube height
    TubeGap = 0.25 * 2.54                                                                                                  # Size of gap above alignment tube
    SysTop = 200                                                                                                           # Location of top of system
    SysBottom = 93.5788                                                                                                    # Location of bottom of system
    RefHSL = 44.15                                                                                                         # Half side length of cross section of reflector
    RefSL = 88.3                                                                                                           # Full side length of cross section of reflector

    # ==========================================================================> Material Parameters <======================================================================== #       
            
    ModNum = str(Moderators[Mod][0])                                                                                       # Determines MCNP input moderator material number
    ModDen = '-' + str(Moderators[Mod][1])                                                                                 # Determines MCNP input moderator material density

    RefNum = str(Reflectors[Ref][0])                                                                                       # Determines MCNP input reflector material number
    RefDen = '-' + str(Reflectors[Ref][1])                                                                                 # Determines MCNP input reflector material density

    # ==============================================================================> Tally Cards <============================================================================ #     

    # Fuel fission tally card written in write file block
    
    MoPNums = np.append(np.where(np.array(CellOrder) == 'mo')[0], np.where(np.array(CellOrder) == 'mo')[0] + len(CellOrder) + 1)
    MoPNums = np.append(MoPNums, np.where(np.array(CellOrder) == 'mo')[0] + 2*(len(CellOrder) + 1)).tolist()
    MolyPlates = ' '.join([str(int(i+1)) for i in MoPNums])
    MolyTallySD = ' '.join(['1'] * (len(MoPNums) + 1))

    # ===============================================================================> Functions <============================================================================= #       

    def align(line):                                                                                                       # Function that aligns the comment at column 80 in
        data = line.split('$')[0]                                                                                          # the MCNP input decks
        comment = line.split('$')[1]
        
        while len(data) < 79:
            data += ' ' 
        
        return data + '$' + comment
        
    def mat_thick(mat):                                                                                           # Determines plate thickness base on the plate
        if mat=='m':                                                                                                       # material
            return ModThick
        elif mat=='f':
            return FuelThick
        elif mat=='mo':
            return MolyThick
            
    def check_inventory(types, used, inventory, unis):
        for i in range(0, len(types)):
            if used[types[i]] < inventory[types[i]]:
                used[types[i]] += 1
                return unis[types[i]], used
            else:
                pass    
            
    # ======================================================================> Determination of Plate Cards <==================================================================== #       

    AllCells = []                                                                                                          # List of all cell numbers in the input deck
    PlateCount = 1                                                                                                         # Initilize count of plates for the MCNP imp card

    def PlateCells(type, PlateCount, AllCells):                                                                       # Function to generate the mobile and stationary
        Plates = []                                                                                                        # sets of plates
        initplatecount = PlateCount
                
        if type == 'solid':                                                                                                # Defines the universes for units with solid, 2.5"
            typeuniverse = 4                                                                                               # hole, and 6" hole fuel plates
            locuniverse = 1
        elif type == '2.5':
            typeuniverse = 5
            locuniverse = 2
        else:
            typeuniverse = 6
            locuniverse = 3
        
        for plate in CellOrder:                                                                                            # Writes the lines determining the plate order based
                                                                                                                           # on the user-defined input
            if PlateCount < 10:
                gap = '    '
            else:
                gap = '   '
                
            if type == 'solid':  
                BottomSurf = PlateCount
                TopSurf = PlateCount - 1
            elif type == '2.5':
                BottomSurf = PlateCount - initplatecount + 1
                TopSurf = PlateCount - initplatecount 
            elif type == '6':
                BottomSurf = PlateCount - initplatecount + 1
                TopSurf = PlateCount - initplatecount

            if plate == 'm' and PlateCount == initplatecount:
                cellline = '{0}{1}{2} {3}  {4}  U={5}  $ Plate {0}\n'.format(PlateCount, gap, ModNum, ModDen, BottomSurf, locuniverse)
            elif plate == 'm' and PlateCount == initplatecount + len(CellOrder) - 1:
                cellline = '{0}{1}{2} {3}  -{4}  U={5}  $ Plate {0}\n'.format(PlateCount, gap, ModNum, ModDen, TopSurf, locuniverse)
            elif plate == 'm':
                cellline = '{0}{1}{2} {3}  {4} -{5}  U={6}  $ Plate {0}\n'.format(PlateCount, gap, ModNum, ModDen, BottomSurf, TopSurf, locuniverse)

            if plate == 'mo' and PlateCount == initplatecount:
                cellline = '{0}{1}4 -10.22  {2}  U={3}  $ Plate {0}\n'.format(PlateCount, gap, BottomSurf, locuniverse)
            elif plate == 'mo' and PlateCount == initplatecount + len(CellOrder) - 1:
                cellline = '{0}{1}4 -10.22  -{2}  U={3}  $ Plate {0}\n'.format(PlateCount, gap, TopSurf, locuniverse)
            elif plate == 'mo':
                cellline = '{0}{1}4 -10.22  {2} -{3}  U={4}  $ Plate {0}\n'.format(PlateCount, gap, BottomSurf, TopSurf, locuniverse)
                
            if plate == 'f' and PlateCount == initplatecount:
                cellline = '{0}{1}0          {2}  U={3}  FILL={4} $ Plate {0}\n'.format(PlateCount, gap, BottomSurf, locuniverse, typeuniverse)
            elif plate == 'f' and PlateCount == initplatecount + len(CellOrder) - 1:
                cellline = '{0}{1}0         -{2}  U={3}  FILL={4} $ Plate {0}\n'.format(PlateCount, gap, TopSurf, locuniverse, typeuniverse)
            elif plate == 'f':
                cellline = '{0}{1}0          {2} -{3}  U={4}  FILL={5} $ Plate {0}\n'.format(PlateCount, gap, BottomSurf, TopSurf, locuniverse, typeuniverse)
                
            Plates.append(align(cellline))                                                                                 # Adds the plate line to list of unit cell cards

            PlateCount += 1
        
        platenotstring = str(PlateCount) + '   0'                                                                          # MCNP input card defining the outside of unit cell
        
        for cell in range(initplatecount, PlateCount+1):                                                                   # Adds cell numbers to the list of all cell numbers
            AllCells.append(str(cell))
            if cell == PlateCount:
                pass
            else:
                platenotstring += ' #' + str(cell)
                
        platenotstring += ' U={0}  $\n'.format(locuniverse)
        Plates.append(align(platenotstring))                                                                               # Adds the outside of unit cell card to MCNP input
        
        return Plates, PlateCount, platenotstring, AllCells
         
    PlatesSolid, PlateCount, Platenotstringsolid, AllCells = PlateCells('solid', PlateCount, AllCells)        # Determines MCNP block for stationary units
    Plates25, PlateCount, Platenotstring25, AllCells = PlateCells('2.5', PlateCount+1, AllCells)     # Determines MCNP block for mobile units
    Plates6, PlateCount, Platenotestring6, AllCells = PlateCells('6', PlateCount+1, AllCells)

    # ===================================================================> Determination of Number of Cells <=================================================================== #

    AllCells = AllCells + ['90', '91', '92', '93', '94', '95', '96', '97', '100', '101', '102', '103', '104', '105',
                           '110', '111', '112', '113', '120', '200', '201', '202', '203', '204', '1000', '1001']

    # ====================================================================> Calculation of Input Parameters <=================================================================== #

    num_m = CellOrder.count('m')                                                                                           # Counts number of moderator plates in the unit cell
    num_f = CellOrder.count('f')                                                                                           # Counts number of fuel plates in the unit cell
    num_mo = CellOrder.count('mo')                                                                                         # Counts number of molybdenum plates in the unit cell

    UCellThick = round(num_m*ModThick + num_f*FuelThick + num_mo*MolyThick, roundval)                                      # Calculates unit cell total thickness
    SDLoc = round(SysBottom + 4*CRThick, roundval)                                                                         # Calculates steel diaphragm position
    CoreTop = round(SDLoc + SDThick + TopUnits*UCellThick, roundval)                                                       # Calculates position of the top of the core
    MobileBase = round(BottomRefThick + PlatenThick + PAPThick, roundval)                                                  # Calculates the position of the bottom reflector top
    LowerPos = round((CoreTop - TopUnits*UCellThick - BottomUnits*UCellThick - SDThick - MobileBase), roundval)            # Calculates the position of the lower unit bottom
    TubeBottom = round((MobileBase + BottomUnits*UCellThick) - TubeHeight, roundval)                                       # Calculates the position of the alignment tube bottom
    TubeTop = round(TubeHeight + TubeBottom - TubeGap, roundval)                                                           # Calculates the position of the alignment tube top
    LowerUnitTop = round(TubeTop + TubeGap, roundval)                                                                      # Calculates the position of the lower unit top
    RefTop = round(SysTop + TopRefThick, roundval)                                                                         # Calculates the position of the reflector top
    RefHeight = round(RefTop - SysBottom, roundval)                                                                        # Calculates the height of reflector
    CellCyl = len(CellOrder)                                                                                               # Calculates the top unit cell plane surface number
    CRTop = CoreTop + 0.75                                                                                                 # Calculates the top of the corner reflector
    TRTop = CRTop + (5.68*2.54)                                                                                            # Calculates the top of the top reflector
    CoreHeight = SysTop - SysBottom

    # ============================================================================> Plate Locations <=========================================================================== #       

    i = 0
    thickness = []
    while i < len(CellOrder):                                                                                              # Determines list of plate thicknesses for surface
        if len(thickness)==0:                                                                                              # positions
            thickness.append(round(UCellThick - mat_thick(CellOrder[i]),roundval))
        else:
            thickness.append(round(thickness[i-1] - mat_thick(CellOrder[i]), roundval))
        i+=1

    # ==============================================================================> Write File <============================================================================== #       

    with open('Base_Input.txt') as ofile:                                                                                  # Open and read base input file
        lines = ofile.readlines()

    fname = '_'.join(['HEUCX', '_'.join(CellOrder), str(ModThick), str(MolyThick), str(TopUnits), str(BottomUnits)]) + '.inp'
        
    with open(fname, 'w') as wfile:                                                                                        # Writes output file with all determined parameters
        
        UCellLocs = []
        UsedPlates = {'solid':0, '2.5':0, '6':0}
        
        for i in range(0, len(lines)):
            if 'ModThick' in lines[i]:
                lines[i] = align(lines[i].replace('ModThick', str(ModThick)))
                
            if 'FuelThick' in lines[i]:
                lines[i] = align(lines[i].replace('FuelThick', str(FuelThick)))
                
            if 'MolyThick' in lines[i]:
                lines[i] = align(lines[i].replace('MolyThick', str(MolyThick)))
                
            if 'UCellThick' in lines[i]:
                lines[i] = align(lines[i].replace('UCellThick', str(UCellThick)))

            if 'ModNum' in lines[i]:
                lines[i] = align(lines[i].replace('ModNum', str(ModNum)))
            
            if 'ModDen' in lines[i]:
                lines[i] = align(lines[i].replace('ModDen', str(ModDen)))
                
            if 'RefNum' in lines[i]:
                lines[i] = align(lines[i].replace('RefNum', str(RefNum)))
            
            if 'RefDen' in lines[i]:
                lines[i] = align(lines[i].replace('RefDen', str(RefDen)))
                
            if 'NumMod' in lines[i]:
                lines[i] = align(lines[i].replace('NumMod', str(num_m)))
            
            if 'NumFuel' in lines[i]:
                lines[i] = align(lines[i].replace('NumFuel', str(num_f)))
            
            if 'NumMoly' in lines[i]:
                lines[i] = align(lines[i].replace('NumMoly', str(num_mo)))
                
            if 'TopUnits' in lines[i]:
                lines[i] = align(lines[i].replace('TopUnits', str(TopUnits)))
            
            if 'BottomUnits' in lines[i]:
                lines[i] = align(lines[i].replace('BottomUnits', str(BottomUnits)))
                
            if 'CoreTop' in lines[i]:
                lines[i] = align(lines[i].replace('CoreTop', str(CoreTop)))
            
            if 'RefTop' in lines[i]:
                lines[i] = align(lines[i].replace('RefTop', str(RefTop)))
                
            if 'LowerPos' in lines[i]:
                lines[i] = align(lines[i].replace('LowerPos', str(LowerPos)))
                
            if 'SDPos' in lines[i]:
                lines[i] = align(lines[i].replace('SDPos', str(SDPos)))
                
            if 'TubeBottom' in lines[i]:
                lines[i] = align(lines[i].replace('TubeBottom', str(TubeBottom)))
                
            if 'TubeTop' in lines[i]:
                lines[i] = align(lines[i].replace('TubeTop', str(TubeTop)))
                
            if 'LowerUnitTop' in lines[i]:
                lines[i] = align(lines[i].replace('LowerUnitTop', str(LowerUnitTop)))
                
            if 'LwrPosRnd' in lines[i]:
                lines[i] = align(lines[i].replace('LwrPosRnd', str(int(LowerPos))))
                
            if 'SDLoc' in lines[i]:
                lines[i] = align(lines[i].replace('SDLoc', str(SDLoc)))
                
            if 'SysTop' in lines[i]:
                lines[i] = align(lines[i].replace('SysTop', str(SysTop)))
            
            if 'SysBottom' in lines[i]:
                lines[i] = align(lines[i].replace('SysBottom', str(SysBottom)))
                
            if 'CoreHeight' in lines[i]:
                lines[i] = align(lines[i].replace('CoreHeight', str(CoreHeight)))
            
            if 'RefHSL' in lines[i]:
                lines[i] = align(lines[i].replace('RefHSL', str(RefHSL)))
                
            if 'RefSL' in lines[i]:
                lines[i] = align(lines[i].replace('RefSL', str(RefSL)))
                
            if 'RefHeight' in lines[i]:
                lines[i] = align(lines[i].replace('RefHeight', str(RefHeight)))
                
            if 'CRTop' in lines[i]:
                lines[i] = align(lines[i].replace('CRTop', str(CRTop)))
                
            if 'TRTop' in lines[i]:
                lines[i] = align(lines[i].replace('TRTop', str(TRTop)))
                
            if 'CellCyl' in lines[i]:
                lines[i] = align(lines[i].replace('CellCyl', str(CellCyl)))
                
            if 'WRITE SOLID PLATES HERE' in lines[i]:
                for line in PlatesSolid:
                    wfile.write(line)
                    
            if 'WRITE 2.5 PLATES HERE' in lines[i]:
                for line in Plates25:
                    wfile.write(line)
            
            if 'WRITE 6 PLATES HERE' in lines[i]:
                for line in Plates6:
                    wfile.write(line)
                    
            if 'WRITE UNIT CELLS HERE' in lines[i]:
                LowerUCellLines = []
                UpperUCellLines = []
                
                statcellnum = 131
                statnotstring = ''
                for j in range(0, TopUnits):
                    uni, UsedPlates = check_inventory(UpperTypes, UsedPlates, PlateTypes, PlateUniverses)
                    AllCells = [str(statcellnum)] + AllCells
                    #statcellloc = round((CoreTop - (j+1)*UCellThick), roundval)
                    statcellloc = round(CoreTop - TopUnits*UCellThick + (j)*UCellThick, roundval)
                    UpperUCellLines.append(align('{0}  0  99 -{3}  FILL={4}  TRCL=(0 0 {1})  U=11  $ Unit cell {2}\n'.format(statcellnum, statcellloc, j+1, CellCyl, uni)))
                    statcellnum += 1
                    UCellLocs.append(statcellloc)
                    if j == TopUnits-1:
                        for k in range(130, statcellnum):
                            statnotstring += '#{0} '.format(int(k))
                        UpperUCellLines.append(align('{0}  0  {1}  U=11  $ Outside of all units\n'.format(statcellnum, statnotstring)))
                        AllCells = [str(statcellnum)] + AllCells
                
                mobcellnum = 114
                mobnotstring = ''
                for j in range(0, BottomUnits):
                    uni, UsedPlates = check_inventory(LowerTypes, UsedPlates, PlateTypes, PlateUniverses)
                    AllCells = [str(mobcellnum)] + AllCells
                    #mobcellloc = round(MobileBase + j*UCellThick, roundval)
                    mobcellloc = round(MobileBase + (BottomUnits-(j+1))*UCellThick, roundval)
                    LowerUCellLines.append(align('{0}   0  99 -{3} 37  FILL={4}  TRCL=(0 0 {1})  U=10  $ Unit cell {2}\n'.format(mobcellnum, mobcellloc, j+1, CellCyl, uni)))
                    mobcellnum += 1
                    UCellLocs.append(round(mobcellloc + LowerPos, roundval))
                    if j == BottomUnits-1:
                        for k in range(110, mobcellnum):
                            mobnotstring += '#{0} '.format(int(k))
                        LowerUCellLines.append(align('{0}   0  {1}  U=10  $ Outside of lower unit\n'.format(mobcellnum, mobnotstring)))
                        AllCells = [str(mobcellnum)] + AllCells
                
                MidLines = ['C =========================> Upper and Lower Unit <============================\n',
                            '130  0  33 -35  FILL=10 TRCL=(0 0 {0})  U=11                              $ Lower unit\n'.format(LowerPos)]

                UCellLines = LowerUCellLines + MidLines + UpperUCellLines
                
                for line in UCellLines:
                    wfile.write(line)
                        
            if 'WRITE PLATE SURFACES HERE' in lines[i]:
                platenum = 1
                for j in range(0, len(thickness)-1):
                    wfile.write(align('{0}    PZ  {1}  $ Plate #{0}\n'.format(j+1, thickness[j])))
                    platenum += 1
                wfile.write(align('{0}   PZ  {1}   $ Stationary unit cell bottom\n'.format(platenum, UCellThick)))
                            
            if 'WRITE KCODE HERE' in lines[i]:
                wfile.write('KCODE  {0}  1.0  {1}  {2}\n'.format(PerCycle, SkipCycle, TotalCycle))
                ptzs = []
                for j in range(0, len(UCellLocs)):
                    ptzs.append(round(UCellLocs[j] + 2*ModThick + 0.5*FuelThick, roundval))
                    ptzs.append(round(UCellLocs[j] + UCellThick - 2*ModThick - 0.5*FuelThick, roundval))
                for k in range(0, len(ptzs)):
                    ptstring = '  0 16.51 {0}  0 -16.51 {0}  16.51 0 {0}  -16.51 0 {0}\n'.format(round(ptzs[k]))
                    if k == 0:
                        wfile.write('KSRC ' + ptstring)
                    else:
                        wfile.write('     ' + ptstring)
                        
            if 'WRITE KOPTS HERE' in lines[i]:
                wfile.write(align('KOPTS  KINETICS=YES  $ Reactor kinetics card\n'))
                
            if 'WRITE KSEN HERE' in lines[i]:
                wfile.write('KSEN1  XS  ISO=42092.80c 42094.80c 42095.80c 42096.80c 42097.80c 42098.80c\n               42100.80c  RXN=-2  ERG=10E-6 200iLOG 100E-3\n')
                wfile.write('KSEN2  XS  ISO=42092.80c 42094.80c 42095.80c 42096.80c 42097.80c 42098.80c\n               42100.80c  RXN=2  ERG=10E-6 200iLOG 100E-3\n')
                
            if 'FuelPlates' in lines[i]:
                FPNums = []
                if UsedPlates['solid'] != 0:
                    FPNums.append('90')
                if UsedPlates['2.5'] != 0:
                    FPNums.append('92')
                if UsedPlates['6'] != 0:
                    FPNums.append('95')
                FuelPlates = ' '.join(FPNums)
                FuelTallySD = ' '.join(['1'] * (len(FPNums) + 1))
           
                lines[i] = align(lines[i].replace('FuelPlates', FuelPlates))
                
            if 'FuelTallySD' in lines[i]:
                lines[i] = align(lines[i].replace('FuelTallySD', FuelTallySD))
                
            if 'MolyPlates' in lines[i]:
                lines[i] = align(lines[i].replace('MolyPlates', MolyPlates))
            
            if 'MolyTallySD' in lines[i]:
                lines[i] = align(lines[i].replace('MolyTallySD', MolyTallySD))    
                        
            if 'Miscellaneous' in lines[i]:
                wfile.write(lines[i])
                wfile.write(align('IMP:N  1  {0}r  0  $ Cell importance\n'.format(len(AllCells)-2)))

            if 'HERE' not in lines[i] and 'Miscellaneous' not in lines[i]:
                wfile.write(lines[i])     

    if flag == 1:
        pass
    elif flag == 2:
        pass
        # USE FOR OPTIMIZATION

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description = 'Generate MCNP Input File')                                             # Parser definition
    parser.add_argument('Mdt', help='Thickness of moderator plates', metavar='<mdt>')                                      # Parser positional variable - moderator thickness
    parser.add_argument('Mot', help='Thickness of molybdenum plates', metavar='<mot>')                                     # Parser positional variable - molybdenum thickness
    parser.add_argument('TopN', help='Number of top unit cells', metavar='<topn>')                                         # Parser positional variable - number of top unit cells
    parser.add_argument('BotN', help='Number of bottom unit cells', metavar='<botn>')                                      # Parser positional variable - number of bottom unit cells 
    args = parser.parse_args()
    
    Generate(int(args.TopN), int(args.BotN), float(args.Mdt), float(args.Mot), 1)