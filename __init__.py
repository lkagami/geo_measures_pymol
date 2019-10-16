'''
G_Measures v.0.8
The "Geometric Measures" script that was developed to carry out geometric analysis on protein structures.

Contributors:

Luciano Porto Kagami, Gustavo Machado das Neves, Luís Fernando Saraiva Macedo Timmers, Rafael Andrade Cáceres and Vera Lucia Eifler-Lima

'''

from __future__ import absolute_import
from __future__ import print_function

# Avoid importing "expensive" modules here (e.g. scipy), since this code is
# executed on PyMOL's startup. Only import such modules inside functions.



'''
See more here: http://www.pymolwiki.org/index.php/Modevectors
'''
from pymol.cgo import *    # get constants
from math import *
from pymol import cmd, preset
import time
import shutil
import tempfile
import math
import subprocess
import errno
import numpy as np
TEMP_PATH = tempfile.mkdtemp()
SHAM_PATH = TEMP_PATH+"/g_sham2.xvg"
TRAJ_PATH = TEMP_PATH+"/trajectory.pdb"
_version_ = str("v.0.8")

def modevectors(first_obj_frame, last_obj_frame, first_state=1, last_state=1, outname="modevectors", head=1.0, tail=0.3, head_length=1.5, headrgb="1.0,1.0,1.0", tailrgb="1.0,1.0,1.0", cutoff=4.0, skip=0, cut=0.5, atom="CA", stat="show", factor=1.0, notail=0):
    """
    Authors Sean Law & Srinivasa
    Michigan State University
    slaw_(at)_msu_dot_edu

    Editor Sacha Yee

    USAGE

    While in PyMOL

    Parameter                Preset            Type    Description
    first_obj_frame          Undefined         String  Object name of the first structure.  The mode vector will propagate from this structure.  Defined by user.
    last_obj_frame           Undefined         String  Object name of the last structure.  The mode vector (arrow head) will end at this structure.  Defined by user.
    first_state              1                 Integer Defines state of first object
    last_state               1                 Integer Defines state of last object
    outname                  modevectors       String  Name of object to store mode vectors in.
    head                     1.0               Float   Radius for the circular base of the arrow head (cone)
    tail                     0.3               Float   Radius for the cylinder of the arrow tail (cylinder)
    head_length              1.5               Float   Length of the arrow head (from the base of the cone to the tip of cone)
    head_rgb                 1.0,1.0,1.0       String  RGB colour for the arrow head.
    tail_rgb                 1.0,1.0,1.0       String  RGB colour for the arrow tail.
    cutoff                   4.0               Float   Skips mode vectors that do not meet the cutoff distance (in Angstroms).
    skip                     0                 Integer Denotes how many atoms to skip.  No arrows will be created for skipped atoms.
    cut                      0.0               Float   Truncates all arrow tail lengths (without disturbing the arrow head) (in Angstroms).
    atom                     CA                String  Designates the atom to derive mode vectors from.
    stat                     show              String  Keeps track and prints statistics (total modevectors, skipped, cutoff).
    factor                   1.0               Float   Multiplies each mode vector length by a specified factor.
                                                       Values between 0 and 1 will decrease the relative mode vector length.
                                                       Values greater than 1 will increase the relative mode vector length.
    notail                   0                 Integer Hides tails and only uses cones (porcupine plot)
    """

    framefirst = cmd.get_model(first_obj_frame, first_state)
    framelast = cmd.get_model(last_obj_frame, last_state)
    objectname = outname
    factor = float(factor)
    arrow_head_radius = float(head)
    arrow_tail_radius = float(tail)
    arrow_head_length = float(head_length)
    cutoff = float(cutoff)
    skip = int(skip)
    cut = float(cut)
    atomtype = atom.strip('"[]()')
    objectname = objectname.strip('"[]()')

    headrgb = headrgb.strip('" []()')
    tailrgb = tailrgb.strip('" []()')
    hr, hg, hb = list(map(float, headrgb.split(',')))
    tr, tg, tb = list(map(float, tailrgb.split(',')))

    version = cmd.get_version()[1]
    arrow = []
    arrowhead = []
    arrowtail = []
    x1 = []
    y1 = []
    z1 = []
    x2 = []
    y2 = []
    z2 = []
    exit_flag = False

##############################################################
#                                                            #
# Define an object called "tail" and store the tail and  a   #
# circular base of the triangle in this object.              #
#                                                            #
##############################################################

    skipcount = 0
    skipcounter = 0
    keepcounter = 0
    atom_lookup = {}
    for atom in framefirst.atom:
        if atom.name == atomtype:
            if skipcount == skip:
                x1.append(atom.coord[0])
                y1.append(atom.coord[1])
                z1.append(atom.coord[2])

                ##########################################
                #                                        #
                # Set atom_lookup for a specific atom    #
                # equal to ONE for the first input set.  #
                # This dictionary will be used as a      #
                # reference for the second set.          #
                #                                        #
                ##########################################

                current_atom = "CHAIN " + atom.chain + " RESID "\
                    + atom.resi + " RESTYPE "\
                    + atom.resn +\
                    " ATMNUM " + str(atom.index)
#               print current_atom
                atom_lookup['current_atom'] = 1

                skipcount = 0
                keepcounter += 1
            else:
                #               print skipcount
                skipcount += 1
                skipcounter += 1

    skipcount = 0
    for atom in framelast.atom:
        if atom.name == atomtype:
            if skipcount == skip:
                x2.append(atom.coord[0])
                y2.append(atom.coord[1])
                z2.append(atom.coord[2])

                #########################################
                #                                       #
                # Get atom information from second set  #
                # and compare with first set.  All      #
                # atoms from this second set MUST be    #
                # found in the first set!  Otherwise,   #
                # the script will exit with an error    #
                # since modevectors can only be created #
                # by calculating values from identical  #
                # sources.                              #
                #                                       #
                #########################################

                current_atom = "CHAIN " + atom.chain + " RESID "\
                    + atom.resi + " RESTYPE "\
                    + atom.resn +\
                    " ATMNUM " + str(atom.index)
#               print current_atom
                if 'current_atom' not in atom_lookup:
                    print("\nError: " + current_atom + " from \""\
                          + last_obj_frame +\
                          " \"is not found in \"" + first_obj_frame + "\".")
                    print("\nPlease check your input and/or selections and try again.")
                    exit_flag = True
                    break

                skipcount = 0
            else:
                skipcount += 1

    if exit_flag == 1:
        ###########################################
        #                                         #
        # Exit script because an atom cannot be   #
        # found in both input files               #
        #                                         #
        ###########################################
        return

    cutoff_counter = 0  # Track number of atoms failing to meet the cutoff

    ###################################################
    #                                                 #
    # Check that the two selections/PDB files contain #
    # the same number of atoms.                       #
    #                                                 #
    ###################################################

    if len(x2) != len(x1):
        print("\nError: \"" + first_obj_frame +\
              "\" and \"" + last_obj_frame +\
              "\" contain different number of residue/atoms.")
        print("\nPlease check your input and/or selections and try again.")
        return
    else:
        # Continue with representing modevectors!
        #########################################
        #                                       #
        # Delete old selection or object if it  #
        # exists so that it can be overwritten  #
        #                                       #
        #########################################
        save_view = cmd.get_view(output=1, quiet=1)
        cmd.delete(objectname)
        cmd.hide(representation="everything", selection=first_obj_frame)
        cmd.hide(representation="everything", selection=last_obj_frame)

    ###################################################
    #                                                 #
    # Begin drawing arrow tails                       #
    #                                                 #
    ###################################################

    arrowtail = []
    for mv in range(len(x1)):
        vectorx = x2[mv] - x1[mv]
        vectory = y2[mv] - y1[mv]
        vectorz = z2[mv] - z1[mv]
        length = sqrt(vectorx ** 2 + vectory ** 2 + vectorz ** 2)
        if length < cutoff:
            cutoff_counter += 1
            continue
        t = 1.0 - (cut / length)
        x2[mv] = x1[mv] + factor * t * vectorx
        y2[mv] = y1[mv] + factor * t * vectory
        z2[mv] = z1[mv] + factor * t * vectorz
        vectorx = x2[mv] - x1[mv]
        vectory = y2[mv] - y1[mv]
        vectorz = z2[mv] - z1[mv]
        length = sqrt(vectorx ** 2 + vectory ** 2 + vectorz ** 2)
        d = arrow_head_length  # Distance from arrow tip to arrow base
        t = 1.0 - (d / length)
        if notail:
            t = 0
        tail = [
            # Tail of cylinder
            CYLINDER, x1[mv], y1[mv], z1[mv]\
            , x1[mv] + (t + 0.01) * vectorx, y1[mv] + (t + 0.01) * vectory, z1[mv] + (t + 0.01) * vectorz\
            , arrow_tail_radius, tr, tg, tb, tr, tg, tb  # Radius and RGB for each cylinder tail
        ]
        if notail == 0:
            arrow.extend(tail)

        x = x1[mv] + t * vectorx
        y = y1[mv] + t * vectory
        z = z1[mv] + t * vectorz
        dx = x2[mv] - x
        dy = y2[mv] - y
        dz = z2[mv] - z
        seg = d / 100
        intfactor = int(factor)
        if version < 1.1:  # Version >= 1.1 has cone primitive
            for i in range(100, 0, -1):  # i=100 is tip of cone
                print(i)
                t1 = seg * i
                t2 = seg * (i + 1)
                radius = arrow_head_radius * (1.0 - i / (100.0))  # Radius of each disc that forms cone
                head = [
                    CYLINDER, x + t2 * dx, y + t2 * dy, z + t2 * dz\
                    , x + t1 * dx, y + t1 * dy, z + t1 * dz\
                    , radius, hr, hg, hb, hr, hg, hb  # Radius and RGB for slice of arrow head
                ]
                arrow.extend(head)
        else:
            head = [
                CONE, x, y, z, x + d * dx, y + d * dy, z + d * dz, arrow_head_radius, 0.0, hr, hg, hb, hr, hg, hb, 1.0, 1.0]
            arrow.extend(head)

##############################################################
#                                                            #
# Load the entire object into PyMOL                          #
#                                                            #
# Print statistics if requested by user                      #
#                                                            #
##############################################################

    if stat == "show":
        natoms = skipcounter + keepcounter
        print("\nTotal number of atoms = " + str(natoms))
        print("Atoms skipped = " + str(skipcounter))
        if keepcounter - cutoff_counter > 0:
            print("Atoms counted = " + str(keepcounter - cutoff_counter) + " (see PyMOL object \"" + objectname + "\")")
        else:
            print("Atoms counted = " + str(keepcounter - cutoff_counter) + " (Empty CGO object not loaded)")
        print("Atoms cutoff  = " + str(cutoff_counter))  # Note that cutoff occurs AFTER skipping!
    if keepcounter - cutoff_counter > 0:
        cmd.delete(objectname)
        cmd.load_cgo(arrow, objectname)  # Ray tracing an empty object will cause a segmentation fault.  No arrows = Do not display in PyMOL!!!
    cmd.show(representation="cartoon", selection=first_obj_frame)
    if (first_obj_frame != last_obj_frame):
        cmd.show(representation="cartoon", selection=last_obj_frame)
        cmd.hide(representation="cartoon", selection=last_obj_frame)
    cmd.bg_color(color="white")
    cmd.set_view(save_view)
    return

cmd.extend("modevectors", modevectors)

def gromacs_flag(name):
    try:
        devnull = open(os.devnull)
        subprocess.Popen([name], stdout=devnull, stderr=devnull).communicate()
    except OSError as e:
        if e.errno == errno.ENOENT:
            return False
    return True

def __init_plugin__(app=None):
    '''
    Add an entry to the PyMOL "Plugin" menu
    '''
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('G_Measure v.0.8', run_plugin_gui)


# global reference to avoid garbage collection of our dialog
dialog = None


def run_plugin_gui():
    '''
    Open our custom dialog
    '''
    global dialog

    if dialog is None:
        dialog = make_dialog()

    dialog.show()


def make_dialog():
    # entry point to PyMOL's API
    from pymol import cmd,stored
    from pymol.Qt import QtWidgets
    from pymol.Qt.utils import loadUi
    from pymol.Qt.utils import getSaveFileNameWithExt
    cmd.bg_color("white")
    cmd.remove('solvent')
    
    def showdialog(msgtitle,msgtxt):
        mb = QtWidgets.QMessageBox()
        mb.setIcon(QtWidgets.QMessageBox.Information)
        mb.setWindowTitle(msgtitle)
        mb.setText(msgtxt)
        mb.setStandardButtons(QtWidgets.QMessageBox.Ok)
        mb.exec_()
    try:
        import pandas as pd        
        import matplotlib
        import matplotlib.pyplot as plt
        from scipy.stats import gaussian_kde
        from mpl_toolkits.mplot3d import Axes3D        
        import mdtraj
        from sklearn.decomposition import PCA
        
    except:
    	showdialog('Note', 'Please install python library requeriments: pandas, matplotlib, scipy, mdtraj, sklearn')

    # create a new Window
    dialog = QtWidgets.QDialog()

    # populate the Window from our *.ui file which was created with the Qt Designer
    uifile = os.path.join(os.path.dirname(__file__), 'geo.ui')
    form = loadUi(uifile, dialog)
    form.progressBar.setProperty("value", 0)
    
    # callback for the "Ray" button
    stored.axe_frame = []
    stored.axe1_column_name = []
    stored.axe1_data= []
    stored.axe2_column_name = []
    stored.axe2_data= []
    def openCSVfileAxe1():
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(None,"Select a CSV file.", "","CSV (*.csv)", options=options)
        if fileName:
            form.le_axe1.setText(fileName)
            df = pd.read_csv(fileName)
            for itens in df.columns.tolist():
                stored.axe1_column_name.append(itens)
            for itens in df.iloc[:,-1]:
                stored.axe1_data.append(itens)
            for itens in df.iloc[:,-2]:
                stored.axe_frame.append(itens)
        else:
            return None

    def openCSVfileAxe2():
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getOpenFileName(None,"Select a CSV file.", "","CSV (*.csv)", options=options)
        if fileName:
            form.le_axe2.setText(fileName)
            df = pd.read_csv(fileName)
            for itens in df.columns.tolist():
                stored.axe2_column_name.append(itens)
            for itens in df.iloc[:,-1]:
                stored.axe2_data.append(itens)            
        else:
            return None

    def setupProcess():        
        # Run the process with a given command
        form.plain_status.setPlainText("Please wait.. Sham is running..")
        os.chdir(TEMP_PATH)
        if gromacs_flag('mdrun'):
            cmd ='g_sham -f g_sham2.xvg -ls free-energy-landscape.xpm'
        elif gromacs_flag('gmx'):
            cmd = 'gmx sham -f g_sham2.xvg -ls free-energy-landscape.xpm'
        os.system(cmd)

    def dataFEL():
        os.chdir(TEMP_PATH)
        xpm_file = 'free-energy-landscape.xpm'
        xpm_handle = open(xpm_file)
        xpm_data = []
        x_axis, y_axis = [], []
        letter_to_value = {}
        for line in xpm_handle:
            if line.startswith("/* x-axis"):
                x_ax = map(float, line.split()[2:-2]) # We trim the last value
                x_axis = list(x_ax)
            if line.startswith("/* y-axis"):
                y_ax = map(float, line.split()[2:-2]) # We trim the last value
                y_axis = list(y_ax)
            if line.startswith('"') and x_axis and y_axis: # Read data
                xpm_data.insert(0, line.strip().strip(',')[1:-1])
            if line.startswith('"') and len(line.split()) > 4:
                letter = line.split()[0][1:]
                value = float(line.split()[-2][1:-1])
                letter_to_value[letter] = value
        xpm_handle.close()
        txt_values = []
        data = []
        for y_index, data_value in enumerate(xpm_data):
            y_value = y_axis[y_index]
            for x_index, x_value in enumerate(x_axis):
                txt_values.append([x_value, y_value, letter_to_value[data_value[x_index]]])
            for x, y, z in txt_values:
                data.append ([float(x),float(y),float(z)])
        if len(stored.axe1_data) != 0 and len(stored.axe2_data) != 0:
            labels = [str(stored.axe1_column_name[-1]), str(stored.axe2_column_name[-1]), 'Gb_E (kj/mol)']
        else:
            labels = ['RMSD (nm)', 'RG (nm)', 'Gb_E (kj/mol)']
        df = pd.DataFrame.from_records(data, columns=labels)
        os.chdir(TEMP_PATH)
        df.to_csv('dataFEL.csv')   

    def clear():
        try:
            cmd.delete('not md')
        except:
            pass
        form.cb_res1.clear()
        form.cb_res2.clear()
        form.cb_res3.clear()
        form.cb_res4.clear()
        form.cb_chain.clear()
        stored.axe1_data.clear()
        stored.axe2_data.clear()
        form.le_axe1.clear()
        form.le_axe2.clear()
        form.bt_clear.setVisible(False)
        form.bt_getcsv.setVisible(False)
        form.bt_plot.setVisible(False)
        form.bt_run.setVisible(False)
        form.bt_unset.setVisible(False)
        form.bt_load.setVisible(True)
        form.cb_tool.setEnabled(True)
        form.sb_1pc.setVisible(False)
        form.sb_2pc.setVisible(False)
        form.label_x.setVisible(False)
        form.label_pca.setVisible(False)
        form.cb_selplot.setVisible(False)
        form.label_plot.setVisible(False)
        form.bt_browse_axe1.setEnabled(True)
        form.bt_browse_axe2.setEnabled(True)
        form.le_axe1.setEnabled(True)
        form.le_axe2.setEnabled(True)
        form.plain_status.setPlainText('Ready')
        stored.axe_frame = []
        stored.axe1_column_name = []
        stored.axe1_data= []
        stored.axe2_column_name = []
        stored.axe2_data= []
        try:
            cmd.color('green', 'md')
        except:
            pass
        try:
            shutil.rmtree(TRAJ_PATH)
        except:
            None

    def run():
        if form.cb_tool.currentText() == "Pincer angle":
            form.plain_status.setPlainText('Running Pincer Angle, please wait.')
            dataAngle = []
            dataFrame =[]
            count = 0
            form.progressBar.setMaximum(cmd.count_states('md'))
            for frame in range(cmd.count_states('md')):
                angle=cmd.angle(None, 'res1', 'res2', 'res3',state=frame)
                dataAngle.append(angle)
                dataFrame.append(frame)
                count += 1
                form.progressBar.setValue(count)
            os.chdir(TEMP_PATH)
            g_data =  pd.DataFrame()
            g_data['Frame'] = dataFrame
            g_data['Angle'] = dataAngle
            g_data.to_csv('dataAngle.csv')
            form.bt_run.setVisible(False)
            form.bt_unset.setVisible(False)
            form.bt_plot.setVisible(True)
            form.bt_getcsv.setVisible(True)
            form.bt_clear.setVisible(True)
            form.plain_status.setPlainText('Done.')
        elif form.cb_tool.currentText() == "Triangle area":
            form.plain_status.setPlainText('Running Triangle Area, please wait.')
            dataArea = []
            dataFrame =[]
            count = 0
            form.progressBar.setMaximum(cmd.count_states('md'))
            for frame in range(cmd.count_states('md')):
                distAB = cmd.distance(None,"res1","res2", state=frame)
                distAC = cmd.distance(None,"res1","res3", state=frame)
                distBC = cmd.distance(None,"res2","res3", state=frame)                        
                sPerimeter = (distAB + distAC + distBC) / 2                        
                area = (sPerimeter*(sPerimeter-distAB)*(sPerimeter-distAC)*(sPerimeter-distBC)) ** 0.5
                dataArea.append(area)
                dataFrame.append(frame)
                count += 1
                form.progressBar.setValue(count)
            os.chdir(TEMP_PATH)
            g_data =  pd.DataFrame()
            g_data['Frame'] = dataFrame
            g_data['Area'] = dataArea
            g_data.to_csv('dataArea.csv')
            form.bt_run.setVisible(False)
            form.bt_unset.setVisible(False)
            form.bt_plot.setVisible(True)
            form.bt_getcsv.setVisible(True)
            form.bt_clear.setVisible(True)
            form.plain_status.setPlainText('Done.')
        elif form.cb_tool.currentText() == "Dihedral angle":
            form.plain_status.setPlainText('Running Dihedral Angle, please wait.')
            dataDihedral = []
            dataFrame =[]
            count = 0
            form.progressBar.setMaximum(cmd.count_states('md'))
            for frame in range(cmd.count_states('md')):
                dihedral=cmd.get_dihedral('res1', 'res2', 'res3', 'res4', state=frame)
                dataDihedral.append(dihedral)
                dataFrame.append(frame)
                count += 1
                form.progressBar.setValue(count)
            os.chdir(TEMP_PATH)
            g_data =  pd.DataFrame()
            g_data['Frame'] = dataFrame
            g_data['Dihedral'] = dataDihedral
            g_data.to_csv('dataDihedral.csv')
            form.bt_run.setVisible(False)
            form.bt_unset.setVisible(False)
            form.bt_plot.setVisible(True)
            form.bt_getcsv.setVisible(True)
            form.bt_clear.setVisible(True)
            form.plain_status.setPlainText('Done.')
        elif form.cb_tool.currentText() == "RMSD":
            form.plain_status.setPlainText('Running RMSD, please wait.')
            dataRMSD = []
            dataFrame = []
            count = 0
            form.progressBar.setMaximum(cmd.count_states('md'))
            t = mdtraj.load(TRAJ_PATH)
            rmsd = mdtraj.rmsd(t, t, 1)
            for frame in range(cmd.count_states('md')):                
                dataRMSD.append(rmsd[frame])
                dataFrame.append(frame)
                count += 1
                form.progressBar.setValue(count)
            os.chdir(TEMP_PATH)
            g_data =  pd.DataFrame()
            g_data['Frame'] = dataFrame
            g_data['RMSD (nm)'] = dataRMSD
            g_data.to_csv('dataRMSD.csv')
            form.bt_run.setVisible(False)
            form.bt_unset.setVisible(False)
            form.bt_plot.setVisible(True)
            form.bt_getcsv.setVisible(True)
            form.bt_clear.setVisible(True)
            form.plain_status.setPlainText('Done.')
        elif form.cb_tool.currentText() == "RG":
            form.plain_status.setPlainText('Running RG, please wait.')
            dataRG = []
            dataFrame = []
            count = 0
            form.progressBar.setMaximum(cmd.count_states('md'))
            t = mdtraj.load(TRAJ_PATH)
            rg = mdtraj.compute_rg(t)
            for frame in range(cmd.count_states('md')):
                dataRG.append(rg[frame])
                dataFrame.append(frame)
                count += 1
                form.progressBar.setValue(count)
            os.chdir(TEMP_PATH)
            g_data =  pd.DataFrame()
            g_data['Frame'] = dataFrame
            g_data['RG (nm)'] = dataRG
            g_data.to_csv('dataRG.csv')
            form.bt_run.setVisible(False)
            form.bt_unset.setVisible(False)
            form.bt_plot.setVisible(True)
            form.bt_getcsv.setVisible(True)
            form.bt_clear.setVisible(True)
            form.plain_status.setPlainText('Done.')
        elif form.cb_tool.currentText() == "PDF":
            form.plain_status.setPlainText('Running PDF, please wait.')
            if len(stored.axe1_data) != 0 and len(stored.axe2_data) != 0:
                if len(stored.axe1_data) == len(stored.axe2_data):
                    form.progressBar.setMaximum(100)
                    form.progressBar.setValue(100)
                    os.chdir(TEMP_PATH)
                    g_data =  pd.DataFrame()
                    g_data['Frame'] = stored.axe_frame
                    g_data[str(stored.axe1_column_name[-1])] = stored.axe1_data
                    g_data[str(stored.axe2_column_name[-1])] = stored.axe2_data
                    g_data.to_csv('dataPDF.csv')
                    form.bt_run.setVisible(False)
                    form.bt_unset.setVisible(False)
                    form.bt_plot.setVisible(True)
                    form.bt_getcsv.setVisible(True)
                    form.bt_clear.setVisible(True)
                    form.plain_status.setPlainText('Done.')
                else:
                    showdialog('Note','The number of frames must be equal for both CSV files.')
            else:
                count = 0
                form.progressBar.setMaximum(cmd.count_states('md'))
                dataFrame =[]
                dataRMSD = []
                dataRG = []
                t = mdtraj.load(TRAJ_PATH)
                rg = mdtraj.compute_rg(t)
                rmsd = mdtraj.rmsd(t, t, 1)
                for frame in range(cmd.count_states('md')):
                    dataFrame.append(frame)
                    dataRMSD.append(rmsd[frame])
                    dataRG.append(rg[frame])
                    count += 1
                    form.progressBar.setValue(count)
                os.chdir(TEMP_PATH)
                g_data =  pd.DataFrame()
                g_data['Frame'] = dataFrame
                g_data['RMSD (nm)'] = dataRMSD
                g_data['RG (nm)'] = dataRG
                g_data.to_csv('dataPDF.csv')
                form.bt_run.setVisible(False)
                form.bt_unset.setVisible(False)
                form.bt_plot.setVisible(True)
                form.bt_getcsv.setVisible(True)
                form.bt_clear.setVisible(True)
                form.plain_status.setPlainText('Done.')
        elif form.cb_tool.currentText() == "PCA":
            form.plain_status.setPlainText('Running PCA, please wait.')
            dataPC1 =[]
            dataPC2 =[]
            dataPC3 =[]
            dataPC4 =[]
            dataPC5 =[]
            dataPC6 =[]
            dataPC7 =[]
            dataPC8 =[]
            dataPC9 =[]
            dataPC10 =[]
            dataPC1_exp =[]
            dataPC2_exp =[]
            dataPC3_exp =[]
            dataPC4_exp =[]
            dataPC5_exp =[]
            dataPC6_exp =[]
            dataPC7_exp =[]
            dataPC8_exp =[]
            dataPC9_exp =[]
            dataPC10_exp =[]
            dataFrame=[]
            traj = mdtraj.load(TRAJ_PATH)
            pca1 = PCA(n_components=10)
            traj.superpose(traj, 0)
            reduced_cartesian = pca1.fit_transform(traj.xyz.reshape(traj.n_frames, traj.n_atoms * 3))
            count = 0
            form.progressBar.setMaximum(cmd.count_states('md'))
            for frame in range(cmd.count_states('md')):
                dataPC1.append(reduced_cartesian[frame, 0])
                dataPC2.append(reduced_cartesian[frame, 1])
                dataPC3.append(reduced_cartesian[frame, 2])
                dataPC4.append(reduced_cartesian[frame, 3])
                dataPC5.append(reduced_cartesian[frame, 4])
                dataPC6.append(reduced_cartesian[frame, 5])
                dataPC7.append(reduced_cartesian[frame, 6])
                dataPC8.append(reduced_cartesian[frame, 7])
                dataPC9.append(reduced_cartesian[frame, 8])
                dataPC10.append(reduced_cartesian[frame, 9])
                dataFrame.append(frame)
                count += 1
                form.progressBar.setValue(count)
            g_data =  pd.DataFrame()
            g_data['Frame'] = dataFrame
            g_data['PC1'] = dataPC1
            g_data['PC2'] = dataPC2
            g_data['PC3'] = dataPC3
            g_data['PC4'] = dataPC4
            g_data['PC5'] = dataPC5
            g_data['PC6'] = dataPC6
            g_data['PC7'] = dataPC7
            g_data['PC8'] = dataPC8
            g_data['PC9'] = dataPC9
            g_data['PC10'] = dataPC10
            os.chdir(TEMP_PATH)
            g_data.to_csv('dataPCA.csv')
            form.bt_run.setVisible(False)
            form.bt_unset.setVisible(False)
            form.bt_plot.setVisible(True)
            form.bt_getcsv.setVisible(True)
            form.bt_clear.setVisible(True)
            form.sb_1pc.setVisible(True)
            form.sb_2pc.setVisible(True)
            form.label_x.setVisible(True)
            form.label_pca.setVisible(True)
            form.plain_status.setPlainText('Done.')
        elif form.cb_tool.currentText() == "PCA_EXP":
            form.plain_status.setPlainText('Running PCA_EXP, please wait.')
            pc = [0,1,2,3,4,5,6,7,8,9]
            PC = []
            stored.dataPC_exp = []            
            traj = mdtraj.load(TRAJ_PATH)
            pca1 = PCA(n_components=10)
            traj.superpose(traj, 0)
            reduced_cartesian = pca1.fit_transform(traj.xyz.reshape(traj.n_frames, traj.n_atoms * 3))
            count = 0
            form.progressBar.setMaximum(len(pc))
            for item in pc:
                stored.dataPC_exp.append(pca1.explained_variance_ratio_[item])
                PC.append(item + 1)
                count += 1
                form.progressBar.setValue(count)
            os.chdir(TEMP_PATH)
            g_data =  pd.DataFrame()
            g_data['PC'] = PC
            g_data['explained_variance_ratio'] = stored.dataPC_exp
            g_data.to_csv('dataPCE.csv')
            form.bt_run.setVisible(False)
            form.bt_unset.setVisible(False)
            form.bt_plot.setVisible(True)
            form.bt_getcsv.setVisible(True)
            form.bt_clear.setVisible(True)
            form.plain_status.setPlainText('Done.')
        elif form.cb_tool.currentText() == "Ramachandran map":
            form.plain_status.setPlainText('Running Ramachandran map, please wait.')
            frame = form.sb_frame_init.value()
            cmd.create('frame_'+str(frame),'md',source_state=frame,target_state=frame)
            cmd.frame(frame)
            angle = cmd.phi_psi('frame_'+str(frame))
            os.chdir(TEMP_PATH)
            g_data =  pd.DataFrame()
            phi = []
            psi = []
            count = 0
            form.progressBar.setMaximum(len(angle.items()))
            for md, phi_psi in angle.items():
                phi.append(phi_psi[0])
                psi.append(phi_psi[1])
                count += 1
                form.progressBar.setValue(count)
            g_data['phi'] = phi
            g_data['psi'] = psi
            g_data.to_csv('dataRAMA.csv')
            cmd.remove('frame_'+str(frame))
            form.bt_run.setVisible(False)
            form.bt_unset.setVisible(False)
            form.bt_plot.setVisible(True)
            form.bt_getcsv.setVisible(True)
            form.bt_clear.setVisible(True)
            form.plain_status.setPlainText('Done.')
        elif form.cb_tool.currentText() == "RMSF":
            form.plain_status.setPlainText('Running RMSF, please wait.')
            rmsfData = []
            res = []
            traj = mdtraj.load(TRAJ_PATH)
            traj.center_coordinates()
            rmsf = mdtraj.rmsf(traj, traj, 0, precentered=True)
            count = 0
            form.progressBar.setMaximum(len(stored.res))
            for item in rmsf:
                rmsfData.append(item)
                count += 1
                form.progressBar.setValue(count)
            for item in stored.res:
                res.append(item[0])
            os.chdir(TEMP_PATH)
            g_data =  pd.DataFrame()
            g_data['Residue'] = res
            g_data['RMSF'] = rmsfData
            g_data.to_csv('dataRMSF.csv')
            form.bt_run.setVisible(False)
            form.bt_unset.setVisible(False)
            form.bt_plot.setVisible(True)
            form.bt_getcsv.setVisible(True)
            form.bt_clear.setVisible(True)
            form.plain_status.setPlainText('Done.')
        elif form.cb_tool.currentText() == "DSSP":
            form.plain_status.setPlainText('Running DSSP, please wait.')
            stored.frameData = []
            dsspData=[]
            res = []
            DSSP_PATH = TEMP_PATH+"/dssp.pdb"
            cmd.remove('resn hoh')
            cmd.save(DSSP_PATH, 'name c+o+n+ca', state=0)
            traj = mdtraj.load(DSSP_PATH)
            dssp = mdtraj.compute_dssp(traj, simplified=True)
            count = 0
            form.progressBar.setMaximum(cmd.count_states('md'))
            for item in dssp:
                dsspData.append(item)
                count += 1
                form.progressBar.setValue(count)

            for frame in range(cmd.count_states('md')):
                stored.frameData.append(frame)

            for item in stored.res:
                res.append(item[0])
            os.chdir(TEMP_PATH)
            g_data =  pd.DataFrame()
            g_data['Residue Index'] = res
            for i in range(cmd.count_states('md')):
                g_data[str(i)] = dsspData[i]
                g_data.to_csv('dataDSSP.csv')
            form.bt_run.setVisible(False)
            form.bt_unset.setVisible(False)
            form.bt_plot.setVisible(False)
            form.bt_getcsv.setVisible(True)
            form.bt_clear.setVisible(True)
            form.plain_status.setPlainText('Done.')            
        elif form.cb_tool.currentText() == "Distance":
            form.plain_status.setPlainText('Running Distance, please wait.')
            distData = []
            frameData = []
            count = 0
            form.progressBar.setMaximum(cmd.count_states('md'))
            for frame in range(cmd.count_states('md')):
                dist = cmd.distance(None,"res1","res2", state=frame)
                distData.append(dist)
                frameData.append(frame)
                count += 1
                form.progressBar.setValue(count)
            os.chdir(TEMP_PATH)
            g_data =  pd.DataFrame()
            g_data['Frame'] = frameData
            g_data['Distance'] = distData
            g_data.to_csv('dataDist.csv')
            form.bt_run.setVisible(False)
            form.bt_unset.setVisible(False)
            form.bt_plot.setVisible(True)
            form.bt_getcsv.setVisible(True)
            form.bt_clear.setVisible(True)
            form.plain_status.setPlainText('Done.')
        elif form.cb_tool.currentText() == "Modevectors":
            form.plain_status.setPlainText('Running Modevectors, please wait.')
            count = 0
            frame_init = form.sb_frame_init.value()
            frame_final = form.sb_frame_final.value()
            form.progressBar.setMaximum(frame_final-frame_init)
            for frame in range(frame_final-frame_init):
                count += 1
                form.progressBar.setValue(count)            
            modevectors(first_obj_frame='md',last_obj_frame ='md', first_state=frame_init, last_state=frame_final)
            cmd.create('frame_init','md',source_state=frame_init,target_state=frame_init)
            cmd.create('frame_final','md',source_state=frame_final,target_state=frame_final)
            cmd.hide(selection='md')
            cmd.cartoon('tube', 'frame_init')
            cmd.cartoon('tube', 'frame_final')
            cmd.frame(frame_init)
            cmd.show('cgo', 'modevectors')
            cmd.color('marine', 'frame_init')
            cmd.color('purpleblue', 'frame_final')
            form.bt_run.setVisible(False)
            form.bt_unset.setVisible(False)
            form.bt_plot.setVisible(False)
            form.bt_getcsv.setVisible(False)
            form.bt_clear.setVisible(True)
            form.plain_status.setPlainText('Done.')
        elif form.cb_tool.currentText() == "Ligand Distance":
            form.plain_status.setPlainText('Running Ligand Distance, please wait.')
            distData = []
            frameData = []
            count = 0
            form.progressBar.setMaximum(cmd.count_states('md'))
            cmd.pseudoatom('cent_lig', selection='ligand', state=0)
            for frame in range(cmd.count_states('md')):
                dist = cmd.distance(None,"res1","cent_lig", state=frame)
                distData.append(dist)
                frameData.append(frame)
                count += 1
                form.progressBar.setValue(count)
            os.chdir(TEMP_PATH)
            g_data =  pd.DataFrame()
            g_data['Frame'] = frameData
            g_data['Distance'] = distData
            g_data.to_csv('dataLigDist.csv')
            form.bt_run.setVisible(False)
            form.bt_unset.setVisible(False)
            form.bt_plot.setVisible(True)
            form.bt_getcsv.setVisible(True)
            form.bt_clear.setVisible(True)
            form.plain_status.setPlainText('Done.')
        elif form.cb_tool.currentText() == "FEL":
            if gromacs_flag('mdrun') or gromacs_flag('gmx'):
                form.plain_status.setPlainText('Running FEL, please wait.')
                if len(stored.axe1_data) != 0 and len(stored.axe2_data) != 0:
                    if len(stored.axe1_data) == len(stored.axe2_data):
                        os.chdir(TEMP_PATH)
                        g_data =  pd.DataFrame()
                        g_data['a'] = stored.axe_frame
                        g_data['b'] = stored.axe1_data
                        g_data['c'] = stored.axe2_data
                        g_data.to_csv('g_sham2.xvg', sep='\t', index=False,  header=0)
                        with open('g_sham2.xvg', 'r') as fin:
                            data = fin.read().splitlines(True)
                        with open('g_sham2.xvg', 'w') as fout:
                            fout.writelines(data[1:])
                        setupProcess()
                        dataFEL()
                        form.bt_run.setVisible(False)
                        form.bt_unset.setVisible(False)
                        form.bt_plot.setVisible(True)
                        form.bt_getcsv.setVisible(True)
                        form.bt_clear.setVisible(True)
                        form.cb_selplot.setVisible(True)
                        form.label_plot.setVisible(True)
                        form.progressBar.setMaximum(100)
                        form.progressBar.setValue(100)
                        form.plain_status.setPlainText('Done.')
                    else:
                        showdialog('Note','The number of frames must be equal for both CSV files.')
                else:
                    count = 0
                    form.progressBar.setMaximum(cmd.count_states('md'))
                    dataFrame =[]
                    dataRMSD = []
                    dataRG = []
                    t = mdtraj.load(TRAJ_PATH)
                    rg = mdtraj.compute_rg(t)
                    rmsd = mdtraj.rmsd(t, t, 1)
                    for frame in range(cmd.count_states('md')):
                        dataFrame.append(frame)
                        dataRMSD.append(rmsd[frame])
                        dataRG.append(rg[frame])
                        count += 1
                        form.progressBar.setValue(count)                      
                    os.chdir(TEMP_PATH)
                    g_data =  pd.DataFrame()
                    g_data['a'] = dataFrame
                    g_data['b'] = dataRMSD
                    g_data['c'] = dataRG
                    g_data.to_csv('g_sham2.xvg', sep='\t', index=False,  header=0)
                    with open('g_sham2.xvg', 'r') as fin:
                        data = fin.read().splitlines(True)
                    with open('g_sham2.xvg', 'w') as fout:
                        fout.writelines(data[1:])
                    setupProcess()
                    dataFEL()
                    form.bt_run.setVisible(False)
                    form.bt_unset.setVisible(False)
                    form.bt_plot.setVisible(True)
                    form.bt_getcsv.setVisible(True)
                    form.bt_clear.setVisible(True)
                    form.cb_selplot.setVisible(True)
                    form.label_plot.setVisible(True)
                    form.plain_status.setPlainText('Done.')                  
            else:
                showdialog('Notice', 'GROMACS program must be intalled')
                clear()
    def seT():
        res_num1 = form.cb_res1.currentText()
        res_num2 = form.cb_res2.currentText()
        res_num3 = form.cb_res3.currentText()
        res_num4 = form.cb_res4.currentText()
        chain = form.cb_chain.currentText()
        num1 = res_num1.split('_')
        num2 = res_num2.split('_')
        num3 = res_num3.split('_')
        num4 = res_num4.split('_')
        lig = form.cb_lig.currentText()
        lig_num = lig.split('_')
        if form.cb_tool.currentText() == "Pincer angle":
            cmd.select('res1','/md//{0}/{1}/CA'.format(chain,num1[1]))
            cmd.select('res2','/md//{0}/{1}/CA'.format(chain,num2[1]))
            cmd.select('res3','/md//{0}/{1}/CA'.format(chain,num3[1]))
            cmd.show_as('cartoon')
            cmd.show_as('licorice', 'chain {0} and resi {1}+{2}+{3}'.format(chain,num1[1],num2[1],num3[1]))
            cmd.color('yellow', 'chain {0} and resi {1}+{2}+{3}'.format(chain,num1[1],num2[1],num3[1]))
            cmd.zoom('res2')
        elif form.cb_tool.currentText() == "Dihedral angle":
            cmd.select('res1','/md//{0}/{1}/CA'.format(chain,num1[1]))
            cmd.select('res2','/md//{0}/{1}/CA'.format(chain,num2[1]))
            cmd.select('res3','/md//{0}/{1}/CA'.format(chain,num3[1]))
            cmd.select('res4','/md//{0}/{1}/CA'.format(chain,num4[1]))
            cmd.show_as('cartoon')
            cmd.show_as('licorice', 'chain {0} and resi {1}+{2}+{3}+{4}'.format(chain,num1[1],num2[1],num3[1],num4[1]))
            cmd.color('yellow', 'chain {0} and resi {1}+{2}+{3}+{4}'.format(chain,num1[1],num2[1],num3[1],num4[1]))
            cmd.zoom('res2')
        elif form.cb_tool.currentText() == "Triangle area":
            cmd.select('res1','/md//{0}/{1}/CA'.format(chain,num1[1]))
            cmd.select('res2','/md//{0}/{1}/CA'.format(chain,num2[1]))
            cmd.select('res3','/md//{0}/{1}/CA'.format(chain,num3[1]))
            cmd.show_as('cartoon')
            cmd.show_as('licorice', 'chain {0} and resi {1}+{2}+{3}'.format(chain,num1[1],num2[1],num3[1]))
            cmd.color('yellow', 'chain {0} and resi {1}+{2}+{3}'.format(chain,num1[1],num2[1],num3[1]))
            cmd.zoom('res2')
        elif form.cb_tool.currentText() == "PDF":
            cmd.color('yellow', 'md')
            cmd.save(TRAJ_PATH, 'n. CA', state=0)
        elif form.cb_tool.currentText() == "RMSD":
            cmd.color('yellow', 'md')
            cmd.save(TRAJ_PATH, 'n. CA', state=0)
        elif form.cb_tool.currentText() == "RG":
            cmd.color('yellow', 'md')
            cmd.save(TRAJ_PATH, 'n. CA', state=0)
        elif form.cb_tool.currentText() == "FEL":
            cmd.color('yellow', 'md')
            cmd.save(TRAJ_PATH, 'n. CA', state=0)
        elif form.cb_tool.currentText() == "PCA":
            cmd.color('yellow', 'md')
            cmd.save(TRAJ_PATH, 'n. CA', state=0)
        elif form.cb_tool.currentText() == "PCA_EXP":
            cmd.color('yellow', 'md')
            cmd.save(TRAJ_PATH, 'n. CA', state=0)
        elif form.cb_tool.currentText() == "DSSP":
            cmd.color('yellow', 'md')
            cmd.save(TRAJ_PATH, 'n. CA', state=0)
        elif form.cb_tool.currentText() =="Ramachandran map":
            cmd.color('yellow', 'md')
        elif form.cb_tool.currentText() =="RMSF":
            cmd.color('yellow', 'md')
            cmd.save(TRAJ_PATH, 'n. CA', state=0)        
        elif form.cb_tool.currentText() =="Distance":
            cmd.select('res1','/md//{0}/{1}/CA'.format(chain,num1[1]))
            cmd.select('res2','/md//{0}/{1}/CA'.format(chain,num2[1]))
            cmd.show_as('cartoon')
            cmd.show_as('licorice', 'res1')
            cmd.color('yellow', 'res1 and res2')
            cmd.zoom('res2')
        elif form.cb_tool.currentText() =="Ligand Distance":            
            cmd.select('ligand','resi {}'.format(lig_num[1]))
            cmd.select('res1','/md//{0}/{1}/CA'.format(chain,num1[1]))
            cmd.show_as('cartoon')
            cmd.show_as('licorice', 'ligand and res1')
            cmd.color('yellow', 'ligand and res1')
            cmd.zoom('ligand')
        form.bt_run.setVisible(True)
        form.bt_unset.setVisible(True)
        form.bt_set.setVisible(False)
        form.cb_res1.setEnabled(False)
        form.cb_res2.setEnabled(False)
        form.cb_res3.setEnabled(False)
        form.cb_res4.setEnabled(False)
        form.cb_chain.setEnabled(False)
        form.sb_frame_init.setEnabled(False)
        form.sb_frame_final.setEnabled(False)
        form.cb_lig.setEnabled(False)
        form.plain_status.setPlainText('Read. Click on Run.')

          
    
    def res_hide():
        
        if form.cb_tool.currentText() == "PDF":
            form.label_frame1.setVisible(False)
            form.label_frame2.setVisible(False)
            form.sb_frame_init.setVisible(False)
            form.sb_frame_final.setVisible(False)
            form.label_res1.setVisible(False)
            form.cb_res1.setVisible(False)
            form.label_res2.setVisible(False)
            form.cb_res2.setVisible(False)
            form.label_res3.setVisible(False)
            form.cb_res3.setVisible(False)
            form.label_res4.setVisible(False)
            form.cb_res4.setVisible(False)
            form.cb_lig.setVisible(False)
            form.label_lig.setVisible(False)
            form.label_axe1.setVisible(True)
            form.label_axe2.setVisible(True)
            form.le_axe1.setVisible(True)
            form.le_axe2.setVisible(True)
            form.bt_browse_axe1.setVisible(True)
            form.bt_browse_axe2.setVisible(True)
            form.plain_status.setPlainText('The Functional Density Function is calculated using different values ​​of mainchain dihedral angles from the considered residue, the mainchain conformation of the equivalent residue between frames. Click on Load.')

        elif form.cb_tool.currentText() == "RG":
            form.label_frame1.setVisible(False)
            form.label_frame2.setVisible(False)
            form.sb_frame_init.setVisible(False)
            form.sb_frame_final.setVisible(False)
            form.label_res1.setVisible(False)
            form.cb_res1.setVisible(False)
            form.label_res2.setVisible(False)
            form.cb_res2.setVisible(False)
            form.label_res3.setVisible(False)
            form.cb_res3.setVisible(False)
            form.label_res4.setVisible(False)
            form.cb_res4.setVisible(False)
            form.cb_lig.setVisible(False)
            form.label_lig.setVisible(False)
            form.label_axe1.setVisible(False)
            form.label_axe2.setVisible(False)
            form.le_axe1.setVisible(False)
            form.le_axe2.setVisible(False)
            form.bt_browse_axe1.setVisible(False)
            form.bt_browse_axe2.setVisible(False)
            form.plain_status.setPlainText('Computes the radius of gyration of the protein in a function of frames. Click on Load.')

        elif form.cb_tool.currentText() == "DSSP":
            form.label_frame1.setVisible(False)
            form.label_frame2.setVisible(False)
            form.sb_frame_init.setVisible(False)
            form.sb_frame_final.setVisible(False)
            form.label_res1.setVisible(False)
            form.cb_res1.setVisible(False)
            form.label_res2.setVisible(False)
            form.cb_res2.setVisible(False)
            form.label_res3.setVisible(False)
            form.cb_res3.setVisible(False)
            form.label_res4.setVisible(False)
            form.cb_res4.setVisible(False)
            form.cb_lig.setVisible(False)
            form.label_lig.setVisible(False)
            form.label_axe1.setVisible(False)
            form.label_axe2.setVisible(False)
            form.le_axe1.setVisible(False)
            form.le_axe2.setVisible(False)
            form.bt_browse_axe1.setVisible(False)
            form.bt_browse_axe2.setVisible(False)
            form.plain_status.setPlainText('Compute Dictionary of protein secondary structure (DSSP) secondary structure in a function of frames. Click on Load.')

        elif form.cb_tool.currentText() == "Ramachandran map":
            form.label_frame1.setVisible(True)
            form.label_frame2.setVisible(False)
            form.sb_frame_init.setVisible(True)
            form.sb_frame_final.setVisible(False)
            form.label_res1.setVisible(False)
            form.cb_res1.setVisible(False)
            form.label_res2.setVisible(False)
            form.cb_res2.setVisible(False)
            form.label_res3.setVisible(False)
            form.cb_res3.setVisible(False)
            form.label_res4.setVisible(False)
            form.cb_res4.setVisible(False)
            form.cb_lig.setVisible(False)
            form.label_lig.setVisible(False)
            form.label_axe1.setVisible(False)
            form.label_axe2.setVisible(False)
            form.le_axe1.setVisible(False)
            form.le_axe2.setVisible(False)
            form.bt_browse_axe1.setVisible(False)
            form.bt_browse_axe2.setVisible(False)
            form.plain_status.setPlainText('The Ramachandran plot is the 2d plot of the ϕ-ψ torsion angles of the protein backbone. It provides a simple view of the conformation of a selected frame of protein . Click on Load.')

        elif form. cb_tool.currentText() == "RMSD":
            form.label_frame1.setVisible(False)
            form.label_frame2.setVisible(False)
            form.sb_frame_init.setVisible(False)
            form.sb_frame_final.setVisible(False)
            form.label_res1.setVisible(False)
            form.cb_res1.setVisible(False)
            form.label_res2.setVisible(False)
            form.cb_res2.setVisible(False)
            form.label_res3.setVisible(False)
            form.cb_res3.setVisible(False)
            form.label_res4.setVisible(False)
            form.cb_res4.setVisible(False)
            form.cb_lig.setVisible(False)
            form.label_lig.setVisible(False)
            form.label_axe1.setVisible(False)
            form.label_axe2.setVisible(False)
            form.le_axe1.setVisible(False)
            form.le_axe2.setVisible(False)
            form.bt_browse_axe1.setVisible(False)
            form.bt_browse_axe2.setVisible(False)
            form.plain_status.setPlainText('Compares two protein structures frames (first frame the first frame with the consecutive frames) by computing the root mean square deviation (RMSD). Click on Load.')

        elif form. cb_tool.currentText() == "RMSF":
            form.label_frame1.setVisible(False)
            form.label_frame2.setVisible(False)
            form.sb_frame_init.setVisible(False)
            form.sb_frame_final.setVisible(False)
            form.label_res1.setVisible(False)
            form.cb_res1.setVisible(False)
            form.label_res2.setVisible(False)
            form.cb_res2.setVisible(False)
            form.label_res3.setVisible(False)
            form.cb_res3.setVisible(False)
            form.label_res4.setVisible(False)
            form.cb_res4.setVisible(False)
            form.cb_lig.setVisible(False)
            form.label_lig.setVisible(False)
            form.label_axe1.setVisible(False)
            form.label_axe2.setVisible(False)
            form.le_axe1.setVisible(False)
            form.le_axe2.setVisible(False)
            form.bt_browse_axe1.setVisible(False)
            form.bt_browse_axe2.setVisible(False)
            form.plain_status.setPlainText('RMSF stands for root mean square fluctuation. This is a numerical measurement similar to RMSD, but instead of indicating positional differences between entire structures over time, RMSF is a calculation of individual residue flexibility, or how much a particular residue moves (fluctuates) during a simulation. Click on Load.')

        elif form.cb_tool.currentText() == "Pincer angle":
            form.label_frame1.setVisible(False)
            form.label_frame2.setVisible(False)
            form.sb_frame_init.setVisible(False)
            form.sb_frame_final.setVisible(False)
            form.label_res1.setVisible(True)
            form.cb_res1.setVisible(True)
            form.label_res2.setVisible(True)
            form.cb_res2.setVisible(True)
            form.label_res3.setVisible(True)
            form.cb_res3.setVisible(True)
            form.label_res4.setVisible(False)
            form.cb_res4.setVisible(False)
            form.cb_lig.setVisible(False)
            form.label_lig.setVisible(False)
            form.label_axe1.setVisible(False)
            form.label_axe2.setVisible(False)
            form.le_axe1.setVisible(False)
            form.le_axe2.setVisible(False)
            form.bt_browse_axe1.setVisible(False)
            form.bt_browse_axe2.setVisible(False)
            form.plain_status.setPlainText('Measure the residue carbon alpha pincer angle. Click on Load.')

        elif form.cb_tool.currentText() == "Triangle area":
            form.label_frame1.setVisible(False)
            form.label_frame2.setVisible(False)
            form.sb_frame_init.setVisible(False)
            form.sb_frame_final.setVisible(False)
            form.label_res1.setVisible(True)
            form.cb_res1.setVisible(True)
            form.label_res2.setVisible(True)
            form.cb_res2.setVisible(True)
            form.label_res3.setVisible(True)
            form.cb_res3.setVisible(True)
            form.label_res4.setVisible(False)
            form.cb_res4.setVisible(False)
            form.cb_lig.setVisible(False)
            form.label_lig.setVisible(False)
            form.label_axe1.setVisible(False)
            form.label_axe2.setVisible(False)
            form.le_axe1.setVisible(False)
            form.le_axe2.setVisible(False)
            form.bt_browse_axe1.setVisible(False)
            form.bt_browse_axe2.setVisible(False)
            form.plain_status.setPlainText('Measures the area of ​​the triangle formed by the 3 selected alpha carbons. Click on Load.')

        elif form.cb_tool.currentText() == "Dihedral angle":
            form.label_frame1.setVisible(False)
            form.label_frame2.setVisible(False)
            form.sb_frame_init.setVisible(False)
            form.sb_frame_final.setVisible(False)
            form.label_res1.setVisible(True)
            form.cb_res1.setVisible(True)
            form.label_res2.setVisible(True)
            form.cb_res2.setVisible(True)
            form.label_res3.setVisible(True)
            form.cb_res3.setVisible(True)
            form.label_res4.setVisible(True)
            form.cb_res4.setVisible(True)
            form.cb_lig.setVisible(False)
            form.label_lig.setVisible(False)
            form.label_axe1.setVisible(False)
            form.label_axe2.setVisible(False)
            form.le_axe1.setVisible(False)
            form.le_axe2.setVisible(False)
            form.bt_browse_axe1.setVisible(False)
            form.bt_browse_axe2.setVisible(False)
            form.plain_status.setPlainText('Measure the residue carbon alpha dihedral angles. Click on Load.')

        elif form.cb_tool.currentText() == "Distance":
            form.label_frame1.setVisible(False)
            form.label_frame2.setVisible(False)
            form.sb_frame_init.setVisible(False)
            form.sb_frame_final.setVisible(False)
            form.label_res1.setVisible(True)
            form.cb_res1.setVisible(True)
            form.label_res2.setVisible(True)
            form.cb_res2.setVisible(True)
            form.label_res3.setVisible(False)
            form.cb_res3.setVisible(False)
            form.label_res4.setVisible(False)
            form.cb_res4.setVisible(False)
            form.cb_lig.setVisible(False)
            form.label_lig.setVisible(False)
            form.label_axe1.setVisible(False)
            form.label_axe2.setVisible(False)
            form.le_axe1.setVisible(False)
            form.le_axe2.setVisible(False)
            form.bt_browse_axe1.setVisible(False)
            form.bt_browse_axe2.setVisible(False)
            form.plain_status.setPlainText('Measure the residue carbon alpha distance. Click on Load.')

        elif form.cb_tool.currentText() == "Ligand Distance":
            form.label_frame1.setVisible(False)
            form.label_frame2.setVisible(False)
            form.sb_frame_init.setVisible(False)
            form.sb_frame_final.setVisible(False)
            form.label_res1.setVisible(True)
            form.cb_res1.setVisible(True)
            form.label_res2.setVisible(False)
            form.cb_res2.setVisible(False)
            form.label_res3.setVisible(False)
            form.cb_res3.setVisible(False)
            form.label_res4.setVisible(False)
            form.cb_res4.setVisible(False)
            form.cb_lig.setVisible(True)
            form.label_lig.setVisible(True)
            form.label_axe1.setVisible(False)
            form.label_axe2.setVisible(False)
            form.le_axe1.setVisible(False)
            form.le_axe2.setVisible(False)
            form.bt_browse_axe1.setVisible(False)
            form.bt_browse_axe2.setVisible(False)
            form.plain_status.setPlainText('Measure the residue carbon alpha and ligand center of mass distance. Click on Load.')

        elif form.cb_tool.currentText() == "Modevectors":
            form.label_frame1.setVisible(True)
            form.label_frame2.setVisible(True)
            form.sb_frame_init.setVisible(True)
            form.sb_frame_final.setVisible(True)
            form.label_res1.setVisible(False)
            form.cb_res1.setVisible(False)
            form.label_res2.setVisible(False)
            form.cb_res2.setVisible(False)
            form.label_res3.setVisible(False)
            form.cb_res3.setVisible(False)
            form.label_res4.setVisible(False)
            form.cb_res4.setVisible(False)
            form.cb_lig.setVisible(False)
            form.label_lig.setVisible(False)
            form.label_axe1.setVisible(False)
            form.label_axe2.setVisible(False)
            form.le_axe1.setVisible(False)
            form.le_axe2.setVisible(False)
            form.bt_browse_axe1.setVisible(False)
            form.bt_browse_axe2.setVisible(False)
            form.plain_status.setPlainText('Visualize the direction of motion between two specified frames. Click on Load.')

        elif form.cb_tool.currentText() == "FEL":
            form.label_frame1.setVisible(False)
            form.label_frame2.setVisible(False)
            form.sb_frame_init.setVisible(False)
            form.sb_frame_final.setVisible(False)
            form.label_res1.setVisible(False)
            form.cb_res1.setVisible(False)
            form.label_res2.setVisible(False)
            form.cb_res2.setVisible(False)
            form.label_res3.setVisible(False)
            form.cb_res3.setVisible(False)
            form.label_res4.setVisible(False)
            form.cb_res4.setVisible(False)
            form.cb_lig.setVisible(False)
            form.label_lig.setVisible(False)
            form.label_axe1.setVisible(True)
            form.label_axe2.setVisible(True)
            form.le_axe1.setVisible(True)
            form.le_axe2.setVisible(True)
            form.bt_browse_axe1.setVisible(True)
            form.bt_browse_axe2.setVisible(True)
            form.plain_status.setPlainText('FEL represents a mapping of all possible conformations a molecule adopted during a simulation, together with their corresponding energy reported as the Gibbs Free Energy. FEL are represented using two variables that reflect specific properties of the system and measure conformational variability. RG measure the torsion angle around a specific bond or the radius of gyration of the protein, and the the RMSD measure the deviation with respective native state (First frame). Click on Load.')

        elif form.cb_tool.currentText() == "PCA":
            form.label_frame1.setVisible(False)
            form.label_frame2.setVisible(False)
            form.sb_frame_init.setVisible(False)
            form.sb_frame_final.setVisible(False)
            form.label_res1.setVisible(False)
            form.cb_res1.setVisible(False)
            form.label_res2.setVisible(False)
            form.cb_res2.setVisible(False)
            form.label_res3.setVisible(False)
            form.cb_res3.setVisible(False)
            form.label_res4.setVisible(False)
            form.cb_res4.setVisible(False)
            form.cb_lig.setVisible(False)
            form.label_lig.setVisible(False)
            form.label_axe1.setVisible(False)
            form.label_axe2.setVisible(False)
            form.le_axe1.setVisible(False)
            form.le_axe2.setVisible(False)
            form.bt_browse_axe1.setVisible(False)
            form.bt_browse_axe2.setVisible(False)
            form.plain_status.setPlainText('Create ten component PCA model, and project our data down into this reduced dimensional space. Click on Load.')
        
        elif form.cb_tool.currentText() == "PCA_EXP":
            form.label_frame1.setVisible(False)
            form.label_frame2.setVisible(False)
            form.sb_frame_init.setVisible(False)
            form.sb_frame_final.setVisible(False)
            form.label_res1.setVisible(False)
            form.cb_res1.setVisible(False)
            form.label_res2.setVisible(False)
            form.cb_res2.setVisible(False)
            form.label_res3.setVisible(False)
            form.cb_res3.setVisible(False)
            form.label_res4.setVisible(False)
            form.cb_res4.setVisible(False)
            form.cb_lig.setVisible(False)
            form.label_lig.setVisible(False)
            form.label_axe1.setVisible(False)
            form.label_axe2.setVisible(False)
            form.le_axe1.setVisible(False)
            form.le_axe2.setVisible(False)
            form.bt_browse_axe1.setVisible(False)
            form.bt_browse_axe2.setVisible(False)
            form.plain_status.setPlainText('Gives the PCA variance explained solely by the i+1st dimension. Click on Load.')

        else:
            form.label_frame1.setVisible(False)
            form.label_frame2.setVisible(False)
            form.sb_frame_init.setVisible(False)
            form.sb_frame_final.setVisible(False)
            form.label_res1.setVisible(True)
            form.cb_res1.setVisible(True)
            form.label_res2.setVisible(True)
            form.cb_res2.setVisible(True)
            form.label_res3.setVisible(True)
            form.cb_res3.setVisible(True)
            form.label_res4.setVisible(False)
            form.cb_res4.setVisible(False)
            form.cb_lig.setVisible(False)
            form.label_lig.setVisible(False)
            form.label_axe1.setVisible(False)
            form.label_axe2.setVisible(False)
            form.le_axe1.setVisible(False)
            form.le_axe2.setVisible(False)
            form.bt_browse_axe1.setVisible(False)
            form.bt_browse_axe2.setVisible(False)
            form.plain_status.setPlainText('Ready. Select a tool.')    
    
    
    def remove_rep(List):
        l = []
        for i in List:
            if i not in l:
                l.append(i)
        l.sort()
        return l
    def unset():
        res_num1 = form.cb_res1.currentText()
        res_num2 = form.cb_res2.currentText()
        res_num3 = form.cb_res3.currentText()
        res_num4 = form.cb_res4.currentText()
        chain = form.cb_chain.currentText()
        num1 = res_num1.split('_')
        num2 = res_num2.split('_')
        num3 = res_num3.split('_')
        num4 = res_num4.split('_')
        lig = form.cb_lig.currentText()
        lig_num = lig.split('_')
        if form.cb_tool.currentText() == "Pincer angle":
            cmd.delete('res1 and res2 and res3')
            cmd.color('green', 'md')
        elif form.cb_tool.currentText() == "Dihedral angle":
            cmd.delete('res1 and res2 and res3 and res4')
            cmd.color('green', 'md')
        elif form.cb_tool.currentText() == "Triangle area":
            cmd.delete('res1 and res2 and res3')
            cmd.color('green', 'md')
        elif form.cb_tool.currentText() == "PDF":
            cmd.color('green', 'md')
            try:
                shutil.rmtree(TRAJ_PATH)
            except:
                None
        elif form.cb_tool.currentText() == "RMSD":
            cmd.color('green', 'md')
            try:
                shutil.rmtree(TRAJ_PATH)
            except:
                None
        elif form.cb_tool.currentText() == "RG":
            cmd.color('green', 'md')
            try:
                shutil.rmtree(TRAJ_PATH)
            except:
                None
        elif form.cb_tool.currentText() == "FEL":
            cmd.color('green', 'md')
            try:
                shutil.rmtree(TRAJ_PATH)
            except:
                None
        elif form.cb_tool.currentText() == "PCA":
            cmd.color('green', 'md')
            try:
                shutil.rmtree(TRAJ_PATH)
            except:
                None
        elif form.cb_tool.currentText() == "PCA_EXP":
            cmd.color('green', 'md')
            try:
                shutil.rmtree(TRAJ_PATH)
            except:
                None
        elif form.cb_tool.currentText() == "DSSP":
            cmd.color('green', 'md')
            try:
                shutil.rmtree(TRAJ_PATH)
            except:
                None
        elif form.cb_tool.currentText() =="Ramachandran map":
            cmd.color('green', 'md')
        elif form.cb_tool.currentText() =="RMSF":
            cmd.color('green', 'md')
            try:
                shutil.rmtree(TRAJ_PATH)
            except:
                None        
        elif form.cb_tool.currentText() =="Distance":
            cmd.delete('res1 and res2')
            cmd.color('green', 'md')
        elif form.cb_tool.currentText() =="Ligand Distance":            
            cmd.delete('res1 and ligand')
            cmd.color('green', 'md')
        form.bt_run.setVisible(False)
        form.bt_unset.setVisible(False)
        form.bt_set.setVisible(True)
        form.cb_res1.setEnabled(True)
        form.cb_res2.setEnabled(True)
        form.cb_res3.setEnabled(True)
        form.cb_res4.setEnabled(True)
        form.cb_chain.setEnabled(True)
        form.sb_frame_init.setEnabled(True)
        form.sb_frame_final.setEnabled(True)
        form.cb_lig.setEnabled(True)

    def load():
        # get form data
        if len(stored.axe1_data) != 0 and len(stored.axe2_data) != 0:
            if len(stored.axe1_data) != 0 and len(stored.axe2_data) != 0:
                form.bt_load.setVisible(False)
                form.bt_run.setVisible(True)
                form.bt_browse_axe1.setEnabled(False)
                form.bt_browse_axe2.setEnabled(False)
                form.le_axe1.setEnabled(False)
                form.le_axe2.setEnabled(False)
                form.plain_status.setPlainText('Using CSV files. Click on Run')
            else:
                showdialog('Note','The number of frames must be equal for both CSV files.')
        else:
            form.cb_res1.setEnabled(True)
            form.cb_res2.setEnabled(True)
            form.cb_res3.setEnabled(True)
            form.cb_res4.setEnabled(True)
            form.cb_chain.setEnabled(True)
            form.bt_browse_axe1.setEnabled(False)
            form.bt_browse_axe2.setEnabled(False)
            form.le_axe1.setEnabled(False)
            form.le_axe2.setEnabled(False)
            try:
                print('Change current name '+str(cmd.get_object_list(selection='(all)'))+' by "md"')
                obj_name = cmd.get_object_list(selection='(all)')
                cmd.set_name(obj_name[0], 'md')
                stored.res=[]
                cmd.iterate("(name ca)","stored.res.append((resi,resn))")       
                for item in stored.res:
                    form.cb_res1.addItem(item[1]+'_'+item[0])
                    form.cb_res2.addItem(item[1]+'_'+item[0])
                    form.cb_res3.addItem(item[1]+'_'+item[0])
                    form.cb_res4.addItem(item[1]+'_'+item[0])

                stored.lig=[]
                cmd.iterate("(organic)","stored.lig.append((resi,resn))")       
                for item in remove_rep(stored.lig):
                    form.cb_lig.addItem(item[1]+'_'+item[0])
                stored.ch=[]    
                for ch in cmd.get_chains('md'):
                	stored.ch.append(ch)
                for item in remove_rep(stored.ch):
                    form.cb_chain.addItem(item)
                form.sb_frame_init.setMinimum(1)
                form.sb_frame_init.setMaximum(cmd.count_states('md')-1)
                form.sb_frame_init.setValue(1)
                form.sb_frame_final.setMinimum(1)
                form.sb_frame_final.setMaximum(cmd.count_states('md')-1)
                form.sb_frame_final.setValue(cmd.count_states('md'))
                form.bt_set.setVisible(True)
                form.bt_load.setVisible(False)
                form.cb_tool.setEnabled(False)
                form.sb_frame_init.setEnabled(True)
                form.sb_frame_final.setEnabled(True)
                
            except:
                showdialog('Note',"Invalid trajectory loaded")


    
    def save_csv():
        options = QtWidgets.QFileDialog.Options()
        options |= QtWidgets.QFileDialog.DontUseNativeDialog
        fileName, _ = QtWidgets.QFileDialog.getSaveFileName(None,"Save CSV file",""," Save CSV Files (*.csv)", options=options)
        if fileName:
            if form.cb_tool.currentText() == "Pincer angle":
                data = TEMP_PATH+'/dataAngle.csv'
            elif form.cb_tool.currentText() == "Dihedral angle":
                data = TEMP_PATH+'/dataDihedral.csv'
            elif form.cb_tool.currentText() == "Triangle area":
                data = TEMP_PATH+'/dataArea.csv'
            elif form.cb_tool.currentText() == "PDF":
                data = TEMP_PATH+'/dataPDF.csv'
            elif form.cb_tool.currentText() == "RMSD":
                data = TEMP_PATH+'/dataRMSD.csv'
            elif form.cb_tool.currentText() == "RG":
                data = TEMP_PATH+'/dataRG.csv'
            elif form.cb_tool.currentText() == "FEL":
                data = TEMP_PATH+'/dataFEL.csv',
            elif form.cb_tool.currentText() == "PCA":
                data = TEMP_PATH+'/dataPCA.csv'
            elif form.cb_tool.currentText() == "PCA_EXP":
                data = TEMP_PATH+'/dataPCE.csv'
            elif form.cb_tool.currentText() =="Ramachandran map":
                data = TEMP_PATH+'/dataRAMA.csv'
            elif form.cb_tool.currentText() =="DSSP":
                data = TEMP_PATH+'/dataDSSP.csv'
            elif form.cb_tool.currentText() =="RMSF":
                data = TEMP_PATH+'/dataRMSF.csv'        
            elif form.cb_tool.currentText() =="Distance":
                data = TEMP_PATH+'/dataDist.csv'
            elif form.cb_tool.currentText() =="Ligand Distance":
                data = TEMP_PATH+'/dataLigDist.csv'
            df = pd.read_csv(data)
            if fileName.endswith('.csv'):
                df.to_csv(fileName)
            else:
                df.to_csv(fileName+'.csv')
        else:
            pass

    def plottingData():
        if form.cb_tool.currentText() == "Pincer angle":
            option='Angle'
            data = TEMP_PATH+'/dataAngle.csv'
        elif form.cb_tool.currentText() == "Dihedral angle":
            option='Dihedral'
            data = TEMP_PATH+'/dataDihedral.csv'
        elif form.cb_tool.currentText() == "Triangle area":
            option='Area'
            data = TEMP_PATH+'/dataArea.csv'
        elif form.cb_tool.currentText() == "PDF":
            option='PDF'
            data = TEMP_PATH+'/dataPDF.csv'
        elif form.cb_tool.currentText() == "RMSD":
            option='RMSD (nm)'
            data = TEMP_PATH+'/dataRMSD.csv'
        elif form.cb_tool.currentText() == "RG":
            option = 'RG (nm)'
            data = TEMP_PATH+'/dataRG.csv'
        elif form.cb_tool.currentText() == "FEL":
            option='FEL'
            data = TEMP_PATH+'/dataFEL.csv',
        elif form.cb_tool.currentText() == "PCA":
            option='PCA'
            data = TEMP_PATH+'/dataPCA.csv'
        elif form.cb_tool.currentText() == "PCA_EXP":
            option='PCA_EXP'
            data = TEMP_PATH+'/dataPCE.csv'
        elif form.cb_tool.currentText() =="Ramachandran map":
            option="Ramachandran map"
            data = TEMP_PATH+'/dataRAMA.csv'
        elif form.cb_tool.currentText() =="RMSF":
            option="RMSF"
            data = TEMP_PATH+'/dataRMSF.csv'        
        elif form.cb_tool.currentText() =="Distance":
            option="Distance"
            data = TEMP_PATH+'/dataDist.csv'

        elif form.cb_tool.currentText() =="Ligand Distance":
            option="Ligand Distance"
            data = TEMP_PATH+'/dataLigDist.csv'
        if option == 'PDF':
            df = pd.read_csv(data)

            fig, (ax1) = plt.subplots(nrows=1)            
            # Setting data
            c_name = df.columns.tolist()
            x = df[str(c_name[-2])]
            y = df[str(c_name[-1])]

            # Calculate the point density
            xy = np.vstack([x,y])
            z = gaussian_kde(xy)(xy)

            # Sort the points by density, so that the densest points are plotted last
            idx = z.argsort()
            x, y, z = x[idx], y[idx], z[idx]

            # Setting plot type 
            pdf = ax1.scatter(x, y, c = z, s = 50, edgecolor = '', cmap=plt.cm.jet)

            # Plot title
            ax1.set_title(str(c_name[-1])+' by '+str(c_name[-2]))

            # Hide right and top spines
            ax1.spines['right'].set_visible(False)
            ax1.spines['top'].set_visible(False)
            ax1.yaxis.set_ticks_position('left')
            ax1.xaxis.set_ticks_position('bottom')

            # Set x and y limits
            xmin = x.min() - 1
            xmax = x.max() + 1
            ymin = y.min() - 1
            ymax = y.max() + 1        
            plt.xlim(xmin, xmax)
            plt.ylim(ymin, ymax)

            # Set x and y labels
            plt.xlabel(str(c_name[-1]))
            plt.ylabel(str(c_name[-2]))

            # Adding the color bar 
            colbar = plt.colorbar(pdf)
            colbar.set_label('Probability Density Function')     
            plt.show()

        elif option == 'RG (nm)' or option == 'RMSD (nm)' or option == 'Angle' or option == 'Dihedral' or option == 'Area':
            df = pd.read_csv(data)
            fig, (ax1) = plt.subplots(nrows=1)
            ax1.plot(df['Frame'], df[option])
            ax1.set_title(option + ' by Time')
            ax1.spines['right'].set_visible(False)
            ax1.spines['top'].set_visible(False)
            ax1.yaxis.set_ticks_position('left')
            ax1.xaxis.set_ticks_position('bottom')
            plt.xlabel('Frame')
            xmin1 = df['Frame'].min() - 1
            xmax1 = df['Frame'].max() + 1
            plt.xlim(xmin1, xmax1)
            plt.ylabel(option)        
            plt.show()
        
        elif option == 'FEL':
            if form.cb_selplot.currentText() == "3D":
                df = pd.read_csv(TEMP_PATH+'/dataFEL.csv')
                fig = plt.figure(figsize=(15,10))
                fig.suptitle('Free Energy Landscape', fontsize=20)
                ax = fig.gca(projection='3d')
                c_name = df.columns.tolist()
                ax.set_xlabel(str(c_name[-3]), fontsize=15)
                ax.set_ylabel(str(c_name[-2]), fontsize=15)
                ax.set_zlabel('Gibbs Free Energy (kj/mol)', fontsize=15)
                ax = fig.gca(projection='3d')
                ax.plot_trisurf(df[str(c_name[-3])], df[str(c_name[-2])], df['Gb_E (kj/mol)'], cmap=plt.cm.jet, linewidth=0, antialiased=False)
                    
                # to Add a color bar which maps values to colors.
                surf=ax.plot_trisurf(df[str(c_name[-3])], df[str(c_name[-2])], df['Gb_E (kj/mol)'], cmap=plt.cm.jet, linewidth=0, antialiased=False)
                colbar = fig.colorbar( surf, shrink=0.5, aspect=5)
                colbar.set_label('Gibbs Free Energy (kj/mol)')
                ax.tricontourf(df[str(c_name[-3])], df[str(c_name[-2])], df['Gb_E (kj/mol)'], zdir='z', offset=-1, cmap=plt.cm.jet)
                    
                # Rotate it
                ax.view_init(30, 15)
                plt.show()
            else:
                df = pd.read_csv(TEMP_PATH+'/dataFEL.csv')
                c_name = df.columns.tolist()
                z = df['Gb_E (kj/mol)']
                X = df[str(c_name[-1])]
                Y = df[str(c_name[-2])]
                fig, ax = plt.subplots()
                fig.suptitle('Free Energy Landscape', fontsize=20)
                trico = ax.tricontourf(df[str(c_name[-3])], df[str(c_name[-2])], df['Gb_E (kj/mol)'], zdir='z', offset=-1, cmap=plt.cm.jet)
                ax.set_xlabel(str(c_name[-3]), fontsize=15)
                ax.set_ylabel(str(c_name[-2]), fontsize=15)
                colbar = fig.colorbar(trico, shrink=0.5, aspect=5)
                colbar.set_label('Gibbs Free Energy (kj/mol)')
                plt.show()
        elif option == 'PCA':
            df = pd.read_csv(data)
            z = df['Frame']
            X = df['PC{}'.format(form.sb_1pc.value())]
            Y = df['PC{}'.format(form.sb_2pc.value())]
            plt.figure()
            plt.scatter(X, Y, c=z)
            plt.xlabel('PC{}'.format(form.sb_1pc.value()))
            plt.ylabel('PC{}'.format(form.sb_2pc.value()))
            plt.title('Cartesian coordinate PCA')
            cbar = plt.colorbar()
            cbar.set_label('Frame')
            plt.show()
        elif option == 'PCA_EXP':
            df = pd.read_csv(data)
            plt.plot(np.cumsum(stored.dataPC_exp))
            plt.xlabel('number of components')
            plt.ylabel('cumulative explained variance')
            plt.show()
        elif option == "Ramachandran map":
            df = pd.read_csv(data)
            plt.title("Ramachandran map")
            plt.scatter(df['phi'], df['psi'])
            plt.xlim([-180, 180])
            plt.ylim([-180, 180])
            plt.plot([-180, 180], [0, 0], color="black")
            plt.plot([0, 0], [-180, 180], color="black")
            plt.locator_params(axis='x', nbins=7)
            plt.xlabel(r'$\phi$')
            plt.ylabel(r'$\psi$')
            plt.grid()
            plt.show()

        elif option == 'RMSF':
            df = pd.read_csv(data)
            x,y = df['Residue'], df['RMSF']
            fig, (ax1) = plt.subplots(nrows=1)
            ax1.plot(x,y)
            ax1.set_title('RMSF by residues')
            ax1.spines['right'].set_visible(False)
            ax1.spines['top'].set_visible(False)
            ax1.set_xticks(x[::100])
            ax1.set_xticklabels(x[::100], rotation='vertical')
            plt.xlabel('Residue Index')
            plt.ylabel('RMSF')        
            plt.show()
        elif option == 'Distance' or option == 'Ligand Distance':
            df = pd.read_csv(data)
            x,y = df['Frame'], df['Distance']
            fig, (ax1) = plt.subplots(nrows=1)
            ax1.plot(x,y)
            ax1.set_title('Distance')
            ax1.spines['right'].set_visible(False)
            ax1.spines['top'].set_visible(False)
            ax1.set_xticks(x[::100])
            ax1.set_xticklabels(x[::100], rotation='vertical')
            plt.xlabel('Frame')
            plt.ylabel('Distance')        
            plt.show()
        

    form.cb_tool.currentIndexChanged.connect(res_hide)
    form.bt_load.clicked.connect(load)
    form.bt_run.clicked.connect(run)
    form.bt_set.clicked.connect(seT)
    form.bt_unset.clicked.connect(unset)
    form.bt_plot.clicked.connect(plottingData)
    form.bt_getcsv.clicked.connect(save_csv)
    form.bt_clear.clicked.connect(clear)
    form.bt_browse_axe1.clicked.connect(openCSVfileAxe1)
    form.bt_browse_axe2.clicked.connect(openCSVfileAxe2)
    form.le_axe1.setReadOnly(True)
    form.le_axe2.setReadOnly(True)
    form.cb_tool.addItem("Pincer angle")
    form.cb_tool.addItem("Dihedral angle")
    form.cb_tool.addItem("Triangle area")
    form.cb_tool.addItem("PDF")
    form.cb_tool.addItem("RMSD")
    form.cb_tool.addItem("RG")
    form.cb_tool.addItem("FEL")
    form.cb_tool.addItem("PCA")
    form.cb_tool.addItem("PCA_EXP")
    form.cb_tool.addItem("Ramachandran map")
    form.cb_tool.addItem("RMSF")
    form.cb_tool.addItem("DSSP")
    form.cb_tool.addItem("Distance")
    form.cb_tool.addItem("Modevectors")
    form.cb_tool.addItem("Ligand Distance")
    form.cb_selplot.addItem('2D')
    form.cb_selplot.addItem('3D')
    form.sb_1pc.setValue(1)
    form.sb_1pc.setMaximum(10)
    form.sb_1pc.setMinimum(1)
    form.sb_2pc.setValue(2)
    form.sb_2pc.setMaximum(10)
    form.sb_2pc.setMinimum(1)
    form.label_plot.setVisible(False)
    form.cb_selplot.setVisible(False)
    form.label_pca.setVisible(False)
    form.sb_1pc.setVisible(False)
    form.label_x.setVisible(False)
    form.sb_2pc.setVisible(False)
    form.bt_unset.setVisible(False)
    form.bt_clear.setVisible(False)
    form.bt_set.setVisible(False)
    form.bt_run.setVisible(False)
    form.bt_plot.setVisible(False)
    form.bt_getcsv.setVisible(False)
    form.sb_frame_init.setEnabled(False)
    form.sb_frame_final.setEnabled(False)
    form.setWindowTitle("G_Measures {}".format(_version_))
    form.plain_status.setReadOnly(True)
    return dialog
