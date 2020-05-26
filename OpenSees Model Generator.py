from IPython import get_ipython
def __reset__(): get_ipython().magic('reset -sf')

# import OpenSeesPy rendering module
from openseespy.postprocessing.Get_Rendering import *

import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

ops.wipe()

#read building summary spreadsheet
bldg = pd.read_excel (r'bldg_table_removed small area.xlsx')
bldg_arr= bldg.to_numpy()

#read original building spreadpsheet
original_bldg= pd.read_excel (r'bldg_removed no wall.xlsx')

#read original wall spreadsheet
original_wall= pd.read_excel (r'wall revised.xlsx')

#record 5 eigenvalues for each building
num_eigen= 5

#initialize matrix to store period data
period_data= np.zeros((bldg.shape[0],num_eigen))

#loop through all buildings to create linear elastic 3D stick models
for bldg_num in range(bldg.shape[0]):
    #read bldg ID
    bldg_id= bldg.iat[bldg_num,20]
    
    #read processed wall data spreadsheet
    wall= pd.read_excel (r'Bldg Plan Jan 4/Wall Table/wall_table%d.xlsx' %(bldg_id))
    wall.fillna(0, inplace=True)
    wall_arr = wall.to_numpy()
    
    # =============================================================================
    #   Set modelbuilder
    # =============================================================================
    ops.model('basic', '-ndm', 3)
    
    # =============================================================================
    #   Define Material
    # =============================================================================
    #1.material ID 2.elastic modulus (float) 3.Poisson’s ratio (float) 4.mass density (float) (optional)
    conc_mat = ["ElasticIsotropic",0,0,0,0]
   
    conc_mat[1] = 1         #material ID 
    conc_mat[3] = 0.17      #Poisson's ratio
    conc_mat[4] = 2400      #mass density
    #elastic modulus
    original_bldg_row= original_bldg.loc[original_bldg['BuildingID'] == bldg_id]

    if pd.isnull(original_bldg_row.iat[0,47]):
        conc_mat[2] = 4.5*10**9*math.sqrt(4000*0.00689476)
        #conc_mat[2] = (3300*(4000*0.00689476)**0.5+6900)*(conc_mat[4]/2300)**1.5*10**6     #take comp strength as 4000 is no data available
    else:
        conc_mat[2] = 4.5*10**9*math.sqrt(original_bldg_row.iat[0,47]*0.00689476)
        #conc_mat[2] = (3300*(original_bldg_row.iat[0,47]*0.00689476)**0.5+6900)*(conc_mat[4]/2300)**1.5*10**6
        
    # =============================================================================
    #     Define Node
    # =============================================================================
    story_height = bldg.iat[bldg_num,4]
    
    num_story= bldg.iat[bldg_num,3]
    
    #joint nodes
    node_data = np.zeros((wall.shape[0]*(num_story+1),4))
    
    floor = 0
    
    #1.node ID 2.x coord 3.y coord 4.z coord
    for story in range(num_story+1):
        node_data[floor:floor+wall.shape[0],0] = list(range(1,wall.shape[0]+1))+np.ones(wall.shape[0])*(story+1)*1000
        node_data[floor:floor+wall.shape[0],1:3] = wall_arr[:,1:3]
        node_data[floor:floor+wall.shape[0],3] = np.ones(wall.shape[0])*story*story_height
        floor = floor + wall.shape[0]
    
    for i in range(node_data.shape[0]):
        ops.node(int(node_data[i,0]), node_data[i,1], node_data[i,2], node_data[i,3])
    
#    #master node at center of rigidity for rigid diaphragm
#    cor_node_data= np.zeros((num_story,4))
#    
#    cor_node_data[:,0]= [x for x in list(range(1,num_story+1))]                #nodes start at 1 on floor 2
#    cor_node_data[:,1]= bldg_arr[bldg_num,17]                                  #COR x coord
#    cor_node_data[:,2]= bldg_arr[bldg_num,18]                                  #COR y coord
#    cor_node_data[:,3]= [x*story_height for x in list(range(1,num_story+1))]   #height of each story
#    
#    for i in range(cor_node_data.shape[0]):
#        ops.node(int(cor_node_data[i,0]), cor_node_data[i,1], cor_node_data[i,2], cor_node_data[i,3])
    
    #floor plate center of mass calculation
    #1.area 2.mass 3.x*mass 4.y*mass
    #all standard metric
    mass_calc= np.zeros((wall.shape[0],4))
        
    mass_calc[:,0]= np.multiply(wall_arr[:,3],wall_arr[:,4])
    mass_calc[:,1]= conc_mat[4]*story_height*mass_calc[:,0]
    mass_calc[:,2]= np.multiply(wall_arr[:,1], mass_calc[:,1])
    mass_calc[:,3]= np.multiply(wall_arr[:,2], mass_calc[:,1])
    
    overall_x_com= bldg.iat[bldg_num,11]
    overall_y_com= bldg.iat[bldg_num,12]
    
    seismic_mass= bldg.iat[bldg_num,2]*(bldg.iat[bldg_num,8]*conc_mat[4]+1200/9.81);
    #back calculate slab COM given available data
    slab_x_com= (overall_x_com*(seismic_mass+sum(mass_calc[:,1]))- sum(mass_calc[:,2]))/seismic_mass
    slab_y_com= (overall_y_com*(seismic_mass+sum(mass_calc[:,1]))- sum(mass_calc[:,3]))/seismic_mass   

    #mass nodes
    mass_node_data= np.zeros((num_story,10))
    
    #1.node ID 2.x coord 3.y coord 4.z coord 5 to 10:mass in 6 DOFs
    mass_node_data[:,0]= [x*10000 for x in list(range(1,num_story+1))]
    mass_node_data[:,1]= slab_x_com
    mass_node_data[:,2]= slab_y_com
    mass_node_data[:,3]= [x*story_height for x in list(range(1,num_story+1))]
    mass_node_data[:,4:6]= seismic_mass
    mass_node_data[:,9]= (bldg.iat[bldg_num,0]**2+bldg.iat[bldg_num,1]**2)*seismic_mass/12

    for i in range(mass_node_data.shape[0]):
        ops.node(int(mass_node_data[i,0]), mass_node_data[i,1], mass_node_data[i,2], mass_node_data[i,3], '-mass', mass_node_data[i,4], mass_node_data[i,5], mass_node_data[i,6], mass_node_data[i,7], mass_node_data[i,8], mass_node_data[i,9])

    #add mass to node to represent DL on columns
#    dl_data= np.zeros(((num_story)*wall.shape[0],7))
#    
#    wall_data= original_wall.loc[original_wall['BuildingID'] == bldg_id]
#    wall_data= wall_data.loc[wall_data['WallFloor'] > 1]
#    wall_data_arr = wall_data.to_numpy()
#    
#    floor = 0
#    
#    for i in range(num_story):
#        dl_data[floor:floor+wall.shape[0],0]= node_data[floor+wall.shape[0]:floor+wall.shape[0]*2,0]
#        dl_data[floor:floor+wall.shape[0],1]= np.multiply(wall_data_arr[:,18],wall_data_arr[:,19])*original_bldg.iat[bldg_num,38]**2*0.0254**2*(1200/9.81+bldg.iat[bldg_num,8]*conc_mat[4])
#        dl_data[floor:floor+wall.shape[0],2]= dl_data[floor:floor+wall.shape[0],1]
#        dl_data[floor:floor+wall.shape[0],3]= dl_data[floor:floor+wall.shape[0],1]
#        dl_data[floor:floor+wall.shape[0],4]= dl_data[floor:floor+wall.shape[0],1]
#        dl_data[floor:floor+wall.shape[0],5]= dl_data[floor:floor+wall.shape[0],1]
#        dl_data[floor:floor+wall.shape[0],6]= dl_data[floor:floor+wall.shape[0],1]
#        floor = floor + wall.shape[0]
#    
#    for i in range(dl_data.shape[0]):
#        ops.mass(dl_data[i,0], dl_data[i,1], dl_data[i,2], dl_data[i,3], dl_data[i,4], dl_data[i,5], dl_data[i,6])
#    
    # =============================================================================
    #     Set Boundary Condition
    # =============================================================================
    #1.node 2 to 6: 1=fixed, 0=free
    constraint_data= np.ones((wall.shape[0],7))
    
    constraint_data[:,0]= node_data[0:wall.shape[0],0]
    
    for i in range(wall.shape[0]):
        ops.fix(constraint_data[i,0], int(constraint_data[i,1]), int(constraint_data[i,2]), int(constraint_data[i,3]), int(constraint_data[i,4]), int(constraint_data[i,5]), int(constraint_data[i,6]))
    
    #constraints for mass nodes
    mass_constraint_data= np.ones((num_story,7))
    
    mass_constraint_data[:,0]= mass_node_data[0:mass_node_data.shape[0],0]
    
    mass_constraint_data[:,1:3]= 0
    mass_constraint_data[:,6]= 1
    
    for i in range(mass_constraint_data.shape[0]):
        ops.fix(mass_constraint_data[i,0], 0,0,1,1,1,0)
    
    # =============================================================================
    #     Apply Wall Constraints
    # =============================================================================
    
    original_wall_data= original_wall.loc[original_wall['BuildingID'] == bldg_id]
    typical_floor= max(original_wall_data['WallFloor'])
    original_wall_data= original_wall_data.loc[original_wall_data['WallFloor'] == typical_floor]
    
    attached= np.zeros((original_wall_data.shape[0],2))
    coupled= np.zeros((original_wall_data.shape[0],2))
    
    attached_row= 0
    coupled_row= 0

    for i in range(original_wall_data.shape[0]):
        if np.isnan(original_wall_data.iat[i,9])!= 1:
            attached[attached_row,0]= original_wall_data.iat[i,2]
            attached[attached_row,1]= original_wall_data.iat[i,9]
            attached_row= attached_row+ 1
        if np.isnan(original_wall_data.iat[i,10])!= 1: 
            coupled[coupled_row,0]= original_wall_data.iat[i,2]
            coupled[coupled_row,1]= original_wall_data.iat[i,10]
            coupled_row= coupled_row+ 1
    attached= np.delete(attached,np.where(~attached.any(axis=1))[0], axis=0)
    coupled= np.delete(coupled,np.where(~coupled.any(axis=1))[0], axis=0)
    
    #assume coupled walls are attached
    combined= np.concatenate((attached, coupled))
    
    wall_constraint_organizer= np.zeros((attached.shape[0]+coupled.shape[0],original_wall_data.shape[0]+10))
            
    for combined_row in range(combined.shape[0]):
        if combined[combined_row,1] in wall_constraint_organizer[:,:]:
            row = np.where(wall_constraint_organizer == combined[combined_row,1])[0][0]
            column= np.where(wall_constraint_organizer[row,:] == 0)[0][0]
            wall_constraint_organizer[row,column]= combined[combined_row,0]
        else:
            row = np.where(wall_constraint_organizer[:,0] == 0)[0][0]
            wall_constraint_organizer[row,0:2]= combined[combined_row,:]
    
    wall_constraint_organizer= np.delete(wall_constraint_organizer,np.where(~wall_constraint_organizer.any(axis=1))[0], axis=0)
    
    wall_constraint_organizer[wall_constraint_organizer==0] = np.nan
            
    for current_row in range(wall_constraint_organizer.shape[0]):
        if wall_constraint_organizer[current_row,0]!='nan':
            for row in range(wall_constraint_organizer.shape[0]):
                if wall_constraint_organizer[row,0]!='nan':
                    if set(wall_constraint_organizer[current_row,:]).isdisjoint(set(wall_constraint_organizer[row,:])):
                        pass
                    else:
                        temp= np.union1d(wall_constraint_organizer[current_row,:], wall_constraint_organizer[row,:])
                        temp= temp[~np.isnan(temp)]
                        wall_constraint_organizer[current_row,0:temp.shape[0]]= temp
                        if row!=current_row:
                            wall_constraint_organizer[row,:]= 'nan'
                            
    wall_constraint_organizer[np.isnan(wall_constraint_organizer)] = 0    
    
    original_wall_data_arr= original_wall_data.to_numpy()
    
    wall_constraint_data= wall_constraint_organizer
    
    for row in range(wall_constraint_organizer.shape[0]):
        if wall_constraint_organizer[row,0]!= 0:
            end_column= np.where(wall_constraint_organizer[row,:] == 0)[0][0]
            for column in range(end_column):
                wall_constraint_data[row,column]= np.where(original_wall_data_arr[:,2] == wall_constraint_organizer[row,column])[0][0]+1

    
    check= np.zeros((1000,2))
    check_row= 0
    
    for row in range(wall_constraint_data.shape[0]):
        if wall_constraint_data[row,0]!= 0:
            end_column= np.where(wall_constraint_data[row,:] == 0)[0][0]
            for column in range(1,end_column):
                for floor in range(num_story):
                    ops.equalDOF(int(wall_constraint_data[row,0]+1000*(floor+2)), int(wall_constraint_data[row,column]+1000*(floor+2)),3,4,5)
                    check[check_row,0]= wall_constraint_data[row,0]+1000*(floor+2)
                    check[check_row,1]= wall_constraint_data[row,column]+1000*(floor+2)
                    check_row= check_row+1
                        
    # =============================================================================
    #     Define Coordinate Transformation
    # =============================================================================
    ops.geomTransf('Linear', 1, 1, 0, 0)
    
    # =============================================================================
    #     Define Elements
    # =============================================================================
    
    #in local coordinate system: 0.J 1.Iz 2.Iy
    moment_of_inertia= np.zeros((wall.shape[0],3))
    
    for i in range(wall.shape[0]):
        moment_of_inertia[i,0]= wall.iat[i,3]*wall.iat[i,4]**3/3
        if wall.iat[i,7]=='E/W':
            moment_of_inertia[i,1]= wall.iat[i,3]**3*wall.iat[i,4]/12
            moment_of_inertia[i,2]= wall.iat[i,4]**3*wall.iat[i,3]/12
        elif wall.iat[i,7]=='N/S':
            moment_of_inertia[i,1]= wall.iat[i,4]**3*wall.iat[i,3]/12
            moment_of_inertia[i,2]= wall.iat[i,3]**3*wall.iat[i,4]/12 
    
    element_data= np.zeros((wall.shape[0]*num_story,11))

    row= 0
        
    for story in range(num_story):
        element_data[row:row+ wall.shape[0],0]= list(range(1,wall.shape[0]+1))+np.ones(wall.shape[0])*(story+1)*1000    #element ID
        #node 1 at the bottom and node 2 at the top of each floor
        element_data[row:row+ wall.shape[0],1]= node_data[row:row+wall.shape[0],0]                      #node 1 (i)
        element_data[row:row+ wall.shape[0],2]= node_data[row+wall.shape[0]:row+wall.shape[0]*2,0]      #node 2 (j)
        element_data[row:row+ wall.shape[0],3]= np.multiply(wall_arr[:,3],wall_arr[:,4])                #area
        element_data[row:row+ wall.shape[0],4]= conc_mat[2]                                             #elastic modulus
        element_data[row:row+ wall.shape[0],5]= element_data[row:row+ wall.shape[0],4]/(2*(1+0.17))     #shear modulus
        element_data[row:row+ wall.shape[0],6]= moment_of_inertia[:,0]                                  #J
        element_data[row:row+ wall.shape[0],7]= moment_of_inertia[:,1]*0.7                              #70% of Iy
        element_data[row:row+ wall.shape[0],8]= moment_of_inertia[:,2]*0.7                              #70% of Iz
        element_data[row:row+ wall.shape[0],9]= 1                                                       #coordinate transfer flag
        element_data[row:row+ wall.shape[0],10]= conc_mat[4]*element_data[row:row+ wall.shape[0],3]     #mass/length
    
        row= row+ wall.shape[0]
        
    for i in range(element_data.shape[0]):
        ops.element('elasticBeamColumn',int(element_data[i,0]),int(element_data[i,1]),int(element_data[i,2]),element_data[i,3],element_data[i,4],element_data[i,5],element_data[i,6],element_data[i,7],element_data[i,8],int(element_data[i,9]),'-mass',element_data[i,10])
    
#     =============================================================================
#         Define Rigid Diaphragm
#     =============================================================================
    for i in range(num_story):
        ops.rigidDiaphragm(3, (i+1)*10000, *list(node_data[(i+1)*wall.shape[0]:(i+1)*wall.shape[0]+ wall.shape[0],0]))    

    # =============================================================================
    #     Define Gravity Load in Columns
    # =============================================================================
    
#    #define time series
#    ops.timeSeries('Constant', 1)
#    
#    #define pattern
#    ops.pattern('Plain', 1, 1)
#    
#    #define nodal loads
#    wall_data= original_wall.loc[original_wall['BuildingID'] == bldg_id]
#   #note: need to modify line below to ensure the largest floor is taken    
#   wall_data= wall_data.loc[wall_data['WallFloor'] > 1]
#    wall_data_arr = wall_data.to_numpy()
#    
#    gravity_load_data= np.zeros((element_data.shape[0],7))
#    
#    row= 0
#    
#    gravity_load_data[:,0]= node_data[wall.shape[0]:node_data.shape[0],0]
#    
#    for story in range(num_story):
#        gravity_load_data[row:row+ wall.shape[0],3]= -1*np.multiply(wall_data_arr[:,18],wall_data_arr[:,19])*original_bldg.iat[bldg_num,38]**2*0.0254**2*(1200+bldg.iat[bldg_num,8]*conc_mat[4]*9.81)
#        
#        row= row+ wall.shape[0]
#    for i in range(gravity_load_data.shape[0]):
#        ops.load(int(gravity_load_data[i,0]), gravity_load_data[i,1], gravity_load_data[i,2], gravity_load_data[i,3], gravity_load_data[i,4], gravity_load_data[i,5], gravity_load_data[i,6])

    # =============================================================================
    #       Analysis Generation
    # =============================================================================
    #eigen command
    eigen_values= ops.eigen(num_eigen)
    
    period= 2*math.pi/np.sqrt(eigen_values)
    
    period_data[bldg_num,:]= period
    
    # render the model after defining all the nodes and elements
    #plot_model()
            
    ops.wipe()

#plot period vs. story
plt.plot(period_data[:,0], bldg_arr[:,3], 'o', color='black',markersize=5)
plt.ylim(0, 35)
plt.xlim(0, 6)
plt.xlabel('Fundamental period (s)')
plt.ylabel('Storeys above grade')
plt.plot([0,3.5], [0,35])

