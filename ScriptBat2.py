# -*- coding: utf-8 -*-
"""
Created on Fri Dec 14 12:15:13 2018

ForPower

Description: To estimate sediment trapping and accumulation in the reservoir

Scripted by Mohit Kaura, originally written in MatLab, 
translated into Python 3.6 by Josh Benjamin
"""
from openpyxl import load_workbook
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import copy

####Sedimentation##############################################################
#Import variables from Damdata Excel Sheet

XX = "E" #Signifies the E column - Battambang 2 - change this letter to change the dam

wb = load_workbook(filename = 'Damdata.xlsx', data_only=True) #Damdata Excel Sheet required
damdata = wb['DamData']
hvsv_ws = wb['hvsv']

#load cell values into memory
initial_res_vol = damdata["{}2".format(XX)].value #Initial Reservoir Volume (m3)
high_supp_elev = damdata["{}3".format(XX)].value #High Supply Level (m)
high_supply_vol = damdata["{}4".format(XX)].value #High Supply Volume (m3)
low_supply_elev = damdata["{}5".format(XX)].value #Low Supply Level (m)
low_supply_vol = damdata["{}6".format(XX)].value #Low Supply Volume (m3)
design_discharge = damdata["{}7".format(XX)].value #Design Discharge (m3/s)
mean_ann_energy = damdata["{}8".format(XX)].value #Mean Annual Energy (GWh)
low_def_rate= damdata["{}9".format(XX)].value #Low deforestation rate
avg_def_rate= damdata["{}10".format(XX)].value #Avg Deforestation rate
high_def_rate= damdata["{}11".format(XX)].value #High deforestation rate
rated_head= damdata["{}12".format(XX)].value #Rated Head (m)
sed_bulk_den = damdata["{}13".format(XX)].value #Bulk sediment density
alpha_trap = damdata["{}14".format(XX)].value #Alpha (TE)

a = damdata["{}15".format(XX)].value #constant a
b = damdata["{}16".format(XX)].value #constant b
c = damdata["{}17".format(XX)].value #constant c
d = damdata["{}18".format(XX)].value #constant d

g = damdata["{}19".format(XX)].value #Accln due to gravity (m/s2)

forcov1 = damdata["{}20".format(XX)].value #Forest cover: Scenario1 (percentage)
forcov2 = damdata["{}21".format(XX)].value #Forest cover: Scenario2 (percentage)
forcov3 = damdata["{}22".format(XX)].value #Forest cover: Scenario3 (percentage)
forcov4 = damdata["{}23".format(XX)].value #Forest cover: Scenario4 (percentage)
forcov2015 = damdata["{}24".format(XX)].value #Forest cover: 2015 (percentage)

sed1 = damdata["{}25".format(XX)].value # initial sediment Scenario1
sed2 = damdata["{}26".format(XX)].value # initial sediment Scenario2
sed3 = damdata["{}27".format(XX)].value # initial sediment Scenario3
sed4 = damdata["{}28".format(XX)].value # initial sediment Scenario4

den_water = damdata["{}29".format(XX)].value #Density of water (kg/m3)
mu = damdata["{}30".format(XX)].value #dam efficiency
cap_factor = damdata["{}31".format(XX)].value #capacity factor
elec = damdata["{}32".format(XX)].value #Electricty cost per Kwh ($)
disc = damdata["{}33".format(XX)].value #Discount Rate
impacting_area = damdata["{}34".format(XX)].value #Annuity
res_height = damdata["{}35".format(XX)].value #Reservoir Height

hvsv = pd.read_excel("Damdata.xlsx",sheet_name="hvsv", usecols=[6,7], skiprows=[0], nrows=14)
hvsv = np.flipud(hvsv.values) #convert from pandas dataframe to array, flip upside down

water_height = hvsv[:,0] #Water Elevation
res_volume = hvsv[:,1]
res_volume = copy.deepcopy(res_volume) #needed to separate from hvsv object
res_volume1 = copy.deepcopy(res_volume) #Reservoir volume at that Water Elevation
res_volume2 = copy.deepcopy(res_volume)
res_volume3 = copy.deepcopy(res_volume)
res_volume4 = copy.deepcopy(res_volume)


####Scenarios##################################################################
"""
SCENARIO 1: Controlled Deforestation
SCENARIO 2: Unmanaged Deforestation
SCENARIO 3: Excessive Deforestation
SCENARIO 4: Conservation and proper management
Every term ending with 1 relates to SCENARIO 1, 2 for SCENARIO 2 and so forth
"""
year = np.arange(1,122)

n = year.size
m = 1
# n is the length of year

# forcover = 100 shows 100% forest cover
sed_in_tons1 = np.zeros(n)
sed_in_tons2 = np.zeros(n)
sed_in_tons3 = np.zeros(n)
sed_in_tons4 = np.zeros(n)

for z in range(0,n):
    if forcov1 > low_def_rate:
        forcov1 = forcov1 - low_def_rate
    else:
        forcov1 = 0
    if forcov2 > avg_def_rate:
        forcov2 = forcov2 - avg_def_rate
    else:
        forcov2 = 0
    if forcov3 >high_def_rate:
        forcov3 = forcov3 - high_def_rate
    else:
        forcov3 = 0
    if forcov4 >forcov2015:
        forcov4 = forcov4 - low_def_rate
    else:
        forcov4 = forcov2015    
        
    sed_1 = a*math.exp(b*forcov1)
    sed_2 = a*math.exp(b*forcov2)
    sed_3 = a*math.exp(b*forcov3)
    sed_4 = a*math.exp(b*forcov4)
    
    sed_in_tons1[z] = sed_1
    sed_in_tons2[z] = sed_2
    sed_in_tons3[z] = sed_3
    sed_in_tons4[z] = sed_4
    
#Sediment in tons for the four scenarios
T = pd.DataFrame({'sed_in_tons1' : sed_in_tons1,
                  'sed_in_tons2' : sed_in_tons2,
                  'sed_in_tons3' : sed_in_tons3,
                  'sed_in_tons4' : sed_in_tons4},index = year)


#Sediments, in m^3
sed_in_cm1 = (sed_in_tons1 / sed_bulk_den)
sed_in_cm2 = (sed_in_tons2 / sed_bulk_den)
sed_in_cm3 = (sed_in_tons3 / sed_bulk_den)
sed_in_cm4 = (sed_in_tons4 / sed_bulk_den)

#Defining initial residence time, in days
initial_resi_time = initial_res_vol / (design_discharge * 3600 * 24);

#Defining trapping efficiency over 100 years, in %
initial_trap_eff = (1 - 0.05 * alpha_trap / math.sqrt(initial_resi_time)) * 100;

sed_tot1 = 0;
sed_tot2 = 0;
sed_tot3 = 0;
sed_tot4 = 0;

sedtot_1 = []
pcleft_1 = []
resvol_1 = []
resit_1 = []
trapeff_1 = []
         
sedtot_2 = []
pcleft_2 = []
resvol_2 = []
resit_2 = []
trapeff_2 = []
         
sedtot_3 = []
pcleft_3 = []
resvol_3 = []
resit_3 = []
trapeff_3 = []

sedtot_4 = []
pcleft_4 = []
resvol_4 = []
resit_4 = []
trapeff_4 = []

trap_eff1 = 0
trap_eff2 = 0
trap_eff3 = 0
trap_eff4 = 0

for i in range(0,n):
   #Defining total settled sediments every year, in m^3
    if i == 0:
        sed_settled1 = sed_in_cm1[i] * (initial_trap_eff / 100)
        sed_tot1 = sed_settled1 + sed_tot1
        sed_total1 = sed_tot1
        sed_settled2 = sed_in_cm2[i] * (initial_trap_eff / 100)
        sed_tot2 = sed_settled2 + sed_tot2
        sed_total2 = sed_tot2
        sed_settled3 = sed_in_cm3[i] * (initial_trap_eff / 100)
        sed_tot3 = sed_settled3 + sed_tot3
        sed_total3 = sed_tot3
        sed_settled4 = sed_in_cm4[i] * (initial_trap_eff / 100)
        sed_tot4 = sed_settled4 + sed_tot4
        sed_total4 = sed_tot4
    else:
        sed_settled1 = sed_in_cm1[i] * (trap_eff1 / 100)
        sed_tot1 = sed_settled1 + sed_tot1
        sed_total1 = sed_tot1
        sed_settled2 = sed_in_cm2[i] * (trap_eff2 / 100)
        sed_tot2 = sed_settled2 + sed_tot2
        sed_total2 = sed_tot2
        sed_settled3 = sed_in_cm3[i] * (trap_eff3 / 100)
        sed_tot3 = sed_settled3 + sed_tot3
        sed_total3 = sed_tot3
        sed_settled4 = sed_in_cm4[i] * (trap_eff4 / 100)
        sed_tot4 = sed_settled4 + sed_tot4
        sed_total4 = sed_tot4
        
    #Defining percent change in reservoir volume, in %
    if (initial_res_vol - sed_total1) / initial_res_vol < 0:
        pc_left1 = 0.1
    else:
        pc_left1 = (initial_res_vol - sed_total1) / initial_res_vol *100
        
    if (initial_res_vol -sed_total2) / initial_res_vol < 0:
        pc_left2 = 0.1
    else:
        pc_left2 = (initial_res_vol - sed_total2) / initial_res_vol *100
    
    if (initial_res_vol - sed_total3) / initial_res_vol < 0:
        pc_left3 = 0.1
    else:
        pc_left3 = (initial_res_vol - sed_total3) / initial_res_vol *100
    
    if (initial_res_vol - sed_total4) / initial_res_vol < 0:
        pc_left4 = 0.1
    else:
        pc_left4 = (initial_res_vol - sed_total4) / initial_res_vol *100
    
    #Defining reservoir volume change over 100 years, in m^3
    if sed_total1 < low_supply_vol:
        res_vol1 = initial_res_vol
    else:
        res_vol1 = initial_res_vol*pc_left1/100
        
    if sed_total2 < low_supply_vol:
        res_vol2 = initial_res_vol
    else:
        res_vol2 = initial_res_vol*pc_left2/100
    
    if sed_total3 < low_supply_vol:
        res_vol3 = initial_res_vol
    else:
        res_vol3 = initial_res_vol*pc_left3/100
    
    if sed_total4 < low_supply_vol:
        res_vol4 = initial_res_vol
    else:
        res_vol4 = initial_res_vol*pc_left4/100
    
    #Defining residence time over the 100 year timeframe
    if res_vol1/(design_discharge*3600*24) < 0:
        resi_time1 = 0
    else:
        resi_time1 = res_vol1/(design_discharge*3600*24)
     
    if res_vol2/(design_discharge*3600*24) < 0:
        resi_time2 = 0
    else:
        resi_time2 = res_vol2/(design_discharge*3600*24)
     
    if res_vol3/(design_discharge*3600*24) < 0:
        resi_time3 = 0
    else:
        resi_time3 = res_vol3/(design_discharge*3600*24)
     
    if res_vol4/(design_discharge*3600*24) < 0:
        resi_time4 = 0
    else:
        resi_time4 = res_vol4/(design_discharge*3600*24)
    
    #Defining trapping efficincies over the 100 year timeframe
    trap_eff1 =  (1-alpha_trap*0.05/math.sqrt(resi_time1))*100
    trap_eff2 =  (1-alpha_trap*0.05/math.sqrt(resi_time2))*100
    trap_eff3 =  (1-alpha_trap*0.05/math.sqrt(resi_time3))*100
    trap_eff4 =  (1-alpha_trap*0.05/math.sqrt(resi_time4))*100
    
    sedtot_1.append(sed_total1)
    pcleft_1.append(pc_left1)
    resvol_1.append(res_vol1)
    resit_1.append(resi_time1)
    trapeff_1.append(trap_eff1)
         
    sedtot_2.append(sed_total2)
    pcleft_2.append(pc_left2)
    resvol_2.append(res_vol2)
    resit_2.append(resi_time2)
    trapeff_2.append(trap_eff2)
         
    sedtot_3.append(sed_total3)
    pcleft_3.append(pc_left3)
    resvol_3.append(res_vol3)
    resit_3.append(resi_time3)
    trapeff_3.append(trap_eff3)

    sedtot_4.append(sed_total4)
    pcleft_4.append(pc_left4)
    resvol_4.append(res_vol4)
    resit_4.append(resi_time4)
    trapeff_4.append(trap_eff4)
    
sed_tot1 = sedtot_1
pc_left1 = pcleft_1
res_vol_left1 = resvol_1
resi_time1 = resit_1
trap_eff1 = trapeff_1

sed_tot2 = sedtot_2
pc_left2 = pcleft_2
res_vol_left2 = resvol_2
resi_time2 = resit_2
trap_eff2 = trapeff_2

sed_tot3 = sedtot_3
pc_left3 = pcleft_3
res_vol_left3 = resvol_3
resi_time3 = resit_3
trap_eff3 = trapeff_3

sed_tot4 = sedtot_4
pc_left4 = pcleft_4
res_vol_left4 = resvol_4
resi_time4 = resit_4
trap_eff4 = trapeff_4

#Tables showing sediments, percentage of reservoir volume left, actual
#volume of reservoir left, residence time and trapping efficiency over n
#years
T1 = pd.DataFrame({'sed_tot1' : sed_tot1,
                  'pc_left1' : pc_left1,
                  'res_vol_left1' : res_vol_left1,
                  'resi_time1' : resi_time1,
                  'trap_eff1' : trap_eff1},index = year)

T2 = pd.DataFrame({'sed_tot2' : sed_tot2,
                  'pc_left2' : pc_left2,
                  'res_vol_left2' : res_vol_left2,
                  'resi_time2' : resi_time2,
                  'trap_eff2' : trap_eff2},index = year)

T3 = pd.DataFrame({'sed_tot3' : sed_tot3,
                  'pc_left3' : pc_left3,
                  'res_vol_left3' : res_vol_left3,
                  'resi_time3' : resi_time3,
                  'trap_eff3' : trap_eff3},index = year)

T4 = pd.DataFrame({'sed_tot4' : sed_tot4,
                  'pc_left4' : pc_left4,
                  'res_vol_left4' : res_vol_left4,
                  'resi_time4' : resi_time4,
                  'trap_eff4' : trap_eff4},index = year)
    

SCENARIO1_Controlled_Deforestation = T1
SCENARIO2_Unmanaged_Deforestation = T2
SCENARIO3_Excessive_Deforestation = T3
SCENARIO4_Conservation = T4

T5 = pd.DataFrame({'sed_tot1' : sed_tot1,
                  'sed_tot2' : sed_tot2,
                  'sed_tot3' : sed_tot3,
                  'sed_tot4' : sed_tot4},index = year)

T6 = pd.DataFrame({'pc_left1' : pc_left1,
                  'pc_left2' : pc_left2,
                  'pc_left3' : pc_left3,
                  'pc_left4' : pc_left4},index = year)


####Comparing Scenarios########################################################
#Graphing watershed sediment yield
#Figure 1
fig, ax = plt.subplots()
plt.plot(year, pc_left1, year, pc_left2, year, pc_left3, year, pc_left4)
plt.title('Reservoir Volume Percent Change vs Time')
plt.xlabel('Year')
plt.ylabel('Percent Change')
plt.legend(('Controlled Deforestation', 'Unmanaged Deforestation', 'Excessive Deforestation', 'Conservation'))
plt.show()

#Figure 2
fig, ax = plt.subplots()
plt.plot(year, sed_tot1, year, sed_tot2, year, sed_tot3, year, sed_tot4)
plt.title('Sediment Accumulation vs Time')
plt.xlabel('Year')
plt.ylabel('Volume of Sediments in m^3')
plt.legend(('Controlled Deforestation', 'Unmanaged Deforestation', 'Excessive Deforestation', 'Conservation'))
plt.show()

#Power Value Loss Calculations
"""
Defining change in reservoir head in m
As sediment volume surpasses sediment volume at a certain water elevation,
the reservoir head below that elevation is considered lost
"""
displaced_head1 = 0
resvol_change1 = np.zeros((1,len(res_volume1)))
rhead1 = []
for i in range(0,n):
    for j in range(0,len(res_volume1)):
        if sed_tot1[i] > res_volume1[j]:
            res_volume1[j] = 0
        resvol_change1[0,:] = res_volume1    
        
    if resvol_change1.sum() > 0:
        current_resvol1 = np.min(resvol_change1[np.nonzero(resvol_change1)])
        sed_elev1 = c*(current_resvol1)**(d)
        displaced_head1 = sed_elev1 - low_supply_elev
        if displaced_head1 > res_height:
            new_head1 = 0
        else:
            new_head1 = rated_head-displaced_head1
    else:
        new_head1 = 0
    res_volume1 = copy.deepcopy(res_volume)
    rhead1.append(new_head1)
res_head1= rhead1

displaced_head2 = 0
resvol_change2 = np.zeros((1,len(res_volume2)))
rhead2 = []
for i in range(0,n):
    for j in range(0,len(res_volume2)):
        if sed_tot2[i] > res_volume2[j]:
            res_volume2[j] = 0
        resvol_change2[0,:] = res_volume2
        
    if resvol_change2.sum() > 0:
        current_resvol2 = np.min(resvol_change2[np.nonzero(resvol_change2)])
        sed_elev2 = c*(current_resvol2)**(d)
        displaced_head2 = sed_elev2 - low_supply_elev
        if displaced_head2 > res_height:
            new_head2 = 0
        else:
            new_head2 = rated_head-displaced_head2
    else:
        new_head2 = 0
    res_volume2 = copy.deepcopy(res_volume)
    rhead2.append(new_head2)
res_head2= rhead2

displaced_head3 = 0
resvol_change3 = np.zeros((1,len(res_volume3)))
rhead3 = []
for i in range(0,n):
    for j in range(0,len(res_volume3)):
        if sed_tot3[i] > res_volume3[j]:
            res_volume3[j] = 0
        resvol_change3[0,:] = res_volume3    
        
    if resvol_change3.sum() > 0:
        current_resvol3 = np.min(resvol_change3[np.nonzero(resvol_change3)])
        sed_elev3 = c*(current_resvol3)**(d)
        displaced_head3 = sed_elev3 - low_supply_elev
        if displaced_head3 > res_height:
            new_head3 = 0
        else:
            new_head3 = rated_head-displaced_head3
    else:
        new_head3 = 0
    res_volume3 = copy.deepcopy(res_volume)
    rhead3.append(new_head3)
res_head3= rhead3

displaced_head4 = 0
resvol_change4 = np.zeros((1,len(res_volume4)))
rhead4 = []
for i in range(0,n):
    for j in range(0,len(res_volume4)):
        if sed_tot4[i] > res_volume4[j]:
            res_volume4[j] = 0
        resvol_change4[0,:] = res_volume4      
        
    if resvol_change4.sum() > 0:
        current_resvol4 = np.min(resvol_change4[np.nonzero(resvol_change4)])
        sed_elev4 = c*(current_resvol4)**(d)
        displaced_head4 = sed_elev4 - low_supply_elev
        if displaced_head4 > res_height:
            new_head4 = 0
        else:
            new_head4 = rated_head-displaced_head4
    else:
        new_head4 = 0
    res_volume4 = copy.deepcopy(res_volume)
    rhead4.append(new_head4)
res_head4= rhead4

#Plotting reservoir head change in all four scenarios
T7 = pd.DataFrame({'res_head1' : res_head1,
                  'res_head2' : res_head2,
                  'res_head3' : res_head3,
                  'res_head4' : res_head4},index = year)

#Graphing change in reservoir head
#Figure 3
fig, ax = plt.subplots()
plt.plot(year, res_head1, year, res_head2, year, res_head3, year, res_head4)
plt.title('Reservoir Head Change vs Time')
plt.xlabel('Year')
plt.ylabel('Reservoir Head Change(m)')
plt.legend(('Controlled Deforestation', 'Unmanaged Deforestation', 'Excessive Deforestation', 'Conservation'))
plt.show()

#Defining power in GW and annual energy generation in GWh
power1 = [mu*den_water*g*design_discharge*x/1E9 for x in res_head1]
power2 = [mu*den_water*g*design_discharge*x/1E9 for x in res_head2]
power3 = [mu*den_water*g*design_discharge*x/1E9 for x in res_head3]
power4 = [mu*den_water*g*design_discharge*x/1E9 for x in res_head4]
gener1 = [8760*cap_factor*x for x in power1]
gener2 = [8760*cap_factor*x for x in power2]
gener3 = [8760*cap_factor*x for x in power3]
gener4 = [8760*cap_factor*x for x in power4]

#Defining annual revenue in USD
rev1 = np.asarray([elec*x for x in gener1])
rev2 = np.asarray([elec*x for x in gener2])
rev3 = np.asarray([elec*x for x in gener3])
rev4 = np.asarray([elec*x for x in gener4])

#Graphing change in revenue over time
#Figure 4
fig, ax = plt.subplots()
plt.plot(year,rev1, year, rev2, year, rev3, year, rev4)
plt.title('Annual Revenue vs Time')
plt.xlabel('Year')
plt.ylabel('Annual revenue in $')
plt.legend(('Controlled Deforestation', 'Unmanaged Deforestation', 'Excessive Deforestation', 'Conservation'))
plt.show()

#Defining annual loss in electricity selling revenue due to loss in reservoir storage
ann_loss_value1 = rev4 - rev1
ann_loss_value2 = rev4 - rev2
ann_loss_value3 = rev4 - rev3
ann_loss_value4 = rev4 - rev4

T8 = pd.DataFrame({'rev1' : rev1,
                  'rev2' : rev2,
                  'rev3' : rev3,
                  'rev4' : rev4},index = year)
    
T9 = pd.DataFrame({'ann_loss_value1' : ann_loss_value1,
                  'ann_loss_value2' : ann_loss_value2,
                  'ann_loss_value3' : ann_loss_value3,
                  'ann_loss_value4' : ann_loss_value4},index = year)

#Graphing annual loss value
#Figure 5
fig, ax = plt.subplots()
plt.plot(year,ann_loss_value1, year,ann_loss_value2, year,ann_loss_value3, year,ann_loss_value4)
plt.title('Annual Value Loss vs Time')
plt.xlabel('Year')
plt.ylabel('Annual Loss Value in $')
plt.legend(('Controlled Deforestation', 'Unmanaged Deforestation', 'Excessive Deforestation', 'Conservation'))
plt.show()

#Defining present value of lost energy production in USD
pv1 = []
pv2 = []
pv3 = []
pv4 = []

for i in range(0,n):
    pval1 = ann_loss_value1[i]/((1+disc)**(i+1))
    pval2 = ann_loss_value2[i]/((1+disc)**(i+1))
    pval3 = ann_loss_value3[i]/((1+disc)**(i+1))
    pval4 = ann_loss_value4[i]/((1+disc)**(i+1))
    pv1.append(pval1)
    pv2.append(pval2)
    pv3.append(pval3)
    pv4.append(pval4)

#note: in original code, this was T3 - renamed to T10 since values are different
T10 = pd.DataFrame({'pv1' : pv1,
                  'pv2' : pv2,
                  'pv3' : pv3,
                  'pv4' : pv4},index = year)
    
#note: original plot did not have labels
fig, ax = plt.subplots()
plt.plot(year,pv1,year,pv2,year,pv3,year,pv4)
plt.title('Present Value of lost Energy Production vs Time')
plt.xlabel('Year')
plt.ylabel('Present Value of lost Energy Production in $')
plt.legend(('Controlled Deforestation', 'Unmanaged Deforestation', 'Excessive Deforestation', 'Conservation'))
plt.show()

#Defining net present value in USD
npv1 = sum(pv1)
npv2 = sum(pv2)
npv3 = sum(pv3)
npv4 = sum(pv4)

#Defining value of forests to Hydropower
pes_ann1 = npv1/impacting_area
pes_ann2 = npv2/impacting_area
pes_ann3 = npv3/impacting_area
pes_ann4 = npv4/impacting_area
