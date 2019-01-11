
% ForPower
% Scripted by : Mohit Kaura
% Description: To estimate sediment trapping and accumulation in the reservoir

%% SEDIMENTATION

% Importing variables from DAmdata Excel Sheet
XX = 4;
damdata = xlsread('Damdata.xlsx','DamData');%Damdata Excel Sheet required
initial_res_vol = damdata(1,XX); %Initial Reservoir Volume(m3)
high_supp_elev = damdata(2,XX); %High Supply Level (m)
high_supply_vol = damdata(3,XX); %High Supply Volume (m3)
low_supply_elev = damdata(4,XX); %Low Supply Level (m)
low_supply_vol = damdata(5,XX); %Low Supply Volume (m3)
design_discharge = damdata(6,XX); %Design Discharge (m3/s)
mean_ann_energy = damdata(7,XX); %Mean Annual Energy (GWh)
low_def_rate= damdata(8,XX); %Low deforestation rate
avg_def_rate= damdata(9,XX); %Avg Deforestation rate
high_def_rate= damdata(10,XX); %High deforestation rate
rated_head= damdata(11,XX); %Rated Head (m)
sed_bulk_den = damdata(12,XX); %Bulk sediment density
alpha_trap = damdata(13,XX); %Alpha (TE)
a = damdata(14,XX); %constant a
b = damdata(15,XX); %constant b
c = damdata(16,XX); %constant c
d = damdata(17,XX); %constant d
g = damdata(18,XX); %Accln due to gravity (m/s2)
forcov1 = damdata(19,XX); %Forest cover: Scenario1 (percentage)
forcov2 = damdata(20,XX); %Forest cover: Scenario2 (percentage)
forcov3 = damdata(21,XX); %Forest cover: Scenario3 (percentage)
forcov4 = damdata(22,XX); %Forest cover: Scenario4 (percentage)
forcov2015 = damdata(23,XX); %Forest cover: 2015 (percentage)
sed1 = damdata(24,XX); % initial sediment Scenario1
sed2 = damdata(25,XX); % initial sediment Scenario2
sed3 = damdata(26,XX); % initial sediment Scenario3
sed4 = damdata(27,XX); % initial sediment Scenario4
den_water = damdata(28,XX); %Density of water (kg/m3)
mu = damdata(29,XX); % dam efficiency
cap_factor = damdata(30,XX); %capacity factor
elec = damdata(31,XX); %Electricty cost per Kwh ($)
disc = damdata(32,XX); %Discount Rate
impacting_area = damdata(33,XX); %Annuity
res_height = damdata(34,XX); %Reservoir Height

hvsv = flipud(xlsread('Damdata.xlsx','hvsv','G:H'));
water_height = hvsv(:,1); %Water Elevation
res_volume = hvsv(:,2);
res_volume1 = hvsv(:,2); %Reservoir volume at that Water Elevation
res_volume2 = hvsv(:,2);
res_volume3 = hvsv(:,2);
res_volume4 = hvsv(:,2);

% SCENARIO 1: Controlled Deforestation
% SCENARIO 2: Unmanaged Deforestation
% SCENARIO 3: Excessive Deforestation
% SCENARIO 4: Conservation and proper management
% EVERY TERM ENDING WITH 1 RELATES TO SCENARIO 1, 2 FOR SCENARIO 2 AND SO
% FORTH
year = [1:121]';
[n,m] = size(year);
% n is the length of year

% forcover = 100 shows 100% forest cover
for z = 1:n
    if forcov1 >low_def_rate
        forcov1 = forcov1 - low_def_rate;
    else
        forcov1 = 0;
    end
    if forcov2 >avg_def_rate
        forcov2 = forcov2 - avg_def_rate;
    else
        forcov2 = 0;
    end
    if forcov3 >high_def_rate
        forcov3 = forcov3 - high_def_rate;
    else
        forcov3 = 0;
    end
    if forcov4 >forcov2015
        forcov4 = forcov4 - low_def_rate;
    else
        forcov4 = forcov2015;
    end
    sed_1 = a*exp(b*forcov1);
    sed_2 = a*exp(b*forcov2);
    sed_3 = a*exp(b*forcov3);
    sed_4 = a*exp(b*forcov4);
sed_in_tons1(z,1) = [sed_1];
sed_in_tons2(z,1) = [sed_2];
sed_in_tons3(z,1) = [sed_3];
sed_in_tons4(z,1) = [sed_4];
end
[sed_in_tons1];
[sed_in_tons2];
[sed_in_tons3];
[sed_in_tons4];
T = table(sed_in_tons1, sed_in_tons2, sed_in_tons3, sed_in_tons4); % Sediment in tons for the four scenarios

% Sediments, in m^3
sed_in_cm1 = (sed_in_tons1 / sed_bulk_den);
sed_in_cm2 = (sed_in_tons2 / sed_bulk_den);
sed_in_cm3 = (sed_in_tons3 / sed_bulk_den);
sed_in_cm4 = (sed_in_tons4 / sed_bulk_den);

% Defining initial residence time, in days
initial_resi_time = initial_res_vol / (design_discharge * 3600 * 24);

% Defining trapping efficiency over 100 years, in %
initial_trap_eff = (1 - 0.05 * alpha_trap / sqrt(initial_resi_time)) * 100;

sed_tot1 = 0;
sed_tot2 = 0;
sed_tot3 = 0;
sed_tot4 = 0;
for i = 1:n
   % Defining total settled sediments every year, in m^3
    if i ==1
        sed_settled1 = sed_in_cm1(i) * (initial_trap_eff / 100);
        sed_tot1 = sed_settled1 +sed_tot1;
        sed_total1 = sed_tot1;
        sed_settled2 = sed_in_cm2(i) * (initial_trap_eff / 100);
        sed_tot2 = sed_settled2 +sed_tot2;
        sed_total2 = sed_tot2;
        sed_settled3 = sed_in_cm3(i) * (initial_trap_eff / 100);
        sed_tot3 = sed_settled3 +sed_tot3;
        sed_total3 = sed_tot3;
        sed_settled4 = sed_in_cm4(i) * (initial_trap_eff / 100);
        sed_tot4 = sed_settled4 +sed_tot4;
        sed_total4 = sed_tot4;
    else
        sed_settled1 = sed_in_cm1(i) * (trap_eff1 / 100);
        sed_tot1 = sed_settled1 +sed_tot1;
        sed_total1 = sed_tot1;
        sed_settled2 = sed_in_cm2(i) * (trap_eff2 / 100);
        sed_tot2 = sed_settled2 +sed_tot2;
        sed_total2 = sed_tot2;
        sed_settled3 = sed_in_cm3(i) * (trap_eff3 / 100);
        sed_tot3 = sed_settled3 +sed_tot3;
        sed_total3 = sed_tot3;
        sed_settled4 = sed_in_cm4(i) * (trap_eff4 / 100);
        sed_tot4 = sed_settled4 +sed_tot4;
        sed_total4 = sed_tot4;
    end
    % Defining percent change in reservoir volume, in %
    if (initial_res_vol -sed_total1) / initial_res_vol<0
        pc_left1 = 0.1;
    else
        pc_left1 = (initial_res_vol -sed_total1) / initial_res_vol *100;
    end
    if (initial_res_vol -sed_total2) / initial_res_vol<0
        pc_left2 = 0.1;
    else
        pc_left2 = (initial_res_vol -sed_total2) / initial_res_vol *100;
    end
    if (initial_res_vol -sed_total3) / initial_res_vol<0
        pc_left3 = 0.1;
    else
        pc_left3 = (initial_res_vol -sed_total3) / initial_res_vol *100;
    end
    if (initial_res_vol -sed_total4) / initial_res_vol<0
        pc_left4 = 0.1;
    else
        pc_left4 = (initial_res_vol -sed_total4) / initial_res_vol *100;
    end
    % Defining reservoir volume change over 100 years, in m^3
    if sed_total1<low_supply_vol
        res_vol1 = initial_res_vol;
    else
        res_vol1 = initial_res_vol*pc_left1/100;
    end
    if sed_total2<low_supply_vol
        res_vol2 = initial_res_vol;
    else
        res_vol2 = initial_res_vol*pc_left2/100;
    end
    if sed_total3<low_supply_vol
        res_vol3 = initial_res_vol;
    else
        res_vol3 = initial_res_vol*pc_left3/100;
    end
    if sed_total4<low_supply_vol
        res_vol4 = initial_res_vol;
    else
        res_vol4 = initial_res_vol*pc_left4/100;
    end
    % Defining residence time over the 100 year timeframe
     if res_vol1/(design_discharge*3600*24)<0
        resi_time1 = 0;
     else
         resi_time1 = res_vol1/(design_discharge*3600*24);
     end
     if res_vol2/(design_discharge*3600*24)<0
        resi_time2 = 0;
     else
         resi_time2 = res_vol2/(design_discharge*3600*24);
     end
     if res_vol3/(design_discharge*3600*24)<0
        resi_time3 = 0;
     else
         resi_time3 = res_vol3/(design_discharge*3600*24);
     end
     if res_vol4/(design_discharge*3600*24)<0
        resi_time4 = 0;
     else
         resi_time4 = res_vol4/(design_discharge*3600*24);
     end
    % Defining trapping efficincies over the 100 year timeframe
trap_eff1 =  (1-alpha_trap*0.05./sqrt(resi_time1))*100;
trap_eff2 =  (1-alpha_trap*0.05./sqrt(resi_time2))*100;
trap_eff3 =  (1-alpha_trap*0.05./sqrt(resi_time3))*100;
trap_eff4 =  (1-alpha_trap*0.05./sqrt(resi_time4))*100;

sedtot_1(:,i) = [sed_total1];
pcleft_1(:,i) = [pc_left1];
resvol_1(i,1) = [res_vol1];
resit_1(1,i)= [resi_time1];
trapeff_1(:,i)=[trap_eff1];
sedtot_2(:,i) = [sed_total2];
pcleft_2(:,i) = [pc_left2];
resvol_2(i,1) = [res_vol2];
resit_2(1,i)= [resi_time2];
trapeff_2(:,i)=[trap_eff2];
sedtot_3(:,i) = [sed_total3];
pcleft_3(:,i) = [pc_left3];
resvol_3(i,1) = [res_vol3];
resit_3(1,i)= [resi_time3];
trapeff_3(:,i)=[trap_eff3];
sedtot_4(:,i) = [sed_total4];
pcleft_4(:,i) = [pc_left4];
resvol_4(i,1) = [res_vol4];
resit_4(1,i)= [resi_time4];
trapeff_4(:,i)=[trap_eff4];
end
sed_tot1 = sedtot_1';
pc_left1 = pcleft_1';
res_vol_left1 = resvol_1;
resi_time1 = resit_1';
trap_eff1 = trapeff_1';
sed_tot2 = sedtot_2';
pc_left2 = pcleft_2';
res_vol_left2 = resvol_2;
resi_time2 = resit_2';
trap_eff2 = trapeff_2';
sed_tot3 = sedtot_3';
pc_left3 = pcleft_3';
res_vol_left3 = resvol_3;
resi_time3 = resit_3';
trap_eff3 = trapeff_3';
sed_tot4 = sedtot_4';
pc_left4 = pcleft_4';
res_vol_left4 = resvol_4;
resi_time4 = resit_4';
trap_eff4 = trapeff_4';
% Table showing sediments, percentage of reservoir volume left, actual
% volume of reservoir left, residence time and traping efficiency over n
% years
T1 = table(year, sed_tot1 , pc_left1, res_vol_left1, resi_time1, trap_eff1 );
T2 = table(year, sed_tot2 , pc_left2, res_vol_left2, resi_time2, trap_eff2 );
T3 = table(year, sed_tot3 , pc_left3, res_vol_left3, resi_time3, trap_eff3 );
T4 = table(year, sed_tot4 , pc_left4, res_vol_left4, resi_time4, trap_eff4 );
SCENARIO1_Controlled_Deforestation = T1;
SCENARIO2_Unmanaged_Deforestation = T2;
SCENARIO3_Excessive_Deforestation = T3;
SCENARIO4_Conservation = T4;
T5 = table(sed_tot1 , sed_tot2,sed_tot3, sed_tot4 )
T6 = table(pc_left1, pc_left2, pc_left3, pc_left4 )
% COMPARING SCENARIOS 

% Graphing watershed sedment yield
figure(1)
plot(year, pc_left1,year, pc_left2, year, pc_left3, year, pc_left4);
title('Reservoir Volume Percent Change vs Time')
xlabel('Year')
ylabel('Percent Change')
legend('Controlled Deforestation', 'Unmanaged Deforestation', 'Excessive Deforestation', 'Conservation')
figure(2)
plot(year, sed_tot1, year, sed_tot2, year, sed_tot3, year, sed_tot4);
title('Sediment Accumulation vs Time')
xlabel('Year')
ylabel('Volume of Sediments in m^3')
legend('Controlled Deforestation', 'Unmanaged Deforestation', 'Excessive Deforestation', 'Conservation')

% Power Value Loss Calculations

% Defining change in reservoir head in m
% As sediment volume surpasses sediment volume at a certain water elvation,
% the reservoir head below that elevation is considered lost
displaced_head1 = 0;
resvol_change1 = zeros(1,length(res_volume1));
for i = 1:n
    for j = 1:length(res_volume1)
        if sed_tot1(i)>res_volume1(j)
            res_volume1(j)=0;
        end
        resvol_change1(1,:)= res_volume1;        
    end
    if sum(resvol_change1)>0
        current_resvol1 = min(resvol_change1(resvol_change1>0));
        sed_elev1 = c*(current_resvol1).^(d);
        displaced_head1 = sed_elev1-low_supply_elev;
        if displaced_head1 > res_height
            new_head1 = 0;
        else
            new_head1 = rated_head-displaced_head1;
        end
    else
        new_head1 = 0;
    end
    res_volume1 = res_volume;
rhead1(i,1)=new_head1;
end
res_head1= rhead1;
displaced_head2 = 0;
resvol_change2 = zeros(1,length(res_volume2));
for i =1:n
    for j = 1:length(res_volume2)
        if sed_tot2(i)>res_volume2(j)
            res_volume2(j)=0;
        end
        resvol_change2(1,:)= res_volume2;
    end
    if sum(resvol_change2)>0
        current_resvol2 = min(resvol_change2(resvol_change2>0));
        sed_elev2 = c*(current_resvol2).^(d);
        displaced_head2 = sed_elev2-low_supply_elev;
        if displaced_head2 > res_height
            new_head2 = 0;
        else
            new_head2 = rated_head-displaced_head2;
        end
    else
        new_head2 = 0;
    end
    res_volume2 = res_volume;
rhead2(i,1)=new_head2;
end
res_head2= rhead2;
displaced_head3 = 0;
resvol_change3 = zeros(1,length(res_volume3));
for i = 1:n
    for j = 1:length(res_volume3)
        if sed_tot3(i)>res_volume3(j)
            res_volume3(j)=0;
        end
        resvol_change3(1,:)= res_volume3;        
    end
    if sum(resvol_change3)>0
        current_resvol3 = min(resvol_change3(resvol_change3>0));
        sed_elev3 = c*(current_resvol3).^(d);
        displaced_head3 = sed_elev3-low_supply_elev;
        if displaced_head3 > res_height
            new_head3 = 0;
        else
            new_head3 = rated_head-displaced_head3;
        end
    else
        new_head3 = 0;
    end
    res_volume3 = res_volume;
rhead3(i,1)=new_head3;
end
res_head3= rhead3;
displaced_head4 = 0;
resvol_change4 = zeros(1,length(res_volume4));
for i = 1:n
    for j = 1:length(res_volume4)
        if sed_tot4(i)>res_volume4(j)
            res_volume4(j)=0;
        end
        resvol_change4(1,:)= res_volume4;        
    end
    if sum(resvol_change4)>0
        current_resvol4 = min(resvol_change4(resvol_change4>0));
        sed_elev4 = c*(current_resvol4).^(d);
        displaced_head4 = sed_elev4-low_supply_elev;
        if displaced_head4 > res_height
            new_head4 = 0;
        else
            new_head4 = rated_head-displaced_head4;
        end
    else
        new_head4 = 0;
    end
    res_volume4 = res_volume;
rhead4(i,1)=new_head4;
end
res_head4= rhead4;
T7 = table(res_head1, res_head2, res_head3, res_head4 ) %Plotting reservoir head change in all four scenarios

% Graphing change in reservoir head
figure(3)
plot(year, res_head1,year, res_head2, year, res_head3, year, res_head4);
title('Reservoir Head Change vs Time')
xlabel('Year')
ylabel('Reservoir Head Change(m)')
legend('Controlled Deforestation', 'Unmanaged Deforestation', 'Excessive Deforestation', 'Conservation')


% Defining power in GW and annual energy generation in GWh
power1= mu*den_water*g*design_discharge*res_head1/1E9;
power2= mu*den_water*g*design_discharge*res_head2/1E9;
power3= mu*den_water*g*design_discharge*res_head3/1E9;
power4= mu*den_water*g*design_discharge*res_head4/1E9;
gener1 = power1*8760*cap_factor;
gener2 = power2*8760*cap_factor;
gener3 = power3*8760*cap_factor;
gener4 = power4*8760*cap_factor;
% Defining annual revenue in USD
rev1 = elec*gener1;
rev2 = elec*gener2;
rev3 = elec*gener3;
rev4 = elec*gener4;
% Graphing change in revenue over time
year;
figure(4)
plot(year,rev1, year, rev2, year, rev3, year, rev4)

title('Annual Revenue vs Time')
xlabel('Year')
ylabel('Annual revenue in $')
legend('Controlled Deforestation', 'Unmanaged Deforestation', 'Excessive Deforestation', 'Conservation')

% Defining annual loss in electricity selling revenues due to loss in
% reservoir storage
ann_loss_value1 = rev4-rev1;
ann_loss_value2 = rev4-rev2;
ann_loss_value3 = rev4-rev3;
ann_loss_value4 = rev4-rev4;
T8 = table(rev1 , rev2,rev3, rev4 )
T9 = table(ann_loss_value1, ann_loss_value2, ann_loss_value3, ann_loss_value4 )

% Graphing annual loss value
year;
figure(5)
plot(year,ann_loss_value1, year,ann_loss_value2, year,ann_loss_value3, year,ann_loss_value4)
title('Annual Value Loss vs Time')
xlabel('Year')
ylabel('Annual Loss Value in $')
legend('Controlled Deforestation', 'Unmanaged Deforestation', 'Excessive Deforestation', 'Conservation')

% Defining present value of lost energy production in USD
for i = 1:n
pval1 = ann_loss_value1(i)/((1+disc)^(i));
pval2 = ann_loss_value2(i)/((1+disc)^(i));
pval3 = ann_loss_value3(i)/((1+disc)^(i));
pval4 = ann_loss_value4(i)/((1+disc)^(i));
pv1(i,1)= pval1;
pv2(i,1)= pval2;
pv3(i,1)= pval3;
pv4(i,1)= pval4;
end
T3 = table(pv1,pv2,pv3, pv4 )
plot(year,pv1,year,pv2,year,pv3,year,pv4)

% Defining net present value in USD
npv1 = sum(pv1)
npv2 = sum(pv2)
npv3 = sum(pv3)
npv4 = sum(pv4)

% Defining value of forests to Hydropower
pes_ann1 = npv1/impacting_area;
pes_ann2 = npv2/impacting_area;
pes_ann3 = npv3/impacting_area;
pes_ann4 = npv4/impacting_area;













