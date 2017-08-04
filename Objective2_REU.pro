;
;AUTHOR: Bailey Donaldson
;
;
;PURPOSE: Objective 2 in the REU packet for SIMA in the VIS
;         
;

function calculate_prism_deg, Inst_mode_ID
 
;********************STEP 1: Read in everything and put data in the correct format********************
  
  ;contains caluclated f_function, SIM_data, NRLSSI_data, wavelength grid, and Sorce Days
  if Inst_mode_ID eq 41 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/VIS_data_f_function.sav'
    data = vis_data
    restore, '/Users/bado8599/Desktop/REU/sav_files/VIS_daily_kappa.sav'
  endif
  
  if Inst_mode_ID eq 43 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/UV_data_f_function.sav'
    data = uv_data
    restore, '/Users/bado8599/Desktop/REU/sav_files/UV_daily_kappa.sav'
  endif
  
  if Inst_mode_ID eq 44 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/IR_data_f_function.sav'
    data = ir_data
    restore, '/Users/bado8599/Desktop/REU/sav_files/IR_daily_kappa.sav'
  endif
  
  ;solar exposure
  restore, '/Users/bado8599/Desktop/REU/sav_files/solarexp_54_55_raw.sav'
  solar54.solar_exp /= 86400d
  
  ;kappa
  readcol, '/Users/bado8599/Desktop/REU/Bailey/ESR_kappa_453_1570_2011_mod_407nm.txt',$
     kappa_w, kappa_k, format = '(d,d)', /silent
  kappa = {wavelength: kappa_w, k:kappa_k}

  ;attenuation (a)
  readcol, '/Users/bado8599/Desktop/REU/sav_files/overlap_vis1_mod_407nm.txt',$
     ray_w, ray_a, format = '(d,d)', /silent
  attenuation = {wavelength: ray_w, a: ray_a}
  
  ;adjust values of T1 to SD from GPS time in micro seconds
  solar_exp_SD = gps2sd(solar54.T1/1d6)
 
 
;*************************STEP 2: Calculate new prism degredation**************************
  
  ;creates empty pointer 
  prism_deg = ptrarr(n_elements(data))

  for j=0, n_elements(data)-1 do begin

    ;is interpolated within the loop because the amount of elements in *data changes every measurement
    attenuation1 = interpol(attenuation.a, attenuation.wavelength, (*data[j]).wavelength)
    kappa1 = interpol(kappa.k, kappa.wavelength, (*data[j]).wavelength)

    ;finds where timestamp in solar_exp is closest to SIM timestamp to find which solar_exp to use in the equation
    hello = min(abs(solar_exp_SD - (*data[j]).sim_timestamp[0]), pos_time)

    ;this is initialized within the first for loop so that at every Sorce Day, a new calculated prism_deg is created,
    ;otherwise it will create one one-dimmensional array of prism_deg for every Sorce Day at every wavelength
    calc_prismdeg_avg_kappa =[]
    calc_prismdeg_daily_kappa = []
    
    for i=0, n_elements((*data[j]).wavelength)-1 do begin
      
      ;contains prism degredation at one SD (the current value of j) for every wavelength using the avg kappa and my calculated f function
      calc_prismdeg_avg_kappa = [calc_prismdeg_avg_kappa, (1d -attenuation1[i])*exp(-kappa1[i]*solar54[pos_time].solar_exp*((*data[j]).f_function[i])) +$
          attenuation1[i]*exp(-kappa1[i]*solar54[pos_time].solar_exp*((*data[j]).f_function[i])/2d)]
          
      ;contains prism degredation at one SD (the current value of j) for every wavelength using my calculated daily kappa and a f function of 1.0
      calc_prismdeg_daily_kappa = [calc_prismdeg_daily_kappa, (1d - attenuation1[i])*exp(-(*daily_kappa[j]).kappa[i]*solar54[pos_time].solar_exp*1d) +$
          attenuation1[i]*exp(-(*daily_kappa[j]).kappa[i]*solar54[pos_time].solar_exp*(1d/2d))]
         
      prism_deg[j] = ptr_new({prism_deg_avg_kappa: calc_prismdeg_avg_kappa, prism_deg_daily_kappa: calc_prismdeg_daily_kappa})
      
    endfor
  endfor
  
  if Inst_mode_ID eq 41 then begin
    save, prism_deg, filename = '/Users/bado8599/Desktop/REU/sav_files/prism_deg_VIS.sav'
  endif
  if Inst_mode_ID eq 43 then begin
    save, prism_deg, filename = '/Users/bado8599/Desktop/REU/sav_files/prism_deg_UV.sav'
  endif  
  if Inst_mode_ID eq 44 then begin
    save, prism_deg, filename = '/Users/bado8599/Desktop/REU/sav_files/prism_deg_IR.sav'
  endif  
end
 
 
;********************STEP 3: plot prism degredation (objective 2a)*************************

;need to have run convolve_nrlssi.pro, calculate_f.pro, and calculate_prism_deg.pro to have the
  ;necessary files to run this function
  
function plot_prism_deg, Inst_mode_ID

  if Inst_mode_ID eq 41 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/VIS_data_f_function.sav'
    data = vis_data
    title_name = 'VIS'
    restore, '/Users/bado8599/Desktop/REU/sav_files/prism_deg_VIS.sav'
  endif
  
  if Inst_mode_ID eq 43 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/UV_data_f_function.sav'
    data = uv_data
    title_name = 'UV'
    restore, '/Users/bado8599/Desktop/REU/sav_files/prism_deg_UV.sav'    
  endif
  if Inst_mode_ID eq 44 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/IR_data_f_function.sav'
    data = ir_data
    title_name = 'IR'
    restore, '/Users/bado8599/Desktop/REU/sav_files/prism_deg_IR.sav'    
  endif
  
  
  choice = ''
  READ, choice, PROMPT = 'Plot prism degredation with constant WL or SD? '

  if strlowcase(choice) eq 'wl' then begin
    ;user_wl = #
    READ, user_wl, PROMPT = 'Choose a Wavelength: '
    wl = 0
    SIM_timestamp = []
    prism_deg_at_chosen_wl = []
    for k=0, n_elements(data)-1 do begin
      ;puts the timestamp into a format we can use to graph (takes it out of the pointer)
      SIM_timestamp = [SIM_timestamp, (*data[k]).SIM_timestamp]
      
      ;finds the data that has the position in the array is closest to the wavelength the user chose
      ;it is looped because the wavelength's position will be different for each Sorce Day
      hello = min(abs((*data[k]).wavelength - user_wl), pos_wl)
      ;creates an array of the prism degradation at each SD at a constant wavelength
      prism_deg_at_chosen_wl = [prism_deg_at_chosen_wl, (*prism_deg[k]).prism_deg_avg_kappa[pos_wl]]
    endfor   
    
    prismdeg_no_outliers = []
    SIM_timestamp_no_outliers = []
    ;try to remove the data outliers
    prismdeg_mean = total(prism_deg_at_chosen_wl)/n_elements(prism_deg_at_chosen_wl)
    ;the numbers below may need to be adjusted depending on the wavelength graphed
    for m=0, n_elements(prism_deg_at_chosen_wl)-1 do begin
      ;the distance from the mean that is permitted may need to be altered
      if prism_deg_at_chosen_wl[m] ge prismdeg_mean-(prismdeg_mean*0.5) and prism_deg_at_chosen_wl[m] le prismdeg_mean+(prismdeg_mean*0.87) then begin
        ;the difference permitted between each number may have to be altered
        if prism_deg_at_chosen_wl[m] ge prism_deg_at_chosen_wl[m-1]-1 and prism_deg_at_chosen_wl[m] le prism_deg_at_chosen_wl[m-1]+1 then begin
          prismdeg_no_outliers = [prismdeg_no_outliers, prism_deg_at_chosen_wl[m]]
          SIM_timestamp_no_outliers = [SIM_timestamp_no_outliers, SIM_timestamp[m]]
        endif
      endif
      
    endfor
    
    ;turns user_wl into a string so it can be used in the title
    hello = min(abs((*data[100]).wavelength - user_wl), pos_wl)
    
    wl_string = strtrim((*data[100]).wavelength[pos_wl], 1)
    
    ;keeps graph from turing red
    device, decompose = 0
   
    
    ;plot of all of the prismdeg values
    lineplot, SIM_timestamp, prism_deg_at_chosen_wl, xtitle = 'Sorce Day', ytitle = 'Prism Degradation', $
      title = 'Sorce Day vs. Prism Degradation at Wavelength: '+wl_string + 'nm', $
      ptitle = 'Sorce Day vs. Prism Degradation at Constant Wavelength for SIMA in the ' + title_name
      
    ;plot of the prism deg values without the outliers
    lineplot, SIM_timestamp_no_outliers, prismdeg_no_outliers, title = 'Sorce Day vs. Prism Degradation at Wavelength: '$
      +wl_string + 'nm without outliers'
    
    
  endif

  if strlowcase(choice) eq 'sd' then begin  
    READ, user_sd, PROMPT = 'Choose a Mission Timestamp (Sorce Day): ' 
    
    ;takes SIM timestamps out of the pointer and puts it into an array
    SIM_timestamp = []
    for k=0, n_elements(data)-1 do begin
      SIM_timestamp = [SIM_timestamp, (*data[k]).SIM_timestamp]
    endfor 
    
    ;user_sd = #   
    ;find timestamp closest to the timestamp the user chose  
    hello = min(abs(SIM_timestamp - user_sd), pos_SD)
    
    ;turns SD into a string so that it can be used in the title
    SD_string = strtrim((*data[pos_SD]).SIM_timestamp, 1)
    
    device, decompose = 0
    lineplot, (*data[pos_SD]).wavelength, (*prism_deg[pos_SD]).prism_deg1, xtitle = 'Wavelength', $
      ytitle = 'Prism Degradation', title = 'Wavelength vs. Prism Degradation at SD ' + SD_string, $
      ptitle = 'Wavelength vs. Prism Degradation at SD ' + SD_string + 'for SIM A in the ' + title_name

  endif

end

;***************STEP 4: plot the Sorce Day vs. solar exposure record (objective 2b)**************
;****I inverted the order from Objective 2*****

function plot_solarexp

  ;solar exposure (solar54 is SIMA)
  restore, '/Users/bado8599/Desktop/REU/sav_files/solarexp_54_55_raw.sav'
  
  ;adjusts values of solar_exposure to be in days instead of seconds
  solar54.solar_exp /= 86400d
  
  ;adjust values of T1 to SD from GPS time in micro seconds
  solar_exp_SD = gps2sd(solar54.T1/1d6)
  
  
  ;contains the SIM timestamps
  ;the VIS, IR, and UV, have different amounts of timestamps because I eliminated the timestamps containing data that was insufficient
    ;VIS had the most amount of measurements that were retainted (after I filtered them in convolve_nrlssi.pro) so it will be 
    ;used in this code
;  restore, '/Users/bado8599/Desktop/REU/sav_files/VIS_data_f_function.sav'
;  data = vis_data

;use this code if you want to use the first graph  
;  SIM_timestamp = []
;  solar_exp = []
;  pos_time = 0
;  for k=0, n_elements(data)-1 do begin
;    
;    ;takes SIM timestamps from the pointer/structure and puts it into an array so it can be graphed
;    SIM_timestamp = [SIM_timestamp, (*data[k]).SIM_timestamp]
;    
;    ;finds where timestamp in solar_exp is closest to SIM timestamp
;    hello = min(abs(solar_exp_SD - (*data[k]).SIM_timestamp[0]), pos_time)
;    
;    ;creates an array containing the solar_exposures at each of the Sorce Days in SIM_timestamp
;    solar_exp = [solar_exp, solar54[pos_time].solar_exp]
;    
;  endfor

  device, decompose = 0
  ;this plot only graphs the solar exposures at each time SIM recorded data
;  lineplot, SIM_timestamp, solar_exp, ytitle = 'Solar Exposure (Days)', $
;    xtitle = 'SIM timestamp (Sorce Day)', title = 'SIM timestamp (Sorce Days) vs. Solar Exposure (Days)', $
;    ptitle = 'SIM timestamp (Sorce Days) vs. Solar Exposure (Days) for SIMA'
   
  ;this plot records all solar esposure measurements   
  lineplot, solar_exp_SD, solar54.solar_exp, title = 'Sorce Day vs. Solar Exposure (Days)', $
    ytitle = 'Solar Exposure (Days)', xtitle = 'Sorce Day', ptitle = 'Sorce Day vs. Solar Exposure'
  
    
end

;********************STEP 5: plot the Sorce day vs. f' for SIMA (Objective 2c) ***********************
;********NOTE: the objective instructed to plot f vs. SD but I inverted the order************

function plot_f, Inst_mode_ID

  ;contains SIM timestamps and f_function values
  if Inst_mode_ID eq 41 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/VIS_data_f_function.sav'
    data = vis_data
    type = 'VIS'
  endif
  
  if Inst_mode_ID eq 43 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/UV_data_f_function.sav'
    data = uv_data
    type = 'UV'    
  endif
  
  if Inst_mode_ID eq 44 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/IR_data_f_function.sav'
    data = ir_data
    type = 'IR'    
  endif
  
  READ, user_wl, PROMPT = 'Which wavelength do you want to plot the f function? '
  ;user_wl = #
  ;puts wl into a string so we can use it in the graph= title
  user_wl_string = strtrim(user_wl, 1)
  
  SIM_timestamp = []
  f_function_at_chosen_wl = []
  pos_wl = 0
  for l=0, n_elements(data)-1 do begin
    
    ;takes SIM timestamp out of the pointer and puts it into an array that can more easily be graphed
    SIM_timestamp = [SIM_timestamp, (*data[l]).SIM_timestamp]
    
    ;finds the data that has the position in the array is closest to the wavelength the user chose
    ;it is looped because the wavelength's position will be different for each Sorce Day
    hello = min(abs((*data[l]).wavelength - user_wl), pos_wl)
    
    ;creates an array of the f function at each SD at a constant wavelength
    f_function_at_chosen_wl = [f_function_at_chosen_wl, (*data[l]).f_function[pos_wl]]
    
  endfor

  f_function_no_outliers = []
  SIM_timestamp_no_outliers = []
;try to remove the data outliers
  f_mean = total(f_function_at_chosen_wl)/n_elements(f_function_at_chosen_wl)
  if f_mean lt 0 then f_mean = (-1)*f_mean
  if f_mean lt 0.0 then f_mean = f_mean*(-1.0)
 ;the numbers below may need to be adjusted depending on the wavelength graphed
  for m=0, n_elements(f_function_at_chosen_wl)-1 do begin
   ;the distance from the mean that is permitted may need to be altered
   if f_function_at_chosen_wl[m] ge f_mean-(f_mean*100) and f_function_at_chosen_wl[m] le f_mean+(f_mean*20) then begin
    ;the difference permitted between each number may have to be altered (usualy is about 0.1)
    if f_function_at_chosen_wl[m] ge f_function_at_chosen_wl[m-1]-.5 and f_function_at_chosen_wl[m] le f_function_at_chosen_wl[m-1]+.5 then begin
      if f_function_at_chosen_wl[m] ne -1 then begin
        f_function_no_outliers = [f_function_no_outliers, f_function_at_chosen_wl[m]]
        SIM_timestamp_no_outliers = [SIM_timestamp_no_outliers, SIM_timestamp[m]]
        endif
      endif
     endif 
  endfor
 
  hello = min(abs((*data[100]).wavelength - user_wl), pos_wl)
  wl_string = strtrim((*data[100]).wavelength[pos_wl], 1)

  
 ;graph of all of the data points
  device, decompose = 0
  lineplot, SIM_timestamp, f_function_at_chosen_wl, xtitle = 'SIM timestamp (Sorce Days)', ytitle = 'F Function',$
    ptitle = "F function vs. Sorce Day in the " + type 
     
  ;graph with the outliers removed
  device, decompose = 0
  lineplot, SIM_timestamp_no_outliers, f_function_no_outliers, title = "F function vs Sorce Day in the " +$
     type + " at " + wl_string + ' without outliers'
 
end


