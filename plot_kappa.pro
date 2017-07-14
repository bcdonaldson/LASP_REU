;
;Author: Bailey Donaldson
;
;Purpose: objective 3 on the REU packet, compare the calculated kappana and the given kappa by
;   graphing kappa vs. Sorce day and kappa vs. wavelength.
;
;Notes: you need to have run the programas "convolve_nrlssi.pro" and "calculate_kappa"
;   before one can use this program

function plot_kappa, Inst_mode_ID


  if Inst_mode_ID eq 41 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/VIS_data.sav'
    data = vis_data  
    restore, '/Users/bado8599/Desktop/REU/sav_files/VIS_daily_kappa.sav'
    type = 'VIS'
  endif
  
  if Inst_mode_ID eq 43 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/UV_data.sav'
    data = uv_data 
    restore, '/Users/bado8599/Desktop/REU/sav_files/UV_daily_kappa.sav'
    type = 'UV'
  endif
  
  if Inst_mode_ID eq 44 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/IR_data.sav'
    data = ir_data
    restore, '/Users/bado8599/Desktop/REU/sav_files/IR_daily_kappa.sav'
    type = 'IR'
  endif
 
 ;provided kappa by LASP, average kappa, varies over wavelength but not Sorce Day
  readcol, '/Users/bado8599/Desktop/REU/Bailey/ESR_kappa_453_1570_2011_mod_407nm.txt', $
    kappa_w, kappa_k, format = '(d,d)', /silent

 
  choice = ''
  READ, choice, PROMPT = 'Plot prism degredation with constant WL or SD? '

  ;this section plots Sorce Day vs. kappa
  if strlowcase(choice) eq 'wl' then begin
    READ, chosen_wl, PROMPT = 'At which constant wavelength would you like to plot? '
    
    ;given kappa only has one kappa value at each wavelength so finding which kappa value is closest to the chosen_wl
    ;only needs to be done once (so it's out of the for loop)
    nunca = min(abs(kappa_w - chosen_wl), pos_wl_kappa)
    
    SIM_timestamp_array = []
    calculated_kappa_array = []
    given_kappa_array = []
    for i=0, n_elements(data)-1 do begin
            
      ;finds subscript (position) of wavelength in every measurement
      nanimo = min(abs((*data[i]).wavelength - chosen_wl), pos_wl)
      calculated_kappa_array = [calculated_kappa_array, (*daily_kappa[i]).kappa[pos_wl]]
      
      SIM_timestamp_array = [SIM_timestamp_array, (*data[i]).SIM_timestamp]
      given_kappa_array = [given_kappa_array, kappa_k[pos_wl_kappa]]
      
    endfor
    
    
    ;the following code creates arrays with the outlying data eliminated
    mean_kappa = total(calculated_kappa_array)/n_elements(calculated_kappa_array)
    if mean_kappa lt 0.0 then mean_kappa = -1*mean_kappa
    
    calculated_kappa_no_outliers = []
    SIM_timestamp_no_outliers = []
    given_kappa_no_outliers = []
    for j=0, n_elements(calculated_kappa_array)-1 do begin
      ;removes kappa that are a random isolated outlier
      if calculated_kappa_array[j] ge calculated_kappa_array[j-1]-.001 and calculated_kappa_array[j] le calculated_kappa_array[j-1]+0.001 then begin
        ;removes kappa that is far outside the mean, even if there are consecutive kappa outliers
        if calculated_kappa_array[j] ge mean_kappa-(mean_kappa*5) and calculated_kappa_array[j] le mean_kappa+(mean_kappa*5) then begin
          calculated_kappa_no_outliers = [calculated_kappa_no_outliers, calculated_kappa_array[j]]
          SIM_timestamp_no_outliers = [SIM_timestamp_no_outliers, SIM_timestamp_array[j]]
          given_kappa_no_outliers = [given_kappa_no_outliers, given_kappa_array[j]]
        endif
      endif
        
    endfor
    
    ;turns user_wl into a string so it can be used in the title
    chosen_wl_string = strtrim(chosen_wl, 1)
    
    device, decompose = 0
   
    lineplot, SIM_timestamp_array, calculated_kappa_array, ptitle = 'Given Kappa and Calculated Kappa comparison by Sorce Day at '$
       + chosen_wl_string +'nm', xtitle = 'Sorce Day', ytitle = 'Kappa', title = 'Calculated Daily Kappa'
       
    ;should result in a plot of a horizontal line since given kappa is constant across Sorce Day
    lineplot, SIM_timestamp_array, given_kappa_array, title = 'Given average kappa'

    ;plot without outliers of calculated daily kappa
    lineplot, SIM_timestamp_no_outliers, calculated_kappa_no_outliers, title = 'Calculated Daily Kappa without Outliers'
    print, 'mean', mean_kappa
  endif
  
 
  
  if strlowcase(choice) eq 'sd' then begin  
    READ, chosen_sd, PROMPT = 'Choose a Mission Timestamp (Sorce Day): '

    ;takes SIM timestamps out of the pointer and puts it into an array
    SIM_timestamp_array = []
    for k=0, n_elements(data)-1 do begin
      SIM_timestamp_array = [SIM_timestamp_array, (*data[k]).SIM_timestamp]
    endfor

    ;find timestamp closest to the timestamp the user chose
    nanimo = min(abs(SIM_timestamp_array - chosen_sd), pos_SD)
    given_kappa_array = interpol(kappa_k, kappa_w, (*data[pos_SD]).wavelength)
    
    
    mean_kappa = total((*daily_kappa[pos_SD]).kappa)/n_elements((*daily_kappa[pos_SD]).kappa)
    if mean_kappa lt 0.0 then mean_kappa= (-1)*mean_kappa
    daily_kappa_no_outliers = []
    wavelength_no_outliers = []
    for l=0, n_elements((*data[pos_SD]).wavelength)-1 do begin
      ;these coefficients may need to be adjusted depending on the graph
      if (*daily_kappa[pos_SD]).kappa[l] ge mean_kappa*(-5.0) and (*daily_kappa[pos_SD]).kappa[l] le mean_kappa*5 then begin
        daily_kappa_no_outliers = [daily_kappa_no_outliers, (*daily_kappa[pos_SD]).kappa[l]]
        wavelength_no_outliers = [wavelength_no_outliers, (*data[pos_SD]).wavelength[l]]
      endif  
    endfor
    
    ;turns SD into a string so that it can be used in the title
    SD_string = strtrim((*data[pos_SD]).SIM_timestamp, 1)

    device, decompose = 0   
    
    lineplot, (*data[pos_SD]).wavelength, (*daily_kappa[pos_SD]).kappa, title = 'Wavelength vs. Calculated Daily Kappa'
    
    ;given average kappa
    lineplot, (*data[pos_SD]).wavelength, given_kappa_array, title = 'Wavelength vs. Given Average Kappa', $
      xtitle = 'Wavelength (nm)', ytitle = 'Kappa', $
      ptitle = 'Given Average Kappa and Calculated Daily Kappa Comparison in the ' + type + ' at SD ' + SD_string   
     
    lineplot, wavelength_no_outliers, daily_kappa_no_outliers, title = 'Wavelength vs. Calculated Daily Kappa no Outliers'
     
  endif

  
end