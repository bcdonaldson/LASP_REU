
;Author: Bailey Donaldson

;Purpose: using the prism degradation found in objective 2 (was calculated using the f function value that was found 
      ;in the code calculate_f) and the uncorrected irradiance to see if the corrected irradiance equals the nrlssi irradiance. 

function calculate_corrected_irr, Inst_mode_ID

  if Inst_mode_ID eq 41 then begin
    ;contains SIM uncorrected irradiance interpolated
    restore, '/Users/bado8599/Desktop/REU/sav_files/VIS_data.sav'
    data = vis_data
    ;contains prism degradation on the same wavelength grid as SIM interpolated 
    restore, '/Users/bado8599/Desktop/REU/sav_files/prism_deg_VIS.sav'
    
  endif

  if Inst_mode_ID eq 43 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/UV_data.sav'
    data = uv_data
    restore, '/Users/bado8599/Desktop/REU/sav_files/prism_deg_UV.sav'
  endif

  if Inst_mode_ID eq 44 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/IR_data.sav'
    data = ir_data
    restore, '/Users/bado8599/Desktop/REU/sav_files/prism_deg_IR.sav'
  endif

  SIM_corrected_irr = ptrarr(n_elements(data))
  ;(SIM_uncorr_irr)/(prism_deg) = SIM_corrected_irr
  for i=0, n_elements(data)-1 do begin
    
    calculated_corrected_irr_avg_kappa = []   
    calculated_corrected_irr_daily_kappa = [] 
    
    for j=0, n_elements((*data[i]).wavelength)-1 do begin
      
      calculated_corrected_irr_avg_kappa = [calculated_corrected_irr_avg_kappa, ((*data[i]).irradiance_sim_int[j])/((*prism_deg[i]).prism_deg_avg_kappa[j])]
      calculated_corrected_irr_daily_kappa = [calculated_corrected_irr_daily_kappa, ((*data[i]).irradiance_sim_int[j])/((*prism_deg[i]).prism_deg_daily_kappa[j])]
      
      SIM_corrected_irr[i] = ptr_new({irrad_avg_kappa: calculated_corrected_irr_avg_kappa, irrad_daily_kappa: calculated_corrected_irr_daily_kappa})
      
    endfor
  endfor
  if Inst_mode_ID eq 41 then begin
    save, SIM_corrected_irr, filename = '/Users/bado8599/Desktop/REU/sav_files/SIM_corrected_irr_VIS.sav'
  endif
  if Inst_mode_ID eq 43 then begin
    save, SIM_corrected_irr, filename = '/Users/bado8599/Desktop/REU/sav_files/SIM_corrected_irr_UV.sav'
  endif
  if Inst_mode_ID eq 44 then begin
    save, SIM_corrected_irr, filename = '/Users/bado8599/Desktop/REU/sav_files/SIM_corrected_irr_IR.sav'
  endif
 
end

;***************************************************************************************************************************************

function plot_corrected_irr, Inst_mode_ID

  if Inst_mode_ID eq 41 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/SIM_corrected_irr_VIS.sav'
    title = 'VIS'
    restore, '/Users/bado8599/Desktop/REU/sav_files/VIS_data.sav'
    data = vis_data
  endif
  if Inst_mode_ID eq 43 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/SIM_corrected_irr_UV.sav'
    title = 'UV'
    restore, '/Users/bado8599/Desktop/REU/sav_files/UV_data.sav'
    data = uv_data
  endif
  if Inst_mode_ID eq 44 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/SIM_corrected_irr_IR.sav'
    title = 'IR'
    restore, '/Users/bado8599/Desktop/REU/sav_files/IR_data.sav'
    data = ir_data
  endif
  
  
  READ, user_wl, PROMPT = 'Choose a Wavelength: '
  wl = 0
  SIM_timestamp = []
  sim_corr_irrad_avg_kappa = []
  sim_corr_irrad_daily_kappa = []
  nrlssi_irrad_conv = []
  nrlssi_irrad_int = []
  sim_uncorr_irr = []
  for k=0, n_elements(data)-1 do begin
    ;puts the timestamp into a format we can use to graph (takes it out of the pointer)
    SIM_timestamp = [SIM_timestamp, (*data[k]).SIM_timestamp]

    ;finds the data that has the position in the array is closest to the wavelength the user chose
    ;it is looped because the wavelength's position will be different for each Sorce Day
    hello = min(abs((*data[k]).wavelength - user_wl), pos_wl)
    
    ;creates an array of the irradiance at each SD at a constant wavelength
    sim_corr_irrad_avg_kappa = [sim_corr_irrad_avg_kappa, (*sim_corrected_irr[k]).irrad_avg_kappa[pos_wl]]
    sim_corr_irrad_daily_kappa = [sim_corr_irrad_daily_kappa, (*sim_corrected_irr[k]).irrad_daily_kappa[pos_wl]]
    nrlssi_irrad_conv = [nrlssi_irrad_conv, (*data[k]).irradiance_nrlssi_conv[pos_wl]]
    nrlssi_irrad_int = [nrlssi_irrad_int, (*data[k]).irradiance_nrlssi_int[pos_wl]] 
    sim_uncorr_irr = [sim_uncorr_irr, (*data[k]).irradiance_sim_int[pos_wl]]
  endfor

  no_outliers_avg_kappa = []
  no_outliers_daily_kappa = []
  SIM_timestamp_no_outliers = []
  SIM_timestamp_no_outliers1 = []
  ;try to remove the data outliers
  corr_irr_mean = total(sim_corr_irrad_avg_kappa)/n_elements(sim_corr_irrad_avg_kappa)
  corr_irrad_median_avg_kapa = median(sim_corr_irrad_avg_kappa, /double, /even)
  corr_irrad_median_daily_kappa = median(sim_corr_irrad_daily_kappa, /double, /even)
  
  ;the numbers below may need to be adjusted depending on the wavelength graphed
  for m=0, n_elements(sim_corr_irrad_avg_kappa)-1 do begin
    ;the distance from the mean that is permitted may need to be altered
    if sim_corr_irrad_avg_kappa[m] ge corr_irrad_median_avg_kapa-(corr_irrad_median_avg_kapa*.5) and sim_corr_irrad_avg_kappa[m] le corr_irrad_median_avg_kapa+(corr_irrad_median_avg_kapa) then begin
      ;the difference permitted between each number may have to be altered
      ;if sim_corr_irrad_avg_kappa[m] ge sim_corr_irrad_avg_kappa[m-1]-1 and sim_corr_irrad_avg_kappa[m] le sim_corr_irrad_avg_kappa[m-1]+1 then begin
        no_outliers_avg_kappa = [no_outliers_avg_kappa, sim_corr_irrad_avg_kappa[m]]
        SIM_timestamp_no_outliers = [SIM_timestamp_no_outliers, SIM_timestamp[m]]
      ;endif
    endif
      
    if sim_corr_irrad_daily_kappa[m] ge corr_irrad_median_daily_kappa-(corr_irrad_median_daily_kappa*.5) and sim_corr_irrad_daily_kappa[m] le corr_irrad_median_daily_kappa+(corr_irrad_median_daily_kappa) then begin
      ;if sim_corr_irrad_daily_kappa[m] ge sim_corr_irrad_daily_kappa[m-1]-1 and sim_corr_irrad_daily_kappa[m] le sim_corr_irrad_daily_kappa[m-1]+1 then begin
        no_outliers_daily_kappa = [no_outliers_daily_kappa, sim_corr_irrad_daily_kappa[m]]
        SIM_timestamp_no_outliers1 = [SIM_timestamp_no_outliers1, SIM_timestamp[m]]
      ;endif
    endif

  endfor
 
  ;turns user_wl into a string so it can be used in the title
  hello = min(abs((*data[100]).wavelength - user_wl), pos_wl)

  wl_string = strtrim((*data[100]).wavelength[pos_wl], 1)

  ;keeps graph from turing red
  device, decompose = 0

  lineplot, SIM_timestamp, nrlssi_irrad_conv, title = 'Convolved NRLSSI2 Irradiance'

  lineplot, SIM_timestamp, nrlssi_irrad_int, title = 'Interpolated NRLSSI2 Irradiance at Wavelength: ' + wl_string + 'nm'

  lineplot, SIM_timestamp, sim_uncorr_irr, title = 'SIM uncorrected Irradiance', ytitle = 'Irradiance', xtitle = 'Sorce Day',$
    ptitle = 'Irradiance vs. Sorce Day for SIMA at 280.05nm'

  ;lineplot, SIM_timestamp, sim_corr_irrad_avg_kappa, xtitle = 'Sorce Day', ytitle = 'Irradiance'

  ;lineplot, SIM_timestamp_no_outliers, no_outliers_avg_kappa, title = 'SIM Corrected Irradiance (f function) without outliers'
  
  ;lineplot, SIM_timestamp, sim_corr_irrad_daily_kappa, title = 'SIM corrected Irradiance (daily kappa)'


  ;lineplot, SIM_timestamp_no_outliers1, no_outliers_daily_kappa, title = 'SIM corrected Irradiance (daily kappa) without outliers'



end

;***************************************************************************************************************************************

function plot_by_SD, Inst_mode_ID
  if Inst_mode_ID eq 41 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/SIM_corrected_irr_VIS.sav'
    title = 'VIS'
    restore, '/Users/bado8599/Desktop/REU/sav_files/VIS_data.sav'
    data = vis_data
  endif
  if Inst_mode_ID eq 43 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/SIM_corrected_irr_UV.sav'
    title = 'UV'
    restore, '/Users/bado8599/Desktop/REU/sav_files/UV_data.sav'
    data = uv_data
  endif
  if Inst_mode_ID eq 44 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/SIM_corrected_irr_IR.sav'
    title = 'IR'
    restore, '/Users/bado8599/Desktop/REU/sav_files/IR_data.sav'
    data = ir_data
  endif
  
  READ, chosen_sd, PROMPT = 'Choose a Mission Timestamp (Sorce Day): '

  ;takes SIM timestamps out of the pointer and puts it into an array
  SIM_timestamp_array = []
  for k=0, n_elements(data)-1 do begin
    SIM_timestamp_array = [SIM_timestamp_array, (*data[k]).SIM_timestamp]
  endfor

  ;find timestamp closest to the timestamp the user chose
  nanimo = min(abs(SIM_timestamp_array - chosen_sd), pos_SD)
  ;given_kappa_array = interpol(kappa_k, kappa_w, (*data[pos_SD]).wavelength)

  sd_string = strtrim(SIM_timestamp_array[pos_SD],1)
  device, decompose = 0
  
  lineplot, (*data[pos_SD]).wavelength, (*data[pos_SD]).irradiance_nrlssi_int, title = 'NRLSSI Interpolated',$
    ptitle='NRLSSI & SIM comparison at SD: '+sd_string, xtitle='Wavelength (nm)', ytitle='Irradiance'
  lineplot, (*data[pos_SD]).wavelength, (*data[pos_SD]).irradiance_nrlssi_conv, title = 'NRLSSI Interpolated and Convolved'
  lineplot, (*data[pos_SD]).wavelength, (*data[pos_SD]).irradiance_sim_int, title = "SIM Irradiance Uncorrected'

end

;********************************************************************************************************************************

;plots (nrlssi_irrad - sim_corr_irrad)/(nrlssi_irrad), aka the absolute error
function plot_irrad_error, Inst_mode_ID

  if Inst_mode_ID eq 41 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/SIM_corrected_irr_VIS.sav'
    title = 'VIS'
    restore, '/Users/bado8599/Desktop/REU/sav_files/VIS_data.sav'
    data = vis_data
  endif
  if Inst_mode_ID eq 43 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/SIM_corrected_irr_UV.sav'
    title = 'UV'
    restore, '/Users/bado8599/Desktop/REU/sav_files/UV_data.sav'
    data = uv_data
  endif
  if Inst_mode_ID eq 44 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/SIM_corrected_irr_IR.sav'
    title = 'IR'
    restore, '/Users/bado8599/Desktop/REU/sav_files/IR_data.sav'
    data = ir_data
  endif
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  READ, user_wl, PROMPT = 'Choose a Wavelength: '
  wl = 0
  SIM_timestamp = []
  error_daily_kappa = []
  error_ffunction = []
  
  for i=0, n_elements(data)-1 do begin
    
    SIM_timestamp = [SIM_timestamp, (*data[i]).SIM_timestamp]
    hello = min(abs((*data[i]).wavelength - user_wl), pos_wl)
    error_daily_kappa = [error_daily_kappa, ((*data[i]).irradiance_nrlssi_conv[pos_wl] - (*SIM_corrected_irr[i]).irrad_daily_kappa[pos_wl])/(*data[i]).irradiance_nrlssi_conv[pos_wl]]
    error_ffunction = [error_ffunction, ((*data[i]).irradiance_nrlssi_conv[pos_wl] - (*SIM_corrected_irr[i]).irrad_avg_kappa[pos_wl])/((*data[i]).irradiance_nrlssi_conv[pos_wl])]
  
  endfor
  k1=where(error_daily_kappa lt .9)
  k2=where(error_ffunction ne 1)
  
  wl_string = strtrim((*data[n_elements(data)-1]).wavelength[pos_wl], 1)
  device, decompose = 0
    
  lineplot, SIM_timestamp[k1], error_daily_kappa[k1], title = 'Error of SIM Irradiance Calculated with Daily Kappa', xtitle = 'Sorce Day', ytitle = 'Irradiance'
  
  lineplot, SIM_timestamp[k2], error_ffunction[k2], title = 'Error of SIM Irradiance Calculated with F Function', ptitle = 'Absolute Error of SIM corrected irradiance at 280.05nm'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  READ, user_sd, PROMPT = 'Choose a SD: '
;  
;  hello1 = min(abs(SIM_timestamp - user_sd), pos_sd)
;  error_kappa2 = []
;  error_ffunction2 = []
;  for j=0, n_elements((*data[pos_sd]).wavelength)-1 do begin
;    error_kappa2 = [error_kappa2, ((*data[pos_sd]).irradiance_nrlssi_conv[j] - (*SIM_corrected_irr[pos_sd]).irrad_daily_kappa[j])/(*data[pos_sd]).irradiance_nrlssi_conv[j]]
;    error_ffunction2 = [error_ffunction2, ((*data[pos_sd]).irradiance_nrlssi_conv[j] - (*SIM_corrected_irr[pos_sd]).irrad_avg_kappa[j])/(*data[pos_sd]).irradiance_nrlssi_conv[j]]
;  endfor
;  
;  k=where(error_kappa2 ne 1)
;  k1=where(error_ffunction2 ne 1)
;  sd_string = strtrim(sim_timestamp[pos_sd],1)
;  
;  lineplot, (*data[pos_sd]).wavelength[k], error_kappa2[k], title = 'SIM Irradiance Calculated with Daily Kappa at SD '+sd_string,$
;     ptitle = 'Absolute Error of SIM corrected irradiance vs. Wavelength in the ' + title
;  
;  lineplot, (*data[pos_sd]).wavelength[k1], error_ffunction2[k1], title = 'SIM Irradiance Calculated with F Function at SD '+sd_string, xtitle = 'Wavelength', ytitle = 'Irradiance'
;  
end