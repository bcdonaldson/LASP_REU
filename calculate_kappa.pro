;
;Author: Bailey Donaldson
;
;Purpose: caluclate kappa that minimizes this equation: 0 = PRISM_DEG - SIM_UNCORR_IRR/NRLSSI2_IRR
;   where f' is kept constant as 1.00. This code is for SIMA and can be used for the IR, Vis, and UV.
;   
;NOTE: the final measurement contained in the solar exposure file if for SD 4864.98
;

function amoeba_kappa, kappa

  common amoeba_kappa_common, data, solar54, attenuation, j, i, pos_time
    
    zero = abs((1d - attenuation[i])*exp(-kappa*solar54[pos_time].solar_exp) +$
      attenuation[i]*exp(-kappa*solar54[pos_time].solar_exp*(1d/2d)) - $
      (*data[j]).irradiance_sim_int[i]/(*data[j]).irradiance_nrlssi_conv[i])

  return, zero

end

function calculate_kappa, Inst_mode_ID
  common amoeba_kappa_common
  
  ;must run code 'convolve_nrlssi.pro' to obtain the following 3 files
  ;these three if statements will restore the SIM data, NRLSSI2 data, wavelength grid, and timestamp
  if Inst_mode_ID eq 41 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/VIS_data.sav'
    data = vis_data
  endif
  
  if Inst_mode_ID eq 43 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/UV_data.sav'
    data = uv_data
  endif
  
  if Inst_mode_ID eq 44 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/IR_data.sav'
    data = ir_data
  endif
  
  ;solar exposure for SIMA
  restore, '/Users/bado8599/Desktop/REU/sav_files/solarexp_54_55_raw.sav'
  solar54.solar_exp /= 86400d
  
  ;a value, these values are measured at different wavelengths then the SIM & NRLSSI data so it will be interpolated
  readcol, '/Users/bado8599/Desktop/REU/sav_files/overlap_vis1_mod_407nm.txt', ray_w, ray_a, format = '(d,d)', /silent
  ;attenuation = {wavelength: ray_w, a: ray_a}
  
  ;kappa (will use this to approximate the kappa values in amoeba)
  readcol, '/Users/bado8599/Desktop/REU/Bailey/ESR_kappa_453_1570_2011_mod_407nm.txt', kappa_w, kappa_k, format = '(d,d)', /silent
  ;kappa = {wavelength: kappa_w, k:kappa_k}
  
  ;adjust values of T1 to SD from GPS time in micro seconds
  solar_exp_SD = gps2sd(solar54.T1/1d6)
  
  ;creates an empty pointer that the final kappa values will be saved into
  daily_kappa = ptrarr(n_elements(data))
  
  ;j loops through each measurement taken by SIM (a couple each Sorce Day)
  for j=0,n_elements(data)-1 do begin
    
    ;interpolated within the loop because the amount of elements in *data changes each measurement
    attenuation = interpol(ray_a, ray_w, (*data[j]).wavelength)
    initial_kappa = interpol(kappa_k, kappa_w, (*data[j]).wavelength)
    
    ;finds where the timestamp in solar exposure is closest to SIM timestamp
    hello = min(abs(solar_exp_SD - (*data[j]).sim_timestamp[0]), pos_time)
    
    ;this is initialized within the first for loop so that at every Sorce Day, a new calculated_kappa is created  
    kappa_array = []
    
    for i=0, n_elements((*data[j]).wavelength)-1 do begin
      ;initial kappa (P0) is approximated by using the average kappa values 
      ;contains kappa at one measurement (each value of j is one measurement) for every wavelength
      kappa = amoeba(1d-6, P0=initial_kappa[i], SCALE=1d, FUNCTION_VALUE=fval, FUNCTION_NAME='amoeba_kappa', ncalls=ncalls, nmax=maxiter)
      
      kappa_array = [kappa_array, kappa]
      
      daily_kappa[j] = ptr_new({kappa: kappa_array})  
      
    endfor
    
  endfor

  
  if Inst_mode_ID eq 41 then begin
    save, daily_kappa, filename = '/Users/bado8599/Desktop/REU/sav_files/VIS_daily_kappa.sav'
  endif

  if Inst_mode_ID eq 43 then begin
    save, daily_kappa, filename = '/Users/bado8599/Desktop/REU/sav_files/UV_daily_kappa.sav'
  endif

  if Inst_mode_ID eq 44 then begin
    save, daily_kappa, filename = '/Users/bado8599/Desktop/REU/sav_files/IR_daily_kappa.sav'
  endif 
  
end