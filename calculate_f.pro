;
;Author: Bailey Donaldson
;
;
;Purpose: find the f_function that minimizes this equation: 0 = PRISM_DEG - SIM_UNCORR_IRR/NRLSSI2_IRR
;This is for SIMA. This can be used for VIS, IR, and UV.
;
;Example: calculate_f(41)
;
;NOTE: the final measurement contained in the solar exposure file if for SD 4864.98, so measurements after this date are less accurate
;
function amoeba_f_function, f_function   

  common amoeba_f_function_common, i, j, pos_time, attenuation1, kappa1, data, solar54

      zero = abs((1d -attenuation1[i])*exp(-kappa1[i]*solar54[pos_time].solar_exp*(f_function)) +$
        attenuation1[i]*exp(-kappa1[i]*solar54[pos_time].solar_exp*(f_function)/2d) -$
        ((*data[j]).irradiance_sim_int[i])/((*data[j]).irradiance_nrlssi_conv[i]))   
  
  return, zero 
  
end


function calculate_f, Inst_mode_ID
  common amoeba_f_function_common
    
    
  if Inst_mode_ID eq 41 then begin
    ;SIM data, NRLSSI data, and wavelength grid in VIS light spectrum
    ;can get this file (and ones for UV and IR) by using the program 'convolve_nrlssi.pro'
    restore, '/Users/bado8599/Desktop/REU/sav_files/VIS_data.sav'
    data = vis_data
  endif
  
  if Inst_mode_ID eq 43 then begin
    ;SIM data, NRLSSI data, and wavelength grid in VIS light spectrum
    restore, '/Users/bado8599/Desktop/REU/sav_files/UV_data.sav'
    data = uv_data
  endif
  
  if Inst_mode_ID eq 44 then begin
    ;SIM data, NRLSSI data, and wavelength grid in VIS light spectrum   
    restore, '/Users/bado8599/Desktop/REU/sav_files/ir_data.sav' 
    data = ir_data   
  endif
  
  ;solar exposure for SIMA
  restore, '/Users/bado8599/Desktop/REU/sav_files/solarexp_54_55_raw.sav'
  solar54.solar_exp /= 86400d

  ;a value
  readcol, '/Users/bado8599/Desktop/REU/sav_files/overlap_vis1_mod_407nm.txt', ray_w, ray_a, format = '(d,d)', /silent
  attenuation = {wavelength: ray_w, a: ray_a}
  
  ;kappa value
  readcol, '/Users/bado8599/Desktop/REU/Bailey/ESR_kappa_453_1570_2011_mod_407nm.txt', kappa_w, kappa_k, format = '(d,d)', /silent
  kappa = {wavelength: kappa_w, k:kappa_k}
  
  ;adjust values of T1 to SD from GPS time in micro seconds   
  solar_exp_SD = gps2sd(solar54.T1/1d6)
      

  ;j loops through each measurement taken by SIM (a couple each Sorce Day)
  for j=0, n_elements(data)-1 do begin
 
    ;is interpolated within the loop because the amount of elements in *data changes every measurement 
    attenuation1 = interpol(attenuation.a, attenuation.wavelength, (*data[j]).wavelength)   
    kappa1 = interpol(kappa.k, kappa.wavelength, (*data[j]).wavelength)
    
    ;finds where timestamp in solar_exp is closest to SIM timestamp
    hello = min(abs(solar_exp_SD - (*data[j]).sim_timestamp[0]), pos_time)

    for i=0, n_elements((*data[j]).wavelength)-1 do begin
          
      ;puts (*data[j]).f_function[i] into a format that can be passed through the first function
      f_function = (*data[j]).f_function[i]
      
      ;finds the value of f_function that would minimize the result of the equation created in my_amoeba_fdeg
      f_function = amoeba(1d-6, P0=1d, SCALE=1d, FUNCTION_VALUE=fval, FUNCTION_NAME='amoeba_f_function', ncalls = ncalls)  
        
      ;fills the empty structure with the result of the minimization
      (*data[j]).f_function[i] = f_function
     
    endfor
  endfor
  
  if Inst_mode_ID eq 41 then begin
    VIS_data = data
    save, VIS_data, filename = '/Users/bado8599/Desktop/REU/sav_files/VIS_data_F_function.sav'
  endif

  if Inst_mode_ID eq 43 then begin
    UV_data = data
    save, UV_data, filename = '/Users/bado8599/Desktop/REU/sav_files/UV_data_F_function.sav'
  endif

  if Inst_mode_ID eq 44 then begin
    IR_data = data
    save, IR_data, filename = '/Users/bado8599/Desktop/REU/sav_files/IR_data_F_function.sav'
  endif  
  
  return, 0
end
