;+
;Author: Laura Sandoval, Bailey Donaldson
;
;
;STEPS:
;(1) reads in the NRLSSI data and the SIM data for the chosen instrument mode ID
;(2) adjusts the NRLSSI data to the correct SD units
;(3) finds the kernel for the convolution
;(3) extracts the data for the SD from NRLSSI data
;(4) the NRLSSI data and SIM data are interpolated and so they are on the same wvl grid 
;     the NRLSSI data is also convoluted
;
;
;PURPOSE: to interpolate the SIM and NRLSSI data onto the same wavelength grid, to convolve
;the NRLSSI data to better match SIM data, and to save the new data into a sav file
;
;*******************************************************

function convolve_nrlssi, Inst_mode_ID 

;***************STEP 1****************

;gets data from NRLSSI, saves it into structure called data
;  read_netCDF, '/Users/bado8599/Desktop/REU/ssi_v02r00_daily_s18820101_e20161231_c20170328.nc', data
  restore, '/Users/bado8599/Desktop/REU/sav_files/original_NRLSSI_data.sav'

;restores the file to get the wavelength grid
  restore, '/Users/bado8599/Desktop/REU/sav_files/SimStandardWavelength.sav'

;get the data from the SIM uncorrected irradiance,
  if Inst_mode_ID eq 41 then begin
    filename='/Users/bado8599/Desktop/REU/uncorr_irr/visa_irrad_23_uncorr.sav'
    titlename='VIS'
    wl_grid = Standard_wave_41.stdwavelength
    ;number of elements the wavelength grid must contain for the data to be added to the final pointer at the end (new_data)
    ;this number was found by finding the mean of the amount of elements, then finding a 5% deviation and subtracting from the mean
    ;in VIS, wl_range must have at least 536 elements to have 2 standard deviation    
    elements= 536.0
    endif

  if Inst_mode_ID eq 43 then begin
    filename='/Users/bado8599/Desktop/REU/uncorr_irr/uva_irrad_23_uncorr.sav'
    titlename='UV'
    wl_grid = Standard_wave_43.stdwavelength
    elements = 1098
    endif

  if Inst_mode_ID eq 44 then begin
    filename= '/Users/bado8599/Desktop/REU/uncorr_irr/ira_irrad_23_uncorr.sav'
    titlename='IR'
    wl_grid = Standard_wave_44.stdwavelength
    elements = 141.0
    endif

;restores the SIM file and renames the variable it contains to the name irraidance 
  restore, filename
  title_name = FILE_BASENAME(filename, '.sav')
  assign_irr = 'irradiance = ' + title_name
  void = execute(assign_irr)


;****************STEP 2****************
;adjusts the NRLSSI data.time to be in SIM sorce days, renames it nrlssi_ts_array
  nrlssi_day0 = ymd2sd(1610, 01, 01)
  data.time += nrlssi_day0
  postime = where(data.time ge 0.0)
  nrlssi_ts_array = data.time[postime]
  nrlssi_irrad = data.ssi[*,postime]


;****************STEP 3****************
; We know the kernel needs to have 13 elements (to do a weighed average
; of 13 points around each point) and the maximum is in the middle at position 6 (triangular profile)
; and has a normalized height of 1.0:
;
;              /\
;             /  \
;            /    \
;           /      \
;          /        \
;    _____/          \_______
;
  val= 1d/6d
  kernel_y = [0d, val, val*2, val*3, val*4, val*5, val*6, val*5, val*4, val*3, val*2, val, 0d]
  kernel_y /= total(kernel_y)
  
  
;****************STEP 4****************
  ;this loop finds out how long the pointer new_data needs to be
  length = 0.0
  for j=0,n_elements(irradiance.spect20)-1 do begin
    pos_sim=where(wl_grid lt (*irradiance.spect20[j]).wavelength[-1] and wl_grid gt (*irradiance.spect20[j]).wavelength[0], count1)
    if count1 lt 1.0 then continue
    if n_elements(pos_sim) ge elements then length++
  endfor

  nrlssi_interpolated = []
  sim_interpolated = []
  new_data = ptrarr(length)
  count = 0.0
  for j=0,n_elements(irradiance.spect20)-1 do begin
    
    ;timestamp is the Sorce Day of SIM, it will go through each SD as the loop appends
    sim_time = gps2sd((*irradiance.spect20[j]).timestamp[0]/1d6)

    ;finds the index of the SD within NRLSSI data that is closest to the timestamp
    nunca = min(abs(nrlssi_ts_array - sim_time), pos_nrlssi_ts)
    
    ;results in array of NRLSSI and SIM interpolated on the same wavelength grid  
    pos_sim=where(wl_grid lt (*irradiance.spect20[j]).wavelength[-1] and wl_grid gt (*irradiance.spect20[j]).wavelength[0], count1)  
    if count1 lt 1.0 then continue
    
    nrlssi_temp = reform(nrlssi_irrad[*, pos_nrlssi_ts])
    nrlssi_interpolated = interpol(nrlssi_temp, data.wavelength, wl_grid)

    ;convolve the interpolated NRLSSI2 irradiance array using the kernel found
    convolved_NRLSSI = convol(nrlssi_interpolated, kernel_y,/edge_mirror)

     if n_elements(pos_sim) ge elements then begin  
       sim_interpolated = interpol((*irradiance.spect20[j]).irradiance, (*irradiance.spect20[j]).wavelength, wl_grid[pos_sim])
       
       new_data[count]=ptr_new({wavelength: wl_grid[pos_sim], irradiance_NRLSSI_conv:convolved_NRLSSI[pos_sim], irradiance_NRLSSI_int:nrlssi_interpolated[pos_sim], $
        irradiance_SIM_int:sim_interpolated, SIM_timestamp:sim_time, NRLSSI_timestamp:nrlssi_ts_array[pos_nrlssi_ts], f_function:dblarr(n_elements(pos_sim))} )
       count++
    endif

  endfor
  
  if Inst_mode_ID eq 41 then begin
    VIS_data = new_data
    save, VIS_data, filename = '/Users/bado8599/Desktop/REU/sav_files/VIS_data.sav'
  endif
  if Inst_mode_ID eq 43 then begin
    UV_data = new_data
    save, UV_data, filename = '/Users/bado8599/Desktop/REU/sav_files/UV_data.sav'
  endif
  if Inst_mode_ID eq 44 then begin
    IR_data = new_data
    save, IR_data, filename = '/Users/bado8599/Desktop/REU/sav_files/IR_data.sav'
  endif
  
  return, 0
end

