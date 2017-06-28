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
;(5) create pointer to save new (interpolated and convoluted) data in
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
    WL_range = where(wl_grid le 950.0 and wl_grid ge 310.0)
    ;number of elements the wavelength grid must contain for the data to be added to the final pointer at the end (new_data)
    ;this number was found by finding the mean of the amount of elements, then finding a 5% deviation and subtracting from the mean
    ;in VIS, wl_range must have at least 536 elements to have 2 standard deviation    
    elements= 536.0
    endif

  if Inst_mode_ID eq 43 then begin
    filename='/Users/bado8599/Desktop/REU/uncorr_irr/uva_irrad_23_uncorr.sav'
    titlename='UV'
    wl_grid = Standard_wave_43.stdwavelength
    WL_range = where(wl_grid le 306.0 and wl_grid ge 240.0)
    ;in UV, wl_range must have at least 1098 elements to have 2SD
    elements = 1098
    endif

  if Inst_mode_ID eq 44 then begin
    filename= '/Users/bado8599/Desktop/REU/uncorr_irr/ira_irrad_23_uncorr.sav'
    titlename='IR'
    wl_grid = Standard_wave_44.stdwavelength
    WL_range = where(wl_grid le 1600.0 and wl_grid ge 950.0)
    ;in IR, wl_range must have at least 141 elements to have 2SD
    elements = 141.0
    endif

;restores the SIM file and renames the variable it contains to the name irraidance 
  restore, filename
  title_name = FILE_BASENAME(filename, '.sav')
  assign_irr = 'irradiance = ' + title_name
  void = execute(assign_irr)


;****************STEP 2****************
;
;adjusts the NRLSSI data.time to be in SIM sorce days, renames it NRLSSI_TS_array
  nrlssi_day0 = ymd2sd(1610, 01, 01)
  data.time += nrlssi_day0
  postime = where(data.time ge 0.0)
  NRLSSI_TS_array = data.time[postime]

;VIS

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
;
  val= 1d/6d
  kernel_y = [0d, val, val*2, val*3, val*4, val*5, val*6, val*5, val*4, val*3, val*2, val, 0d]

  ;divide kernel_y by the sum of its values to create the final kernel_y
  kernel_y /= total(kernel_y)

;this loop finds out how long the pointer new_data needs to be
  length = 0.0
  for j=0,n_elements(irradiance.spect20)-1 do begin
    pos_sim=where(wl_grid lt (*irradiance.spect20[j]).wavelength[-1] and wl_grid gt (*irradiance.spect20[j]).wavelength[0], count1)
    if count1 lt 1.0 then continue
    if n_elements(pos_sim) ge elements then length++
  endfor

;****************STEP 4****************
  NRLSSI_SSI_new = []
  SIM_SSI_new = []
  new_data = ptrarr(length)
  count = 0.0
  for j=0,n_elements(irradiance.spect20)-1 do begin
    
  ;timestamp is the Sorce Day of SIM, it will go through each SD as the loop appends
    timestamp_1st = gps2sd((*irradiance.spect20[j]).timestamp[0]/1d6)

   ;finds the index of the SD within NRLSSI data that is closest to the timestamp
    TS_index_NRLSSI = min(abs(NRLSSI_TS_array - timestamp_1st), pos_NRLSSI_TS)

    ;results in array of NRLSSI and SIM interpolated on the same wavelength grid  
    pos_sim=where(wl_grid lt (*irradiance.spect20[j]).wavelength[-1] and wl_grid gt (*irradiance.spect20[j]).wavelength[0], count1)  
    if count1 lt 1.0 then continue
    NRLSSI_SSI_NEW = interpol(data.ssi[*,pos_NRLSSI_TS], data.wavelength, wl_grid)

    ;convolve the interpolated NRLSSI2 irradiance array using the kernel found
    convolved_NRLSSI = convol(NRLSSI_SSI_new, kernel_y,/edge_mirror)

     if n_elements(pos_sim) ge elements then begin  
      ;with SPlINE -> Program caused arithmetic error: Floating divide by 0/Floating illegal operand   
       SIM_SSI_new = interpol((*irradiance.spect20[j]).irradiance, (*irradiance.spect20[j]).wavelength, wl_grid[pos_sim])
       
       new_data[count]=ptr_new({wavelength: wl_grid[pos_sim], irradiance_NRLSSI_conv:convolved_NRLSSI[pos_sim], irradiance_NRLSSI_int:NRLSSI_SSI_new[pos_sim], $
        irradiance_SIM_int:SIM_SSI_new, SIM_timestamp:timestamp_1st, NRLSSI_timestamp:NRLSSI_TS_array[pos_NRLSSI_TS], f_function:dblarr(n_elements(pos_sim))} )
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

