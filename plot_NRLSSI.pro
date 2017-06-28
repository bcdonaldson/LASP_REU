;+
;Author: Laura Sandoval, Bailey Donaldson
;
;
;STEPS:
;(1) reads in all of the files needed
;(2) adjusts the NRLSSI data to the correct SD units
;(3) extracts the data for the SD the user chooses from SIM data
;(4) the NRLSSI data is interpolated so it is on the same wvl grid as SIM
;(5) The kernel is found
;(6) The interpolated NRLSSI data is convolved
;(7) Everything is plotted
;(8) Plots Sorce day vs. Irradiance at constant Wavelength 
;
;***********************
;
function plot_NRLSSI, Inst_mode_ID

  constant_var = ''
  read, constant_var, PROMPT = 'Which is kept constant: Wl or SD? '


;***************STEP 1****************
;
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
  endif

  if Inst_mode_ID eq 43 then begin
    filename='/Users/bado8599/Desktop/REU/uncorr_irr/uva_irrad_23_uncorr.sav'
    titlename='UV'
    wl_grid = Standard_wave_43.stdwavelength
    WL_range = where(wl_grid le 306.0 and wl_grid ge 240.0)
  endif

  if Inst_mode_ID eq 44 then begin
    filename= '/Users/bado8599/Desktop/REU/uncorr_irr/ira_irrad_23_uncorr.sav'
    titlename='IR'
    wl_grid = Standard_wave_44.stdwavelength
    WL_range = where(wl_grid le 1600.0 and wl_grid ge 950.0)
  endif

;restores the SIM file and renames the variable it contains to the name irraidance
  restore, filename
  title_name = FILE_BASENAME(filename, '.sav')
  assign_irr = 'irradiance = ' + title_name
  void = execute(assign_irr)


;****************STEP 2****************

;adjusts the NRLSSI data.time to be in SIM sorce days, renames it NRLSSI_TS_array
  nrlssi_day0 = ymd2sd(1610, 01, 01)
  data.time += nrlssi_day0
  postime = where(data.time ge 0.0)
  NRLSSI_TS_array = data.time[postime]



  ;****************STEP 3****************
if strlowcase(constant_var) eq 'sd' then begin
  
;for now the Sorce Day will be kept constant at 453.6784
  read, user_time, PROMPT = 'Choose a timestamp (Missions Days): '
;  user_time = 453.6784

;forms an array containing the Sorce Days used in SIM measurements
  SIM_TS_array = []
  for i=0, n_elements(irradiance.spect20) - 1 do begin
    SIM_TS_array = [SIM_TS_array, gps2sd((*irradiance.spect20[i]).timestamp[0]/1d6)]
  endfor

;finds the index of the SD within SIM data that is closest to the SD the user chose$
;and calls it pos_SIM_TS (position of SIM time stamp)
  TS_index_SIM = min(abs(SIM_TS_array - user_time), pos_SIM_TS)

;finds the index of the SD within NRLSSI data that is closest to the SD the user chose$
;and calls it pos_NRLSSI_TS (position of SIM time stamp)
  TS_index_NRLSSI = min(abs(NRLSSI_TS_array - user_time), pos_NRLSSI_TS)



;***************STEP 4****************

;will be an array containing the NRLSSI irradiance measurements on the SIM wavelength grid
  NRLSSI_SSI_new = []
  SIM_SSI_new = []

  NRLSSI_SSI_NEW = interpol(data.ssi[*,pos_NRLSSI_TS], data.wavelength, wl_grid)
  SIM_SSI_new = interpol((*irradiance.spect20[pos_SIM_TS]).irradiance, (*irradiance.spect20[pos_SIM_TS]).wavelength, wl_grid, /spline)


;****************STEP 5****************

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


;****************STEP 6****************
;convolve the interpolated NRLSSI2 irradiance array using the kernel found
  convolved_NRLSSI = convol(NRLSSI_SSI_new, kernel_y, /edge_mirror)


;****************STEP 7****************


;keeps the lineplot from turning red;
  device, decompose = 0
  
  plot_range = where((*irradiance.spect20[pos_SIM_TS]).wavelength[0] lt wl_grid and (*irradiance.spect20[pos_SIM_TS]).wavelength[-1] gt wl_grid, pos)

;plot interpolated SIM at a constant SD, wavelength vs irradiance
  lineplot, wl_grid[plot_range], SIM_SSI_new[plot_range], $
    ptitle = titlename + ':NRLSSI2 and SIM Comparison', xtitle = 'Wavelength', ytitle = 'Irradiance (W/m^2*nm)', $
    title = 'Interpolated SIM at SD: ' + strtrim(string(SIM_TS_array[pos_SIM_TS], format = '(F0.2)'),0), charsize = 1.5

;plots interpolated NRLSSI at a constant SD, wavelength vs. irradiance
  lineplot, wl_grid, NRLSSI_SSI_new, title = 'Interpolated NRLSSI at SD: ' + $
    strtrim(string(NRLSSI_TS_array[pos_NRLSSI_TS], format = '(F0.2)'),0)
  
;plots convolved and interpolated NRLSSI, wavelength vs. irradiance
 lineplot, wl_grid, convolved_NRLSSI, title = 'Convolved and Interpolated NRLSSI at SD: ' + $
    strtrim(string(NRLSSI_TS_array[pos_NRLSSI_TS], format = '(F0.2)'),0)
  
endif

;****************STEP 8****************
 if strlowcase(constant_var) eq 'wl' then begin
   wl = 0
   wl_pos=0
   user_wl =0

   ;finds wavelength closest to the one the user chose (user_wl)
   READ, user_wl, PROMPT = 'Choose wavelength: '
   mn = min(abs(data.wavelength - user_wl), pos)

   ;another plotting method, first line keeps the plot from turning red
   device, decompose = 0
   lineplot, NRLSSI_TS_array, data.ssi[pos,*], xtitle='Sorce day', ytitle='Irradiance (W/m^2*nm)', $
     title = 'NRLSSI data Sorce day vs. Irradiance at: ' + user_wl

  endif

  return, 0
end
