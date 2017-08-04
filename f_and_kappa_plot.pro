;Author: Bailey Donaldson
;
;Purpose: graph calculated f function and calculated kappa on a surface and contour plot
;
function save_f_and_kappa, Inst_mode_ID

  if Inst_mode_ID eq 41 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/VIS_data_F_function.sav'
    data = vis_data
    restore, '/Users/bado8599/Desktop/REU/sav_files/VIS_daily_kappa.sav'
    type = 'VIS'
    restore, '/Users/bado8599/Desktop/REU/sav_files/prism_deg_VIS.sav'
  endif

  if Inst_mode_ID eq 43 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/UV_data_F_function.sav'
    data = uv_data
    restore, '/Users/bado8599/Desktop/REU/sav_files/UV_daily_kappa.sav'
    type = 'UV'
    restore, '/Users/bado8599/Desktop/REU/sav_files/prism_deg_UV.sav'
  endif

  if Inst_mode_ID eq 44 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/IR_data_F_function.sav'
    data = ir_data
    restore, '/Users/bado8599/Desktop/REU/sav_files/IR_daily_kappa.sav'
    type = 'IR'
    restore, '/Users/bado8599/Desktop/REU/sav_files/prism_deg_IR.sav'
  endif

  ;puts kappa and f function into a one-dimensional array so it is easier to plot

  nelem=0L
  for i=0, n_elements(data)-1 do nelem+=n_elements((*data[i]).wavelength)

  wavelength = dblarr(nelem)
  sorce_day = dblarr(nelem)
  f_function = dblarr(nelem)
  kappa = dblarr(nelem)
  prism = dblarr(nelem)

  p0=0L

  for i=0, n_elements(data)-1 do begin
    p1= p0 + n_elements((*data[i]).wavelength)-1L
    wavelength[p0:p1] = (*data[i]).wavelength
    sorce_day[p0:p1] = replicate((*data[i]).SIM_timestamp,n_elements((*data[i]).wavelength))
    f_function[p0:p1] = (*data[i]).f_function
    kappa[p0:p1] = (*daily_kappa[i]).kappa
    prism[p0:p1] = (*prism_deg[i]).prism_deg_avg_kappa
    p0=p1+1   
  endfor

  ;f
  ;sset=bspline_iterfit(sorce_day,f_function, maxiter=0, requiren=10, bkspace=10)
  ;yfitA=bspline_valu(sorce_day ,sset)

  ;k
  ;resistant_mean,(irradA-yfit),3.0,mean,goodvec=keep0
  ; sset1=bspline_iterfit(sorce_day, kappa, maxiter=0, requiren=10, bkspace=5)
  ; yfitA1=bspline_valu(kappa, sset1)

  ;  sset=bspline_iterfit(timeA[keep0],irradA[keep0],maxiter=0,requiren=10,bkspace=5)
  ;  yfitA=bspline_valu(timeA,sset)
  ;  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  data_array = {wavelength:wavelength, sorce_day:sorce_day, f_function:f_function, daily_kappa:kappa, prism: prism}

  print, 'max f', max(f_function)
  print, 'max k', max(kappa)
  print, 'min f', min(f_function)
  k2 = where(kappa ne -1)
  print, 'min k', min(kappa[k2])
  
  if Inst_mode_ID eq 41 then begin
    save, data_array, filename = '/Users/bado8599/Desktop/REU/sav_files/VIS_3D_plot_array.sav'
  endif
  if Inst_mode_ID eq 43 then begin
    save, data_array, filename = '/Users/bado8599/Desktop/REU/sav_files/UV_3D_plot_array.sav'
  endif
  if Inst_mode_ID eq 44 then begin
    save, data_array, filename = '/Users/bado8599/Desktop/REU/sav_files/IR_3D_plot_array.sav'
  endif

end
;started vis at 3:27pm

function plot_f_kappa1, Inst_mode_ID

  if Inst_mode_ID eq 41 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/VIS_3D_plot_array.sav'
    k1 = where(data_array.f_function gt -1 and data_array.f_function lt 7 and data_array.wavelength gt 315 and data_array.wavelength lt 920)
    k2 = where(data_array.daily_kappa lt 0.0013 and data_array.daily_kappa gt -.0003 and data_array.wavelength gt 315 and data_array.wavelength lt 920)
    k3 = where(data_array.wavelength gt 315 and data_array.wavelength lt 920 and data_array.prism lt 1.2 and data_array.prism gt .7)
    type = 'VIS'
  endif
  if Inst_mode_ID eq 43 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/UV_3D_plot_array.sav'
    k1=where(data_array.f_function lt 4 and data_array.f_function gt -3 and data_array.wavelength lt 305 and data_array.wavelength gt 240)
    k2=where(data_array.daily_kappa lt .005 and data_array.daily_kappa gt -.003 and data_array.wavelength lt 305 and data_array.wavelength gt 240)
    k3 = where(data_array.wavelength lt 305 and data_array.wavelength gt 240 and data_array.prism lt 1.5)
    type = 'UV'
  endif
  if Inst_mode_ID eq 44 then begin
    restore, '/Users/bado8599/Desktop/REU/sav_files/IR_3D_plot_array.sav'
    k1=where(data_array.f_function lt 80 and data_array.f_function gt 0 and data_array.wavelength lt 1600 and data_array.wavelength gt 950)
    k2=where(data_array.daily_kappa lt .0017 and data_array.daily_kappa gt -0.0001 and data_array.wavelength lt 1600 and data_array.wavelength gt 950)
    type = 'IR'
    k3=where(data_array.wavelength lt 1600 and data_array.wavelength gt 950 and data_array.prism lt 1.05 and data_array.prism gt .70)
  endif

  ;contour plots
  ct = colortable(72, /reverse)
  
;;;;;;;;;;;;;;;;;;;;;;;;;;F FUNCTION;;;;;;;;;;;;;;;;;;;;;;;;
  w1 = data_array.wavelength[k1]
  sd1 = data_array.sorce_day[k1]
  f1= data_array.f_function[k1]
 
  ;epsilon in vis is 1.0, uv is 0.55, and ir is
 ; grid_input, w1, sd1, f1, x1, y1, z1, duplicates='First', epsilon=1.4
 ; contour1 = contour(z1, x1, y1, /irr, /fill, rgb_table =ct, xtitle = 'Wavelength (nm)', ytitle = 'Sorce Day', title = 'F function in the '+type)
  ;cb = colorbar(title = 'F function')
  ;surface1 = surface(z1, x1, y1, /irr, texture_image = colors, color='orange red', xtitle = 'Wavelength (nm)', ytitle='Sorce Day', title='F function in the '+type)
  
;;;;;;;;;;;;;;;;KAPPA;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
  w2 = data_array.wavelength[k2]
  sd2 = data_array.sorce_day[k2]
  dk2 = data_array.daily_kappa[k2]

  ;epsilon in uv is 0.8, vis is 1.2, ir is .85
  grid_input, w2, sd2, dk2, x2, y2, z2, duplicates='First', epsilon=.85
 ; contour2 = contour(z2, x2, y2, /irr, /fill, rgb_table=ct, xtitle = 'Wavelength (nm)', ytitle = 'Sorce Day', title = 'Daily Kappa in the '+type)
 ; cb = colorbar(title = 'Kappa')
  surface2= surface(z2, x2, y2, /irregular, color='olive drab', xtitle='Wavelength (nm)', ytitle='Sorce Day', title='Daily Kappa in the '+type)

;;;;;;;;;;;;;;;;PRISM DEG::::::::::::::::::
;w1 = data_array.wavelength[k3]
;sd1 = data_array.sorce_day[k3]
;p1= data_array.prism[k3]
;
;grid_input, w1, sd1, p1, x1, y1, z1, duplicates='First', epsilon=2.0
;contour1 = contour(z1, x1, y1, /irr, /fill, rgb_table =ct, xtitle = 'Wavelength (nm)', ytitle = 'Sorce Day', title = 'Prism deg in the '+type)
;cb = colorbar(title = 'Prsim deg')
end