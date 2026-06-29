# =============================================================================
# MaxentRSystemGUI.R - Main GUI using tcltk
# =============================================================================

library(tcltk)
library(tcltk2)

create_maxent_gui <- function() {
  
  main_window <- tktoplevel()
  tkwm.title(main_window, "R Automated Maxent-R System")
  tkwm.geometry(main_window, "750x700")
  
  title_frame <- tkframe(main_window)
  title_label <- tklabel(title_frame, text = "R Automated Maxent-R System",
                         font = tkfont.create(family = "Arial", size = 16, weight = "bold"))
  tkpack(title_label, pady = 10)
  tkpack(title_frame, side = "top", fill = "x")
  
  notebook <- tkwidget(main_window, "ttk::notebook")
  
  # TAB 1: Input Data
  tab1 <- tkframe(notebook)
  tkadd(notebook, tab1, text = "  Input Data  ")
  
  occ_frame <- ttklabelframe(tab1, text = "Occurrence Data", padding = 15)
  occ_path_frame <- tkframe(occ_frame)
  occ_label <- tklabel(occ_path_frame, text = "Occurrence File:")
  tkpack(occ_label, side = "left", padx = 5)
  occ_entry <- tkentry(occ_path_frame, width = 40)
  tkpack(occ_entry, side = "left", padx = 5, fill = "x", expand = TRUE)
  browse_btn <- tkbutton(occ_path_frame, text = "Browse...", command = function() {
    fp <- tclvalue(tkgetOpenFile())
    if (nchar(fp) > 0) {
      tkdelete(occ_entry, 0, "end")
      tkinsert(occ_entry, 0, fp)
    }
  })
  tkpack(browse_btn, side = "left", padx = 5)
  tkpack(occ_path_frame, fill = "x", pady = 5)
  
  use_raw_gbif <- tclVar("0")
  raw_check <- tkcheckbutton(occ_frame, text = "Use raw GBIF data", variable = use_raw_gbif)
  tkpack(raw_check, anchor = "w", pady = 5)
  
  use_iucn <- tclVar("0")
  iucn_check <- tkcheckbutton(occ_frame, text = "Use IUCN data (SHP)", variable = use_iucn)
  tkpack(iucn_check, anchor = "w", pady = 5)
  
  shp_frame <- tkframe(occ_frame)
  shp_label <- tklabel(shp_frame, text = "SHP path:")
  tkpack(shp_label, side = "left", padx = 5)
  shp_entry <- tkentry(shp_frame, width = 35)
  tkpack(shp_entry, side = "left", padx = 5, fill = "x", expand = TRUE)
  shp_browse <- tkbutton(shp_frame, text = "Browse...", command = function() {
    fp <- tclvalue(tkgetOpenFile(filetypes = "{{Shapefile} {.shp}}"))
    if (nchar(fp) > 0) {
      tkdelete(shp_entry, 0, "end")
      tkinsert(shp_entry, 0, fp)
    }
  })
  tkpack(shp_browse, side = "left", padx = 5)
  tkpack(shp_frame, fill = "x", pady = 5)
  
  tkpack(occ_frame, fill = "both", expand = TRUE, padx = 20, pady = 20)
  
  # TAB 2: Geographic Extent
  tab2 <- tkframe(notebook)
  tkadd(notebook, tab2, text = "  Geographic Extent  ")
  
  clip_frame <- ttklabelframe(tab2, text = "Clipping Strategy", padding = 15)
  clip_to_extent <- tclVar("1")
  
  coord_frame <- tkframe(clip_frame)
  max_row <- tkframe(coord_frame)
  tkpack(tklabel(max_row, text = "Max Lat:"), side = "left", padx = 5)
  max_lat_entry <- tkentry(max_row, width = 10, state = "disabled")
  tkpack(max_lat_entry, side = "left", padx = 5)
  tkpack(tklabel(max_row, text = "Max Lon:"), side = "left", padx = 15)
  max_lon_entry <- tkentry(max_row, width = 10, state = "disabled")
  tkpack(max_lon_entry, side = "left", padx = 5)
  tkpack(max_row, pady = 5)
  
  min_row <- tkframe(coord_frame)
  tkpack(tklabel(min_row, text = "Min Lat:"), side = "left", padx = 5)
  min_lat_entry <- tkentry(min_row, width = 10, state = "disabled")
  tkpack(min_lat_entry, side = "left", padx = 5)
  tkpack(tklabel(min_row, text = "Min Lon:"), side = "left", padx = 15)
  min_lon_entry <- tkentry(min_row, width = 10, state = "disabled")
  tkpack(min_lon_entry, side = "left", padx = 5)
  tkpack(min_row, pady = 5)
  
  clip_check <- tkcheckbutton(clip_frame, text = "Clip to extent of points",
                               variable = clip_to_extent,
                               command = function() {
                                 state <- if(tclvalue(clip_to_extent) == "0") "normal" else "disabled"
                                 tkconfigure(max_lat_entry, state = state)
                                 tkconfigure(max_lon_entry, state = state)
                                 tkconfigure(min_lat_entry, state = state)
                                 tkconfigure(min_lon_entry, state = state)
                               })
  tkpack(clip_check, anchor = "w", pady = 5)
  tkpack(coord_frame, fill = "x", pady = 10)
  tkpack(clip_frame, fill = "x", padx = 20, pady = 10)
  
  # TAB 3: Model Parameters
  tab3 <- tkframe(notebook)
  tkadd(notebook, tab3, text = "  Model Parameters  ")
  
  var_frame <- ttklabelframe(tab3, text = "Variable Selection", padding = 15)
  
  vi_frame <- tkframe(var_frame)
  tkpack(tklabel(vi_frame, text = "Variable Importance Threshold:"), side = "left", padx = 5)
  var_imp_val <- tclVar("0.05")
  vi_slider <- tkscale(vi_frame, from = 0, to = 1, orient = "horizontal",
                       length = 250, variable = var_imp_val, resolution = 0.01, showvalue = FALSE)
  tkpack(vi_slider, side = "left", padx = 5)
  vi_label <- tklabel(vi_frame, textvariable = var_imp_val, font = tkfont.create(weight = "bold"))
  tkpack(vi_label, side = "left", padx = 5)
  tkpack(vi_frame, fill = "x", pady = 5)
  
  nv_frame <- tkframe(var_frame)
  tkpack(tklabel(nv_frame, text = "Number of Important Variables:"), side = "left", padx = 5)
  num_vars_entry <- tkentry(nv_frame, width = 10)
  tkinsert(num_vars_entry, 0, "10")
  tkpack(num_vars_entry, side = "left", padx = 5)
  tkpack(nv_frame, fill = "x", pady = 5)
  
  rt_frame <- tkframe(var_frame)
  tkpack(tklabel(rt_frame, text = "Required TIFs (comma-sep, optional):"), side = "left", padx = 5)
  req_tif_entry <- tkentry(rt_frame, width = 30)
  tkpack(req_tif_entry, side = "left", padx = 5)
  tkpack(rt_frame, fill = "x", pady = 5)
  
  tkpack(var_frame, fill = "x", padx = 20, pady = 10)
  
  ms_frame <- ttklabelframe(tab3, text = "Model Selection", padding = 15)
  
  crit_frame <- tkframe(ms_frame)
  tkpack(tklabel(crit_frame, text = "Selection Criterion:"), side = "left", padx = 5)
  selection_criterion <- tclVar("AICc")
  tkpack(tkradiobutton(crit_frame, text = "AIC", variable = selection_criterion, value = "AIC"), side = "left", padx = 5)
  tkpack(tkradiobutton(crit_frame, text = "AICc", variable = selection_criterion, value = "AICc"), side = "left", padx = 5)
  tkpack(tkradiobutton(crit_frame, text = "BIC", variable = selection_criterion, value = "BIC"), side = "left", padx = 5)
  tkpack(crit_frame, fill = "x", pady = 5)
  
  corr_frame <- tkframe(ms_frame)
  tkpack(tklabel(corr_frame, text = "Correlation Threshold (r):"), side = "left", padx = 5)
  corr_val <- tclVar("0.80")
  corr_slider <- tkscale(corr_frame, from = 0, to = 1, orient = "horizontal",
                         length = 250, variable = corr_val, resolution = 0.01, showvalue = FALSE)
  tkpack(corr_slider, side = "left", padx = 5)
  corr_label <- tklabel(corr_frame, textvariable = corr_val, font = tkfont.create(weight = "bold"))
  tkpack(corr_label, side = "left", padx = 5)
  tkpack(corr_frame, fill = "x", pady = 5)
  
  beta_frame <- tkframe(ms_frame)
  tkpack(tklabel(beta_frame, text = "Max Beta:"), side = "left", padx = 5)
  max_beta_entry <- tkentry(beta_frame, width = 8)
  tkinsert(max_beta_entry, 0, "5.0")
  tkpack(max_beta_entry, side = "left", padx = 5)
  tkpack(tklabel(beta_frame, text = "Beta Increment:"), side = "left", padx = 15)
  beta_inc_entry <- tkentry(beta_frame, width = 8)
  tkinsert(beta_inc_entry, 0, "0.5")
  tkpack(beta_inc_entry, side = "left", padx = 5)
  tkpack(beta_frame, fill = "x", pady = 5)
  
  rep_frame <- tkframe(ms_frame)
  tkpack(tklabel(rep_frame, text = "Replicates (final model):"), side = "left", padx = 5)
  reps_entry <- tkentry(rep_frame, width = 10)
  tkinsert(reps_entry, 0, "10")
  tkpack(reps_entry, side = "left", padx = 5)
  tkpack(rep_frame, fill = "x", pady = 5)
  
  tkpack(ms_frame, fill = "x", padx = 20, pady = 10)
  
  # TAB 4: Step Selection
  tab4 <- tkframe(notebook)
  tkadd(notebook, tab4, text = "  Step Selection  ")
  
  info_lab <- tklabel(tab4, 
    text = "Select which steps to execute. Note: Some steps depend on previous steps.",
    justify = "left", foreground = "blue")
  tkpack(info_lab, anchor = "w", padx = 20, pady = 10)
  
  step_names <- c(
    "Step 1-2: Get location data and clip to extent",
    "Step 3: Remove duplicates",
    "Step 4-5: Define spatial extent / clip rasters",
    "Step 6: Run global model",
    "Step 7: Identify top variables",
    "Step 8: Calculate correlation",
    "Step 9: Generate permutations",
    "Step 10: Run all permutations",
    "Step 11: Identify top model",
    "Step 12: Optimize regularization",
    "Step 13: Run final model"
  )
  step_keys <- c("step1","step3","step4","step6","step7","step8",
                 "step9","step10","step11","step12","step13")
  step_vars <- list()
  
  btn_frame <- tkframe(tab4)
  tkpack(tkbutton(btn_frame, text = "Select All", command = function() {
    for (v in step_vars) tclvalue(v) <- "1"
  }), side = "left", padx = 5)
  tkpack(tkbutton(btn_frame, text = "Deselect All", command = function() {
    for (v in step_vars) tclvalue(v) <- "0"
  }), side = "left", padx = 5)
  tkpack(btn_frame, anchor = "w", padx = 20, pady = 10)
  
  for (i in seq_along(step_names)) {
    step_vars[[step_keys[i]]] <- tclVar("1")
    chk <- tkcheckbutton(tab4, text = step_names[i], variable = step_vars[[step_keys[i]]])
    tkpack(chk, anchor = "w", padx = 30, pady = 2)
  }
  
  tkpack(notebook, fill = "both", expand = TRUE, padx = 10, pady = 10)
  
  # Progress
  prog_frame <- tkframe(main_window)
  prog_label <- tklabel(prog_frame, text = "Ready", font = tkfont.create(size = 9))
  tkpack(prog_label, anchor = "w", padx = 20)
  progress_bar <- tkwidget(prog_frame, "ttk::progressbar", length = 700, mode = "determinate", maximum = 14)
  tkpack(progress_bar, padx = 20, pady = 5)
  tkpack(prog_frame, side = "top", fill = "x")
  
  update_progress <- function(step, name) {
    tkconfigure(progress_bar, value = step)
    tkconfigure(prog_label, text = sprintf("Step %d: %s", step, name))
    tcl("update")
  }
  
  # Buttons
  button_panel <- tkframe(main_window)
  
  process_btn <- tkbutton(button_panel, text = "Process", 
                          font = tkfont.create(size = 11, weight = "bold"),
                          command = function() {
    occ_file <- tclvalue(tkget(occ_entry))
    if (nchar(occ_file) == 0) {
      tkmessageBox(title = "Error", message = "Please specify the occurrence file.", icon = "error")
      return()
    }
    
    config <- WorkflowConfiguration()
    config$step1 <- StepConfig$Step1Config(
      enabled = tclvalue(step_vars$step1) == "1",
      occurrence_file = occ_file,
      use_raw_gbif = tclvalue(use_raw_gbif) == "1",
      use_iucn_shp = tclvalue(use_iucn) == "1",
      shp_path = tclvalue(tkget(shp_entry))
    )
    config$step3 <- StepConfig$Step3Config(enabled = tclvalue(step_vars$step3) == "1")
    
    clip_ext <- tclvalue(clip_to_extent) == "1"
    config$step4 <- StepConfig$Step4Config(
      enabled = tclvalue(step_vars$step4) == "1",
      clip_to_extent = clip_ext,
      max_lat = if (clip_ext) NA else as.numeric(tclvalue(tkget(max_lat_entry))),
      max_lon = if (clip_ext) NA else as.numeric(tclvalue(tkget(max_lon_entry))),
      min_lat = if (clip_ext) NA else as.numeric(tclvalue(tkget(min_lat_entry))),
      min_lon = if (clip_ext) NA else as.numeric(tclvalue(tkget(min_lon_entry)))
    )
    config$step6 <- StepConfig$Step6Config(enabled = tclvalue(step_vars$step6) == "1")
    config$step7 <- StepConfig$Step7Config(
      enabled = tclvalue(step_vars$step7) == "1",
      threshold = as.numeric(tclvalue(var_imp_val)),
      top_count = as.integer(tclvalue(tkget(num_vars_entry)))
    )
    config$step8 <- StepConfig$Step8Config(enabled = tclvalue(step_vars$step8) == "1")
    config$step9 <- StepConfig$Step9Config(
      enabled = tclvalue(step_vars$step9) == "1",
      correlation_threshold = as.numeric(tclvalue(corr_val))
    )
    
    req_tifs_str <- tclvalue(tkget(req_tif_entry))
    req_tifs <- if (nchar(req_tifs_str) > 0) trimws(strsplit(req_tifs_str, ",")[[1]]) else character(0)
    
    config$step10 <- StepConfig$Step10Config(
      enabled = tclvalue(step_vars$step10) == "1", required_tifs = req_tifs)
    config$step11 <- StepConfig$Step11Config(
      enabled = tclvalue(step_vars$step11) == "1",
      selection_criterion = tclvalue(selection_criterion))
    config$step12 <- StepConfig$Step12Config(
      enabled = tclvalue(step_vars$step12) == "1",
      max_beta = as.numeric(tclvalue(tkget(max_beta_entry))),
      beta_increment = as.numeric(tclvalue(tkget(beta_inc_entry))),
      selection_criterion = tclvalue(selection_criterion))
    config$step13 <- StepConfig$Step13Config(
      enabled = tclvalue(step_vars$step13) == "1",
      replicates = as.integer(tclvalue(tkget(reps_entry))))
    
    answer <- tkmessageBox(title = "Confirm",
                           message = "Begin processing with these settings?",
                           icon = "question", type = "okcancel")
    if (tclvalue(answer) != "ok") return()
    
    result <- tryCatch({
      ProjectManager$runMethodology(config, update_progress)
      "success"
    }, error = function(e) paste("Error:", e$message))
    
    if (result == "success") {
      tkmessageBox(title = "Complete", message = "Workflow completed successfully!", icon = "info")
    } else {
      tkmessageBox(title = "Error", message = result, icon = "error")
    }
  })
  tkpack(process_btn, side = "left", padx = 5)
  
  help_btn <- tkbutton(button_panel, text = "Help", command = function() {
    help_text <- paste(
      "R Automated Maxent-R System Help\n\n",
      "INPUT DATA TAB:\n",
      "  - Specify your species occurrence CSV file\n",
      "  - Check Use raw GBIF data if it is untrimmed GBIF data\n",
      "  - Optionally clip to IUCN SHP boundary\n\n",
      "GEOGRAPHIC EXTENT TAB:\n",
      "  - Auto-clip to point extent, or specify bounds manually\n\n",
      "MODEL PARAMETERS TAB:\n",
      "  - Variable importance threshold (0-1)\n",
      "  - Number of top variables to consider\n",
      "  - Correlation threshold for combinations\n",
      "  - Selection criterion: AIC, AICc, or BIC\n",
      "  - Beta multiplier optimization range\n\n",
      "STEP SELECTION TAB:\n",
      "  - Toggle individual workflow steps on/off\n",
      sep = "")
    
    hw <- tktoplevel()
    tkwm.title(hw, "Help")
    ht <- tktext(hw, width = 70, height = 25, wrap = "word")
    tkinsert(ht, "end", help_text)
    tkconfigure(ht, state = "disabled")
    tkpack(ht, padx = 10, pady = 10)
  })
  tkpack(help_btn, side = "left", padx = 5)
  
  exit_btn <- tkbutton(button_panel, text = "Exit", command = function() tkdestroy(main_window))
  tkpack(exit_btn, side = "left", padx = 5)
  
  tkpack(button_panel, side = "bottom", pady = 15)
  tkfocus(main_window)
  invisible(main_window)
}

