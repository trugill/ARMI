\# MaxEnt Command-Line Flags Reference



This document describes the command-line flags available in MaxEnt.



| Category | Flag | Abbrv | Type | Default | Description |

|----------|------|:----:|------|---------|-------------|

| \*\*Input\*\* | samplesfile | s | File | — | CSV containing species occurrence records used to train the model. |

| \*\*Input\*\* | environmentallayers | e | File/Directory | — | Directory of environmental raster layers or an SWD CSV file. |

| \*\*Input\*\* | projectionlayers | j | File/Directory | — | Environmental layers used when projecting the trained model. |

| \*\*Input\*\* | testsamplesfile | T | File | — | Independent occurrence records used only for evaluation. |

| \*\*Input\*\* | biasfile |  | File | — | Raster describing sampling bias for selecting background points. |

| \*\*Output\*\* | outputdirectory | o | Directory | — | Directory where all output files are written. |

| \*\*Output\*\* | outputformat |  | String | cloglog | Transformation used for prediction values (cloglog, logistic, cumulative, raw). |

| \*\*Output\*\* | outputfiletype |  | String | asc | Raster file format used for prediction grids. |

| \*\*Output\*\* | outputgrids | x | Boolean | true | Write raster prediction grids. |

| \*\*Output\*\* | pictures |  | Boolean | true | Create PNG images of prediction maps. |

| \*\*Output\*\* | plots |  | Boolean | true | Generate plots included in the HTML report. |

| \*\*Output\*\* | responsecurves | P | Boolean | false | Generate response curve figures. |

| \*\*Output\*\* | responsecurvesexponent |  | Boolean | false | Plot the raw exponent instead of transformed predictions. |

| \*\*Output\*\* | writeplotdata |  | Boolean | false | Export numerical response curve data. |

| \*\*Output\*\* | writebackgroundpredictions |  | Boolean | false | Export predictions at all background locations. |

| \*\*Output\*\* | writeclampgrid |  | Boolean | true | Produce a raster showing the effects of clamping. |

| \*\*Output\*\* | writemess |  | Boolean | true | Produce a MESS map identifying novel environments. |

| \*\*Output\*\* | perspeciesresults |  | Boolean | false | Write a separate `maxentResults.csv` for each species. |

| \*\*Output\*\* | appendtoresultsfile |  | Boolean | false | Append new results to the existing `maxentResults.csv` instead of overwriting it. |

| \*\*Output\*\* | logfile |  | String | maxent.log | Name of the diagnostic log file. |

| \*\*Modeling\*\* | autofeature | A | Boolean | true | Automatically choose feature classes based on sample size. |

| \*\*Modeling\*\* | linear | l | Boolean | true | Enable linear features. |

| \*\*Modeling\*\* | quadratic | q | Boolean | true | Enable quadratic features. |

| \*\*Modeling\*\* | hinge | h | Boolean | true | Enable hinge features. |

| \*\*Modeling\*\* | product | p | Boolean | true | Enable product (interaction) features. |

| \*\*Modeling\*\* | threshold |  | Boolean | false | Enable threshold features. |

| \*\*Modeling\*\* | betamultiplier | b | Double | 1.0 | Multiplies all automatic regularization parameters. Higher values produce smoother, less complex models. |

| \*\*Modeling\*\* | beta\_lqp |  | Double | -1.0 | Regularization multiplier for linear, quadratic, and product features. Negative values use automatic estimation. |

| \*\*Modeling\*\* | beta\_hinge |  | Double | -1.0 | Regularization multiplier for hinge features. |

| \*\*Modeling\*\* | beta\_threshold |  | Double | -1.0 | Regularization multiplier for threshold features. |

| \*\*Modeling\*\* | beta\_categorical |  | Double | -1.0 | Regularization multiplier for categorical variables. |

| \*\*Modeling\*\* | maximumiterations | m | Integer | 500 | Maximum number of optimization iterations during model fitting. |

| \*\*Modeling\*\* | convergencethreshold | c | Double | 1.0E-5 | Stop optimization when improvement in log loss falls below this threshold. |

| \*\*Modeling\*\* | defaultprevalence |  | Double | 0.5 | Assumed prevalence used when converting raw predictions to logistic or cloglog output. |

| \*\*Modeling\*\* | l2lqthreshold |  | Integer | 10 | Minimum sample size before quadratic features are automatically enabled. |

| \*\*Modeling\*\* | hingethreshold |  | Integer | 15 | Minimum sample size before hinge features are automatically enabled. |

| \*\*Modeling\*\* | lq2lqptthreshold |  | Integer | 80 | Minimum sample size before product and threshold features are automatically enabled. |

| \*\*Evaluation\*\* | randomtestpoints | X | Integer | 0 | Percentage of occurrence records randomly reserved for testing. |

| \*\*Evaluation\*\* | replicates |  | Integer | 1 | Number of replicate model runs. |

| \*\*Evaluation\*\* | replicatetype |  | String | crossvalidate | Replicate strategy: `crossvalidate`, `bootstrap`, or `subsample`. |

| \*\*Evaluation\*\* | jackknife | J | Boolean | false | Measure variable importance by training models with variables omitted and in isolation. |

| \*\*Evaluation\*\* | randomseed |  | Boolean | false | Use a different random seed for each run, resulting in different train/test splits and background samples. |

| \*\*Projection\*\* | extrapolate |  | Boolean | true | Allow predictions outside the environmental range encountered during training. |

| \*\*Projection\*\* | doclamp |  | Boolean | true | Clamp environmental values outside the training range during projection. |

| \*\*Projection\*\* | fadebyclamping |  | Boolean | false | Reduce predictions in regions where clamping substantially affects results. |

| \*\*Projection\*\* | applythresholdrule |  | String | — | Apply a named threshold rule to generate a binary prediction raster. |

| \*\*Data Processing\*\* | removeduplicates |  | Boolean | true | Remove duplicate occurrence records occupying the same raster cell (or identical coordinates for SWD data). |

| \*\*Data Processing\*\* | addsamplestobackground | d | Boolean | true | Add occurrence environmental combinations not already represented in the background sample. |

| \*\*Data Processing\*\* | addallsamplestobackground |  | Boolean | false | Add every occurrence record to the background, regardless of duplication. |

| \*\*Data Processing\*\* | maximumbackground | MB | Integer | 10000 | Maximum number of background points sampled during model training. |

| \*\*Data Processing\*\* | allowpartialdata |  | Boolean | false | Allow occurrence records with missing environmental values during training. |

| \*\*Data Processing\*\* | nodata | n | Integer | -9999 | Value interpreted as missing data in SWD files. |

| \*\*Performance\*\* | threads |  | Integer | 1 | Number of processor threads used during model training. |

| \*\*Performance\*\* | cache |  | Boolean | true | Cache ASCII raster files as `.mxe` files for faster subsequent access. |

| \*\*Interface\*\* | visible | z | Boolean | true | Display the MaxEnt graphical user interface. |

| \*\*Interface\*\* | autorun | a | Boolean | false | Automatically begin model training when MaxEnt starts. |

| \*\*Interface\*\* | tooltips |  | Boolean | true | Display explanatory tooltips in the graphical interface. |

| \*\*Interface\*\* | warnings |  | Boolean | true | Display warning dialogs during execution. |

| \*\*Interface\*\* | askoverwrite | r | Boolean | true | Ask before overwriting existing output files. |

| \*\*Interface\*\* | skipifexists | S | Boolean | false | Skip modeling species whose output files already exist. |

| \*\*Interface\*\* | logscale |  | Boolean | true | Use a logarithmic color scale in prediction map images. |

| \*\*Interface\*\* | adjustsampleradius |  | Integer | 0 | Increase or decrease the size of occurrence points displayed on prediction maps. |

| \*\*Advanced\*\* | togglelayertype | t | String | — | Toggle selected environmental layers between continuous and categorical types. |

| \*\*Advanced\*\* | togglespeciesselected | E | String | — | Toggle selection of species whose names match a specified prefix. |

| \*\*Advanced\*\* | togglelayerselected | N | String | — | Toggle selection of environmental layers whose names match a specified prefix. |

| \*\*Advanced\*\* | prefixes |  | Boolean | true | Allow toggle commands to match prefixes instead of exact names. |

| \*\*Advanced\*\* | verbose | v | Boolean | false | Print detailed diagnostic information useful for debugging. |

