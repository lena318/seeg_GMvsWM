#!/bind/bash

image_directory=/Users/andyrevell/mount/USERS/arevell/papers/paper005/data_raw/figure1_data
FIGURE_DIRECTORY=/Users/andyrevell/mount/USERS/arevell/papers/paper005/paper005/figures/figure1
#coordinates
x=0.24419277906417847
y=-5.40368777513504
z=24.70941925048828

zoom=113.17829558699817

#fsleyes --scene 3d --worldLoc $x $y $z --displaySpace world --cameraRotation -70.97  0.80 -0.04 --zoom $zoom --hideLegend --lightPos 250.24408847093582 244.5962079167366 24.70941925048828 --offset 0.0 0.0 --hideCursor --bgColour 1 1 1 --fgColour 0.0039215686274509665 0.0039215686274509665 0.0039215686274509665 --cursorColour 0.0 1.0 0.0 --colourBarLocation top --colourBarLabelSide top-left --performance 3 --movieSync ${image_directory}/sub-RID0365_T01_CT_to_T00_mprageANTs.nii.gz --name "sub-RID0365_T01_CT_to_T00_mprageANTs_1" --overlayType volume --alpha 100.0 --brightness 35.5468750082764 --contrast 50.97807649312962 --cmap greyscale --negativeCmap greyscale --displayRange 199.76316921582065 4214.658704428504 --clippingRange 199.76316921582065 3111.95 --gamma 0.0 --cmapResolution 256 --interpolation none --numSteps 500 --blendFactor 0.5496261433725306 --smoothing 0 --resolution 100 --numInnerSteps 10 --clipMode intersection --volume 0
fsleyes render --outfile $FIGURE_DIRECTORY/sub-RID0365_T01_CT_to_T00_mprageANTs.png -sz 400 400 --scene 3d --worldLoc 0.24419277906417847 -5.40368777513504 24.70941925048828 --displaySpace world --cameraRotation -37.00 -0.08  0.85 --zoom 113.17829558699817 --hideLegend --lightPos 250.24408847093582 244.5962079167366 24.70941925048828 --offset -0.12 0.0 --hideCursor --bgColour 1 1 1 --fgColour 0.0039215686274509665 0.0039215686274509665 0.0039215686274509665 --cursorColour 0.0 1.0 0.0 --colourBarLocation top --colourBarLabelSide top-left --performance 3 --movieSync ${image_directory}/sub-RID0365_T01_CT_to_T00_mprageANTs.nii.gz --name "sub-RID0365_T01_CT_to_T00_mprageANTs" --overlayType volume --alpha 100.0 --brightness 35.5468750082764 --contrast 50.97807649312962 --cmap greyscale --negativeCmap greyscale --displayRange 199.76316921582065 4214.658704428504 --clippingRange 199.76316921582065 3111.95 --gamma 0.3497009 --cmapResolution 1024 --interpolation spline --invert --numSteps 500 --blendFactor 0.4091417688050172 --smoothing 0 --resolution 100 --numInnerSteps 10 --clipMode intersection --volume 0





