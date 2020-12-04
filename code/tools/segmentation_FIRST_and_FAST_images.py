#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 11:23:28 2020

@author: arevell
"""

#combine subcortical segmentation (FIST) and regular segmentation (FAST)
import os

def first_and_fast_segmentation(ifname_preop3T, ifname_T00, ofpath_segmentation_sub_ID):
    
    #copy files to temporary
    cmd = "cp {0}.nii.gz {1}".format(ifname_preop3T, ofpath_segmentation_sub_ID); print(cmd); os.system(cmd)
    
    #copy files to temporary
    cmd = "cp {0}.nii.gz {1}".format(ifname_T00, ofpath_segmentation_sub_ID); print(cmd); os.system(cmd)
    
    
    ifname_preop3T =  os.path.join(ofpath_segmentation_sub_ID, os.path.basename(ifname_preop3T))
    ifname_T00 =  os.path.join(ofpath_segmentation_sub_ID, os.path.basename(ifname_T00))
    
    
    #Orient all images to standard RAS
    ##Orient the preop3T image
    cmd = "fslreorient2std {0}.nii.gz {1}.nii.gz".format(ifname_preop3T, "{0}_std".format(ifname_preop3T)); print(cmd); os.system(cmd)
    ##Orient the T00  image
    cmd = "fslreorient2std {0}.nii.gz {1}.nii.gz".format(ifname_T00, "{0}_std".format(ifname_T00)); print(cmd); os.system(cmd)
   

    
    #brain extraction of preop. Brain Extraction parameters matter to get FIRST to work properly
    #setting g and f parameters for bet
    f = 0.5; g = -0.3
    if (os.path.basename(ifname_preop3T) == "sub-RID0194_ses-preop3T_acq-3D_T1w"):
        f = 0.5; g= -0.3
    if (os.path.basename(ifname_preop3T) == "sub-RID0278_ses-preop3T_acq-3D_T1w"):
        f = 0.3; g= -0.4
    if (os.path.basename(ifname_preop3T) == "sub-RID0320_ses-preop3T_acq-3D_T1w"):
        f = 0.5; g= -0.3
    if (os.path.basename(ifname_preop3T) == "sub-RID0420_ses-preop3T_acq-3D_T1w"):
        f = 0.5; g= -0.5
    if (os.path.basename(ifname_preop3T) == "sub-RID0502_ses-preop3T_acq-3D_T1w"):
        f = 0.3; g= -0.4
    if (os.path.basename(ifname_preop3T) == "sub-RID0508_ses-preop3T_acq-3D_T1w"):
        f = 0.4; g= -0.5
    
    #default:f=0.5, g=-0.3 #RID0194:f ,g  #RID0320:f=0.5'g=-0.3    #RID0420:f=0.5,g=-0.5    #RID502:f=0.3,g=-0.4  #RID0508 f=0.4 g=-0.5;
    cmd = "bet {0}_std.nii.gz {0}_std_bet.nii.gz -f {1} -g {2}".format(ifname_preop3T, f, g); print(cmd); os.system(cmd)
    
    #brain extraction of T00
    cmd = "bet {0}_std.nii.gz {0}_std_bet.nii.gz -f 0.15".format(ifname_T00); print(cmd); os.system(cmd)
    
    #subcortical segmentation (FIRST)
    cmd = "run_first_all -i {0}_std_bet.nii.gz -o {0}_std_bet_subcort.nii.gz -b -v".format(ifname_preop3T); print(cmd); os.system(cmd)
    
    #seg of preop img  (FAST)
    cmd = "fast -n 3 -H 0.25 -t 1 -v {0}_std_bet.nii.gz".format(ifname_preop3T); print(cmd); os.system(cmd)
    
    
def register_preop3T_to_T00(ifname_preop3T, ifname_T00, ofpath_segmentation_sub_ID, ofbase_flirt, ofbase_fnirt):
    
    ifname_preop3T =  os.path.join(ofpath_segmentation_sub_ID, os.path.basename(ifname_preop3T))
    ifname_T00 =  os.path.join(ofpath_segmentation_sub_ID, os.path.basename(ifname_T00))
    
    #linear reg of preop to T00 space
    cmd = "flirt -in {0}_std_bet.nii.gz -ref {1}_std_bet.nii.gz -dof 12 -out {2} -omat {2}.mat -v".format(ifname_preop3T, ifname_T00, ofbase_flirt)
    print(cmd)
    os.system(cmd)
    
    #non linear reg of preop to T00 space
    cmd = "fnirt --in={0}_std.nii.gz --ref={1}_std.nii.gz --aff={2}.mat --iout={3} -v --cout={3}_coef --fout={3}_warp".format(ifname_preop3T, ifname_T00, ofbase_flirt, ofbase_fnirt)
    print(cmd)
    os.system(cmd)
    
    
def applywarp_to_combined_first_fast(ofname_FIRST_FAST_COMBINED, ifname_T00, ofbase_fnirt, ofpath_segmentation_sub_ID, ofname_FIRST_FAST_COMBINED_to_T00):
    ifname_T00 =  os.path.join(ofpath_segmentation_sub_ID, os.path.basename(ifname_T00))
    #warp combined_FIRST_and_FAST image to T00 space
    cmd = "applywarp -i {0} -r {1}_std.nii.gz -w {2}_warp.nii.gz --interp=nn -o {3}".format(ofname_FIRST_FAST_COMBINED, ifname_T00, ofbase_fnirt, ofname_FIRST_FAST_COMBINED_to_T00)
    print(cmd)
    os.system(cmd)

"""


#For RID0309 because images are not in correct orientation

##Orient the preop3T image
#fslreorient2std sub-RID0309_ses-preop3T_acq-3D_T1w.nii.gz sub-RID0309_ses-preop3T_acq-3D_T1w_std.nii.gz
#fslreorient2std sub-RID0309_ses-preop3T_acq-3D_T1w_bet.nii.gz sub-RID0309_ses-preop3T_acq-3D_T1w_bet_std.nii.gz
##Orient the T00  image
#fslreorient2std sub-RID0309_T00_mprage.nii.gz sub-RID0309_T00_mprage_std.nii.gz
#fslreorient2std sub-RID0309_T00_mprage_bet.nii.gz sub-RID0309_T00_mprage_bet_std.nii.gz
#flirt -in sub-RID0309_ses-preop3T_acq-3D_T1w_bet_std.nii.gz -ref sub-RID0309_T00_mprage_bet_std.nii.gz -dof 12 -out sub-RID0309_preop3T_to_T00 -omat sub-RID0309_preop3T_to_T00.mat -v
#ifname_preop3T_std = "{0}_std".format(ifname_preop3T)
#ifname_T00_std = "{0}_std".format(ifname_T00)
#cmd = "fnirt --in={0}.nii.gz --ref={1}.nii.gz --aff={2}.mat --iout={3} -v --cout={3}_coef --fout={3}_warp".format(ifname_preop3T_std, ifname_T00_std, ofbase_flirt, ofbase_fnirt)
#print(cmd)
#os.system(cmd)
#cmd = "applywarp -i {0} -r {1}.nii.gz -w {2}_warp.nii.gz --interp=nn -o {3}".format(ofname_FIRST_FAST_COMBINED, ifname_T00_std, ofbase_fnirt, ofname_FIRST_FAST_COMBINED_to_T00)
#print(cmd)
#os.system(cmd)


#For RID0194

##Orient the preop3T image
fslreorient2std sub-RID0194_ses-preop3T_acq-3D_T1w.nii.gz sub-RID0194_ses-preop3T_acq-3D_T1w_std.nii.gz
fslreorient2std sub-RID0194_ses-preop3T_acq-3D_T1w_bet.nii.gz sub-RID0194_ses-preop3T_acq-3D_T1w_bet_std.nii.gz
##Orient the T00  image
fslreorient2std sub-RID0194_T00_mprage.nii.gz sub-RID0194_T00_mprage_std.nii.gz
fslreorient2std sub-RID0194_T00_mprage_bet.nii.gz sub-RID0194_T00_mprage_bet_std.nii.gz
#Orient GM_WM
fslreorient2std sub-RID0194_ses-preop3T_acq-3D_T1w_bet_GM_WM_CSF.nii.gz sub-RID0194_ses-preop3T_acq-3D_T1w_bet_GM_WM_CSF_std.nii.gz
flirt -in sub-RID0194_ses-preop3T_acq-3D_T1w_bet_std.nii.gz -ref sub-RID0194_T00_mprage_bet_std.nii.gz -dof 12 -out sub-RID0194_preop3T_to_T00 -omat sub-RID0194_preop3T_to_T00.mat -v
ifname_preop3T_std = "{0}_std".format(ifname_preop3T)
ifname_T00_std = "{0}_std".format(ifname_T00)
cmd = "fnirt --in={0}.nii.gz --ref={1}.nii.gz --aff={2}.mat --iout={3} -v --cout={3}_coef --fout={3}_warp".format(ifname_preop3T_std, ifname_T00_std, ofbase_flirt, ofbase_fnirt)
print(cmd)
os.system(cmd)
ofname_FIRST_FAST_COMBINED = "{0}_std.nii.gz".format(os.path.splitext(os.path.splitext(ofname_FIRST_FAST_COMBINED)[0])[0])
cmd = "applywarp -i {0} -r {1}.nii.gz -w {2}_warp.nii.gz --interp=nn -o {3}".format(ofname_FIRST_FAST_COMBINED, ifname_T00_std, ofbase_fnirt, ofname_FIRST_FAST_COMBINED_to_T00)
print(cmd)
os.system(cmd)



#For RID0320

##Orient the preop3T image
fslreorient2std sub-RID0320_ses-preop3T_acq-3D_T1w.nii.gz sub-RID0320_ses-preop3T_acq-3D_T1w_std.nii.gz
fslreorient2std sub-RID0320_ses-preop3T_acq-3D_T1w_bet.nii.gz sub-RID0320_ses-preop3T_acq-3D_T1w_bet_std.nii.gz
##Orient the T00  image
fslreorient2std sub-RID0320_T00_mprage.nii.gz sub-RID0320_T00_mprage_std.nii.gz
fslreorient2std sub-RID0320_T00_mprage_bet.nii.gz sub-RID0320_T00_mprage_bet_std.nii.gz
flirt -in sub-RID0320_ses-preop3T_acq-3D_T1w_bet_std.nii.gz -ref sub-RID0320_T00_mprage_bet_std.nii.gz -dof 12 -out sub-RID0320_preop3T_to_T00 -omat sub-RID0320_preop3T_to_T00.mat -v
ifname_preop3T_std = "{0}_std".format(ifname_preop3T)
ifname_T00_std = "{0}_std".format(ifname_T00)
cmd = "fnirt --in={0}.nii.gz --ref={1}.nii.gz --aff={2}.mat --iout={3} -v --cout={3}_coef --fout={3}_warp".format(ifname_preop3T_std, ifname_T00_std, ofbase_flirt, ofbase_fnirt)
print(cmd)
os.system(cmd)
cmd = "applywarp -i {0} -r {1}.nii.gz -w {2}_warp.nii.gz --interp=nn -o {3}".format(ofname_FIRST_FAST_COMBINED, ifname_T00_std, ofbase_fnirt, ofname_FIRST_FAST_COMBINED_to_T00)
print(cmd)
os.system(cmd)




"""










