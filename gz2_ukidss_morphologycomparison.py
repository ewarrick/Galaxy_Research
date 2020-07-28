##Creating Plots comparing GZ2 and UKIDSS Morphologies##

#Import necessary packages
import astropy
import numpy as np
import cv2
from astropy.visualization import hist
#Coleman's plotting notebook
from astropy.utils.data import get_pkg_data_filename
from astropy.io import fits
import matplotlib as mpl
from matplotlib import pyplot as plt
import numpy as np
plt.style.use(astropy.visualization.astropy_mpl_style)

#Open fits file of cross-matched sample
hdul=fits.open("UKIDS_GZ2_NSAv1_0_1match.fits")

#Get fits file data
data=hdul[1].data
cols=hdul[1].columns

#######Building the Sample Selection for GZ2 and UKIDSS########

#sample selection gz2
spiralsgz2=data[data["t01_smooth_or_features_a02_features_or_disk_debiased_gz2"]>0.430]
notedgegz2=spiralsgz2[spiralsgz2["t02_edgeon_a05_no_debiased_gz2"]>0.715]
Nnotedgegz2=notedgegz2[notedgegz2['t02_edgeon_a05_no_count']>20] #make a histogram in the notedge subsample
visiblespiralgz2=Nnotedgegz2[Nnotedgegz2["t04_spiral_a08_spiral_debiased_gz2"]>0.5] #take out Nnotedge
petromaggz2=visiblespiralgz2[visiblespiralgz2["PETROMAG_MR"]<-19]
petromag1gz2=petromaggz2[petromaggz2["REDSHIFT"]<0.035]
g=notedgegz2['SERSIC_ABSMAG'][:,3]
r=notedgegz2['SERSIC_ABSMAG'][:,4]
g_rgz2=g-r #optical color of g-r #optical color
M_rgz2=notedgegz2['PETROMAG_MR'] #petro mag k-corrected
zgz2=notedgegz2['REDSHIFT'] #red shift
pbargz2=notedgegz2['t03_bar_a06_bar_debiased_gz2']
#galaxies and barred gals
discs_gz2=np.where((M_rgz2 < -20.2) & (zgz2 > 0.01) & (zgz2 < 0.06)) #disc/galaxies filtering
bars_gz2=np.where((M_rgz2 < -20.2) & (zgz2 > 0.01) & (zgz2 < 0.06) & (pbargz2 > 0.4)) #bar filtering
#set up bins and plot for gz2 data
min=0.2
max=1.0
step=0.06
bins=np.arange(min,max,step)
plotbins=np.arange(min,max-step,step)
#set up the histograms
galhist_gz2,outbins=np.histogram(g_rgz2[discs_gz2],bins)
barhist_gz2,outbins=np.histogram(g_rgz2[bars_gz2],bins)
#Bar Fraction of optical data
frac_gz2=(barhist_gz2/galhist_gz2) #gz2 bar fraction
e_fbar_gz2=np.sqrt(barhist_gz2)/galhist_gz2 #error bars for gz2 data
#arm winding and bulge size gz2
w_mag_gz2=(0.5)*petromag1gz2["t10_arms_winding_a29_medium_debiased_gz2"]+1.0*petromag1gz2["t10_arms_winding_a28_tight_debiased_gz2"]
b_avg_gz2=0.2*petromag1gz2["t05_bulge_prominence_a11_just_noticeable_debiased_gz2"]+0.8*petromag1gz2["t05_bulge_prominence_a12_obvious_debiased_gz2"]+1.0*petromag1gz2["t05_bulge_prominence_a13_dominant_debiased_gz2"]

#sample selection ukidss
spirals=data[data["t00_smooth_or_features_a1_features_debiased_rh_ukidss"]>0.430]
notedge=spirals[spirals["t01_disk_edge_on_a1_no_debiased_rh_ukidss"]>0.715]
Nnotedge=notedge[notedge['t01_disk_edge_on_count_weighted_ukidss']>20]#; changed on 7/5/20
visiblespiral=Nnotedge[Nnotedge["t03_spiral_a0_spiral_debiased_rh_ukidss"]>0.5]
petromag=visiblespiral[visiblespiral["PETROMAG_MR"]<-19]
petromag1=petromag[petromag["REDSHIFT"]<0.035]
#bar fraction of IR data
g=notedge['SERSIC_ABSMAG'][:,3]
r=notedge['SERSIC_ABSMAG'][:,4]
g_r=g-r #optical color of (g-r)
pbar = notedge['t02_bar_a0_bar_debiased_rh_ukidss']
M_r=notedge['PETROMAG_MR'] #petro mag k-corrected
z=notedge['REDSHIFT'] #red shift
#galaxies and barred gals
discs=np.where((M_r < -19.38) & (z < 0.06) & (z > 0.01))
bars=np.where((M_r < -19.38) & (z < 0.06) & (z > 0.01) & (pbar>0.4))
#create histograms/bins
min=0.2
max=1.0
step=0.06
bins=np.arange(min,max,step)
plotbins=np.arange(min,max-step,step)
#setting up the histograms
galhist,outbins=np.histogram(g_r[discs],bins)
barhist,outbins=np.histogram(g_r[bars],bins)
#bar fraction
frac=barhist/galhist #bar fraction ukidss data
e_fbar=np.sqrt(barhist)/galhist #error bars for UKIDSS data

#arm winding and bulge size ukidss
w_mag=(0.5)*petromag1["t09_arms_winding_a1_medium_debiased_rh_ukidss"]+1.0*petromag1["t09_arms_winding_a0_tight_debiased_rh_ukidss"]
b_avg=0.2*petromag1["t04_bulge_prominence_a1_just_noticeable_debiased_rh_ukidss"]+0.8*petromag1["t04_bulge_prominence_a2_obvious_debiased_rh_ukidss"]+1.0*petromag1["t04_bulge_prominence_a3_dominant_debiased_rh_ukidss"]

#######Plotting the Subfigures########
fig, axs = plt.subplots(2, 1, sharex=True)
# Remove horizontal space between axes
fig.subplots_adjust(hspace=0)

#Bar fraction and Errors for Bar Fraction
axs[0].errorbar(plotbins,frac_gz2,yerr=e_fbar_gz2,ls="solid", label="GZ2", marker="s")
axs[0].errorbar(plotbins,frac,yerr=e_fbar,ls="dashed", label="UKIDSS",marker="o")
#Tick marks
axs[0].set_yticks(np.arange(0.05,0.55,0.05))
#Axes limits
axs[0].set_ylim(0.05,0.55)
axs[0].set_xlim(0.2,0.9)
#Y Axis Label
axs[0].set_ylabel("Bar Fraction")
axs[0].grid(None)
#Galaxy Histogram
axs[1].hist(g_rgz2[discs_gz2],plotbins,label="gz2", histtype='step')
axs[1].hist(g_r[discs],plotbins,label="ukidss", histtype='step', ls='dashed')
#X Axis Limits
axs[1].set_xlim(0.2,0.9)
#Axis Labels
axs[1].set_ylabel("Number")
axs[1].set_xlabel("(g-r)")
axs[1].legend(loc='upper left')
axs[1].grid(None)
plt.style.use('seaborn-colorblind')
plt.show()

#Arm Winding vs Bulge Size UKIDSS Plot
plt.figure()
hb = plt.hexbin(b_avg, w_mag, gridsize=10, cmap='viridis')
plt.xlabel('Bulge Size (b_avg)')
plt.ylabel('Arm Winding (w_mag)')
plt.title("Bulge Size vs Arm Winding from UKIDSS")
cb=plt.colorbar(hb)
plt.show()
#Arm Winding vs Bulge Size GZ2 Plot
plt.figure()
hb = plt.hexbin(b_avg_gz2, w_mag_gz2, gridsize=10, cmap='viridis')
plt.xlabel('b_avg')
plt.ylabel('w_mag')
plt.title("Bulge Size vs Arm Winding from GZ2")
cb=plt.colorbar(hb)
plt.show()
