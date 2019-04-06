import os
import sys
sys.path.append('C:\Program Files\GDCM 2.8\lib')
import pydicom
from PIL import *
import subprocess
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm
import numpy as np
from scipy.signal import butter, filtfilt
from scipy.signal import savgol_filter
from scipy.stats import linregress
from scipy import signal
from scipy.signal import find_peaks,peak_prominences,peak_widths
import fitz

# np.set_printoptions(threshold=np.nan)

#*************FILTERING EXAMPLE**************************************
# def butter_lowpass(cutoff, fs, order=5):
#     nyq = 0.5 * fs
#     normal_cutoff = cutoff / nyq
#     b, a = butter(order, normal_cutoff, btype='low', analog=False)
#     return b, a
#
# def butter_lowpass_filtfilt(data, cutoff, fs, order=5):
#     b, a = butter_lowpass(cutoff, fs, order=order)
#     y = filtfilt(b, a, data)
#     return y

    # cutoff = 500
    # fs = 70000
    # amp_filt = butter_lowpass_filtfilt(amplitude, cutoff, fs)

    # N=100
    # amp_filt=np.convolve(amplitude, np.ones((N,)) / N, mode='full')
#*************FILTERING EXAMPLE**************************************

def running_mean(x, N):
    out = np.zeros_like(x, dtype=np.float64)
    dim_len = x.shape[0]
    for i in range(dim_len):
        if N%2 == 0:
            a, b = i - (N-1)//2, i + (N-1)//2 + 2
        else:
            a, b = i - (N-1)//2, i + (N-1)//2 + 1

        #cap indices to min and max indices
        a = max(0, a)
        b = min(dim_len, b)
        out[i] = np.mean(x[a:b])
    return out


#axial visualization and scrolling
def multi_slice_viewer(volume,dx,dy):
    # remove_keymap_conflicts({'j', 'k'})
    fig, ax = plt.subplots()
    ax.volume = volume
    ax.index = volume.shape[2] // 2
    print(ax.index)
    extent = (0, 0 + (volume.shape[1] * dx),
              0, 0 + (volume.shape[0] * dy))
    ax.imshow(volume[:,:,ax.index],extent=extent)
    ax.set_xlabel('x distance [mm]')
    ax.set_ylabel('y distance [mm]')
    ax.set_title("slice="+ str(ax.index))
    fig.suptitle('Axial view', fontsize=16)
    fig.canvas.mpl_connect('key_press_event', process_key_axial)

def process_key_axial(event):
    fig = event.canvas.figure
    ax = fig.axes[0]
    if event.key == 'j':
        previous_slice_axial(ax)
    elif event.key == 'k':
        next_slice_axial(ax)
    ax.set_title("slice="+str(ax.index))
    fig.canvas.draw()

def previous_slice_axial(ax):
    volume = ax.volume
    ax.index = (ax.index - 1) % volume.shape[2]  # wrap around using %
    print(ax.index,volume.shape[2])
    ax.images[0].set_array(volume[:,:,ax.index])

def next_slice_axial(ax):
    volume = ax.volume
    ax.index = (ax.index + 1) % volume.shape[2]
    print(ax.index,volume.shape[2])
    ax.images[0].set_array(volume[:,:,ax.index])




#work here
def minimize_junction(amplitude, peaks,peak_type,dx):
    print('number of peaks=',peaks)
    print('amplitude dimensions=', amplitude.shape[1])
    junctions_fig=[]

    amp_prev=0
    amp_filt_prev=0

    fig = plt.figure(figsize=(10, 6))# create the plot

    kk = 1  # counter for figure generation
    for j in range(0, amplitude.shape[1]-1):
        for k in range(j + 1, amplitude.shape[1]): #looping through remaining images
            amp_base_res = signal.savgol_filter(signal.resample(amplitude[:, j], 60), 21, 5)
            amp_overlay_res = signal.savgol_filter(signal.resample(amplitude[:, k], 60), 21, 5)
            peak1, _ = find_peaks(amp_base_res, prominence=5000)
            peak2, _ = find_peaks(amp_overlay_res, prominence=5000)
            print('peak1,peak2,diff', peak1, peak2, abs(peak2 - peak1))

            if abs(peak2 - peak1) < 18:  # if the two peaks are close together proceeed to analysis
                print(j, k, 'the data is contiguous... analizing', )
                cumsum_prev = 1e6
                # if j==1:
                #     exit(0)
                # print('peak=', j, peaks[j], peak_type[j])
                amp_base_res = amplitude[:, j]
                amp_overlay_res = amplitude[:, k]

                if peak_type[j] == 0:
                    inc = -1
                else:
                    inc = 1
                for i in range(0, inc * 38, inc * 1):
                    # print('i', i)
                    # for i in range(0, 1, 38):
                    x = np.linspace(0, 0 + (len(amp_base_res) * dx), len(amplitude),
                                    endpoint=False)  # definition of the distance axis
                    amp_overlay_res_roll = np.roll(amp_overlay_res, i)

                    # amplitude is the vector to analyze +-500 samples from the center
                    amp_tot = (amp_base_res[peaks[j] - 1000:peaks[j] + 1000] + amp_overlay_res_roll[
                                                                               peaks[j] - 1000:peaks[
                                                                                                   j] + 1000])  # divided by 2 to normalize
                    xsel = x[peaks[j] - 1000:peaks[j] + 1000]

                    amp_filt = running_mean(amp_tot, 281)
                    cumsum = np.sum(np.abs(amp_tot - amp_filt))

                    if cumsum > cumsum_prev:  # then we went too far
                        # print('peak=', j, 'i=', i + 1, "dx=", dx, "delta=", abs(i + 1) * dx, "cumsum=", cumsum,
                        #       "cumsum_prev=", cumsum_prev, '<-final')
                        # junctions_fig.append(plt.figure(figsize=(10, 6)))
                        ax = fig.add_subplot(amplitude.shape[1] - 1, 1, kk)

                        # plt.plot(-1*(amp_base_res+amp_overlay_res_roll))
                        ax.plot(amp_prev)
                        ax.plot(amp_filt_prev)
                        if kk == 1:
                            ax.set_title('Minimization result')
                        if kk == amplitude.shape[1] - 1:  # if we reach the final plot the add the x axis label
                            ax.set_xlabel('x distance [mm]')

                        ax.set_ylabel('amplitude')
                        ax.annotate('delta=' + str(abs(i - inc * 1) * dx) + ' mm', xy=(2, 1), xycoords='axes fraction',
                                    xytext=(.35, .10))

                        # plt.show()

                        kk = kk + 1
                        break
                    else:
                        amp_prev = amp_tot
                        amp_filt_prev = amp_filt
                        cumsum_prev = cumsum
                        # print('i=', i, "dx=", dx, "delta=", abs(i) * dx, "cumsum=", cumsum, "cumsum_prev=", cumsum_prev)

                # plt.show()


            else:
                print(j, k, 'the data is not contiguous finding another curve in dataset')










    return fig











# this subroutine aims to find the peaks
def peak_find(ampl_resamp,dx):
    # print('type=',type(ampl_resamp),'shape=',ampl_resamp.shape[1])
    # peaks=np.zeros(ampl_resamp.shape[1]-1,dtype=int)
    peak_figs=[]
    peaks=[]
    peak_type=[]
    print('number of curves=',ampl_resamp.shape[1])
    for j in range(0, ampl_resamp.shape[1]-1):
        print('j=',j)
        amp_base_res = signal.savgol_filter(signal.resample(ampl_resamp[:,j],60),21,5)
        for k in range(j+1,ampl_resamp.shape[1]):
            amp_overlay_res = signal.savgol_filter(signal.resample(ampl_resamp[:,k],60),21,5)
            peak1, _ = find_peaks(amp_base_res, prominence=5000)
            peak2, _ = find_peaks(amp_overlay_res, prominence=5000)
            print('peak find',peak1,peak2,abs(peak2-peak1))
            # plt.figure()
            # plt.plot(amp_base_res)
            # plt.plot(amp_overlay_res)
            # plt.show()

            if abs(peak2-peak1)<18:  # if the two peaks are separated the two fields are not adjacent.
                amp_peak = (ampl_resamp[:,j] + ampl_resamp[:,k]) / 2
                x = np.linspace(0, 0 + (len(amp_peak) * dx / 10), len(amp_peak),
                                endpoint=False)  # definition of the distance axis
                # plt.figure()
                # plt.plot(amp_peak)
                # plt.show()

                # #now we need to find if we are going to find a peak or a through
                # while True:  # example of infinite loops using try and except to catch only numbers
                #     line = input('Is there a peak or through in the amplitude profile [peak(p)/through(t)]> ')
                #     try:
                #         ##        if line == 'done':
                #         ##            break
                #         poption = str(line.lower())
                #         if poption.startswith(('p', 'pe', 'peak', 't', 'thro', 'through')):
                #             break
                #
                #     except:
                #         print('Please enter a valid option:')
                # if poption.startswith(('p', 'pe', 'peak')):
                #     peak, _ = find_peaks(amp_peak, prominence=2000)
                # elif poption.startswith(('t', 'thro', 'through')):
                #     peak, _ = find_peaks(-amp_peak, prominence=2000)
                peak, _ = find_peaks(amp_peak, prominence=2000)
                if len(peak)!=1:
                    peak, _ = find_peaks(-amp_peak, prominence=2000)
                    peaks.append(peak[0])
                    peak_type.append(0)
                else:
                    peak_type.append(1)
                    peaks.append(peak[0])

                print('peak=', j, peaks)
                print(j,k,'the data is contiguous')
                fig = plt.figure(figsize=(10, 6))
                plt.plot(x, amp_peak, label='Total amplitude profile')
                plt.plot(x[peak], amp_peak[peak], "x", label='Peaks detected')
                plt.ylabel('amplitude [a.u.]')
                plt.xlabel('distance [mm]')
                plt.legend()
                fig.suptitle('Junctions', fontsize=16)
                peak_figs.append(fig)
                # plt.show()

            else:
                print(j,k,'the data is not contiguous finding another curve in dataset')





    print('peaks=',peaks)
    return peaks, peak_type, peak_figs































#this subroutine will merge the two jaws into a single image and display a graph of the overlap
def merge_view_vert(volume,dx,dy):


    #creating merged volume
    merge_vol = np.zeros((volume.shape[0], volume.shape[1]))

    #creating vector for processing along cols (one row)
    amplitude = np.zeros((volume.shape[1],volume.shape[2]))  # 1 if it is vertical 0 if the bars are horizontal

    x = np.linspace(0,0+(volume.shape[1]*dx),volume.shape[1],endpoint=False)#definition of the distance axis
    # x = np.arange(0,)#definition of the distance axis

    #merging the two images together
    ampl_resamp = np.zeros(((volume.shape[1])*10, volume.shape[2]))
    amp_peak = np.zeros((volume.shape[1]) * 10)


    for slice in tqdm(range(0,volume.shape[2])):
        merge_vol = merge_vol + volume[:, :, slice]
        amplitude[:,slice]=volume[int(volume.shape[0] / 2), :, slice]
        ampl_resamp[:,slice]=signal.resample(amplitude[:,slice],int(len(amplitude))*10) #resampling the amplitude vector
        amp_peak=amp_peak+ampl_resamp[:,slice]/volume.shape[2]








    fig, ax = plt.subplots(nrows=2,squeeze=True,sharex=True, figsize=(6,8))

    extent = (0, 0 + (volume.shape[1] * dx),
              0, 0 + (volume.shape[0] * dy))

    print(len(x),merge_vol.shape)
    ax[0].imshow(merge_vol,extent=extent,aspect='auto')
    # ax[0].set_aspect('equal', 'box')
    ax[0].set_xlabel('x distance [mm]')
    ax[0].set_ylabel('y distance [mm]')

    ax[1].plot(x,amplitude)
    ax[1].set_xlabel('x distance [mm]')
    ax[1].set_ylabel('amplitude')
    # ax.set_title("slice=" + str(ax.index))
    fig.suptitle('Merged volume', fontsize=16)




    peaks, peak_type,peak_figs = peak_find(ampl_resamp,dx)
    junctions_figs=minimize_junction(ampl_resamp,peaks,peak_type, dx / 10)


    return fig, peak_figs, junctions_figs













































#this subroutine will merge the 4 jaws and analyze the two upper and two lower pairs
#each jaw is 6cm in length (60mm)
def merge_view_horz(volume,dx,dy):


    # creating merged volume
    merge_vol = np.zeros((volume.shape[0], volume.shape[1]))

    # creating vector for processing along cols (one row)
    amplitude = np.zeros((volume.shape[0],volume.shape[2]))  # 1 if it is vertical 0 if the bars are horizontal

    x = np.linspace(0, 0 + (volume.shape[0] * dy), volume.shape[0], endpoint=False)  # definition of the distance axis
    # x = np.arange(0,) #definition of the distance axis

    # merging the two images together
    ampl_resamp = np.zeros(((volume.shape[0]) * 10, volume.shape[2]))
    amp_peak = np.zeros((volume.shape[0]) * 10)


    for slice in tqdm(range(0, volume.shape[2])):
        merge_vol = merge_vol + volume[:, :, slice]
        amplitude[:, slice] = volume[:,int(volume.shape[1] / 2), slice]
        ampl_resamp[:, slice] = signal.resample(amplitude[:, slice], int(len(amplitude)) * 10)  # resampling the amplitude vector
        amp_peak = amp_peak + ampl_resamp[:, slice] / volume.shape[2]


    fig, ax = plt.subplots(nrows=2, squeeze=True, sharey=True, figsize=(6,8))


    extent = (0, 0 + (volume.shape[1] * dx),
              0, 0 + (volume.shape[0] * dy))

    print(len(x), merge_vol.shape)
    ax[0].imshow(merge_vol, extent=extent, aspect='auto')
    # ax[0].set_aspect('equal', 'box')
    ax[0].set_xlabel('x distance [mm]')
    ax[0].set_ylabel('y distance [mm]')

    ax[1].plot(amplitude,x,label='Amplitude profile')
    # ax[1].set_xlabel('x distance [mm]')
    # ax[1].set_ylabel('amplitude')
    ax[1].set_xlabel('amplitude')
    ax[1].set_ylabel('y distance [mm]')
    ax[1].legend()
    # ax.set_title("slice=" + str(ax.index))
    fig.suptitle('Merged volume', fontsize=16)

    peaks, peak_type,peak_figs = peak_find(ampl_resamp,dx)
    junction_figs = minimize_junction(ampl_resamp,peaks,peak_type, dx / 10)




    return fig,peak_figs, junction_figs

































































#this subroutine will merge the 4 jaws and analyze the two upper and two lower pairs
#each jaw is 6cm in length (60mm)
def merge_view_filtrot(volume,dx,dy):



    # creating merged volume
    merge_vol = np.zeros((volume.shape[0], volume.shape[1]))

    # creating vector for processing along cols (one row)
    amplitude = np.zeros((volume.shape[0],volume.shape[2]))  # 1 if it is vertical 0 if the bars are horizontal

    x = np.linspace(0, 0 + (volume.shape[0] * dy), volume.shape[0], endpoint=False)  # definition of the distance axis
    # x = np.arange(0,) #definition of the distance axis

    # merging the two images together
    ampl_resamp = np.zeros(((volume.shape[0]) * 10, volume.shape[2]))
    amp_peak = np.zeros((volume.shape[0]) * 10)


    for slice in tqdm(range(0, volume.shape[2])):
        merge_vol = merge_vol + volume[:, :, slice]
        amplitude[:, slice] = volume[:,int(volume.shape[1] / 2), slice]
        ampl_resamp[:, slice] = signal.resample(amplitude[:, slice], int(len(amplitude)) * 10)  # resampling the amplitude vector
        amp_peak = amp_peak + ampl_resamp[:, slice] / volume.shape[2]


    fig, ax = plt.subplots(nrows=2, squeeze=True, sharey=True, figsize=(6,8))


    extent = (0, 0 + (volume.shape[1] * dx),
              0, 0 + (volume.shape[0] * dy))

    print(len(x), merge_vol.shape)
    ax[0].imshow(merge_vol, extent=extent, aspect='auto')
    # ax[0].set_aspect('equal', 'box')
    ax[0].set_xlabel('x distance [mm]')
    ax[0].set_ylabel('y distance [mm]')

    ax[1].plot(amplitude,x,label='Amplitude profile')
    # ax[1].set_xlabel('x distance [mm]')
    # ax[1].set_ylabel('amplitude')
    ax[1].set_xlabel('amplitude')
    ax[1].set_ylabel('y distance [mm]')
    ax[1].legend()
    # ax.set_title("slice=" + str(ax.index))
    fig.suptitle('Merged volume', fontsize=16)

    peaks, peak_type,peak_figs = peak_find(ampl_resamp,dx)
    junction_figs = minimize_junction(ampl_resamp,peaks,peak_type, dx / 10)




    return fig,peak_figs, junction_figs













































# this routine anlyzes the volume and autodetect what type of analysis to carry on (x, Y, Field Rot)
def folder_analyze(volume):
    print(volume.shape)

    for slice in range(0,volume.shape[2]):
        # plt.figure()
        # plt.imshow(volume[:,:,slice])
        # plt.show()
        stack1 = np.sum(volume[:,:,slice],axis=0)
        maxstack1= np.max(stack1)

        stack2 = np.sum(volume[:,:,slice],axis=1)
        maxstack2= np.max(stack2)

        if maxstack2/maxstack1>1.5:
            return 2
        elif maxstack2/maxstack1<0.5:
            return 1
        else:
            return 3



































#the data loaded is not in HU we need to rescale it
# def get_pixels_hu(scans):
#     image = np.stack([s.pixel_array for s in scans])
#     # Convert to int16 (from sometimes int16),
#     # should be possible as values should always be low enough (<32k)
#     image = image.astype(np.int16)
#
#     # Set outside-of-scan pixels to 1
#     # The intercept is usually -1024, so air is approximately 0
#     image[image == -2000] = 0
#
#     # Convert to Hounsfield units (HU)
#     intercept = scans[0].RescaleIntercept
#     slope = scans[0].RescaleSlope
#
#     if slope != 1:
#         image = slope * image.astype(np.float64)
#         image = image.astype(np.int16)
#
#     image += np.int16(intercept)
#
#     return np.array(image, dtype=np.int16)

# interpolate the dataset so the spacing is equal in every direction
# def resample(image, scan, new_spacing=[1, 1, 1]): #this will resample the data to 1mmx1mmx1mm
#     # Determine current pixel spacing
#     spacing = map(float, ([scan[0].SliceThickness] + scan[0].PixelSpacing))
#     spacing = np.array(list(spacing))
#
#     resize_factor = spacing / new_spacing
#     new_real_shape = image.shape * resize_factor
#     new_shape = np.round(new_real_shape)
#     real_resize_factor = new_shape / image.shape
#     new_spacing = spacing / real_resize_factor
#
#     image = scipy.ndimage.interpolation.zoom(image, real_resize_factor)
#
#     return image, new_spacing


def read_dicom3D(dirname,poption):
    slice=0
    # lstFilesDCM = [] #empty list to store dicom files
    for subdir, dirs, files in os.walk(dirname):
        k = 0
        for file in tqdm(sorted(files)):
            print('filename=',file)
            if poption.startswith(('y', 'yeah', 'yes')):
                subprocess.call(["gdcmconv", "-w", dirname + file, os.path.splitext(dirname + file)[0] + "_decomp" + ".dcm"])
                dataset = pydicom.dcmread(os.path.splitext(dirname + file)[0] + "_decomp" + ".dcm")
                if k==0:
                    ArrayDicom = np.zeros((dataset.Rows, dataset.Columns), dtype=dataset.pixel_array.dtype)
                    ArrayDicom = np.dstack((ArrayDicom, dataset.pixel_array))
                    # print("slice thickness [mm]=",dataset.SliceThickness)
                    SID=dataset.RTImageSID
                    dx=1/(SID*(1/dataset.ImagePlanePixelSpacing[0])/1000)
                    dy=1/(SID*(1/dataset.ImagePlanePixelSpacing[1])/1000)
                    print("pixel spacing row [mm]=", dx)
                    print("pixel spacing col [mm]=", dy)
                else:
                    ArrayDicom=np.dstack((ArrayDicom, dataset.pixel_array))
            elif poption.startswith(('n', 'no', 'nope')):
                dataset = pydicom.dcmread(dirname + file)
                if k==0:
                    ArrayDicom = np.zeros((dataset.Rows,dataset.Columns,0), dtype=dataset.pixel_array.dtype)
                    ArrayDicom = np.dstack((ArrayDicom, dataset.pixel_array))
                    # print("slice thickness [mm]=", dataset.SliceThickness)
                    SID = dataset.RTImageSID
                    dx=1/(SID*(1/dataset.ImagePlanePixelSpacing[0])/1000)
                    dy=1/(SID*(1/dataset.ImagePlanePixelSpacing[1])/1000)
                    print("pixel spacing row [mm]=", dx)
                    print("pixel spacing col [mm]=", dy)
                else:
                    ArrayDicom = np.dstack((ArrayDicom, dataset.pixel_array))
            print(k)
            k=k+1




    # we have loaded the image but the spacing is different in the x,y and z directions
    #spacing = map(float, ([scan[0].SliceThickness] + scan[0].PixelSpacing))
    #spacing = np.array(list(spacing))
    print(type(ArrayDicom),np.shape(ArrayDicom))


    doption = folder_analyze(ArrayDicom)
    # print('doption=',doption)
    # exit(0)



    multi_slice_viewer(ArrayDicom,dx,dy)
    if doption==1:
        fig,peak_figs, junctions_figs=merge_view_vert(ArrayDicom, dx, dy)

    elif doption==2:
        fig,peak_figs,junctions_figs= merge_view_horz(ArrayDicom, dx, dy)
        # print('peak_figs********************************************************=', len(peak_figs),peak_figs)
    else:
        fig, peak_figs, junctions_figs = merge_view_filtrot(ArrayDicom, dx, dy)


    with PdfPages('multipage_pdf.pdf') as pdf:
        pdf.savefig(fig)
        for i in range(0,len(peak_figs)):
            pdf.savefig(peak_figs[i])
        # for j in range(0,len(junctions_figs)):
        pdf.savefig(junctions_figs)
        plt.close()



    if sys.platform=='linux':
        # Uncomment for Linux machine
        doc = fitz.open("/mnt/home/peter/Dropbox/PhDMedPhysi/scripts-medphys/multipage_pdf.pdf")  # open the PDF
        rect = fitz.Rect(0, 0, 100, 32)  # where to put image: use upper left corner

        for page in doc:
            page.insertImage(rect, filename="/mnt/home/peter/Dropbox/PhDMedPhysi/scripts-medphys/ahs-logo.png")

        doc.saveIncr()

    elif sys.platform == 'windows':
        # Uncomment for Windows machine
        doc = fitz.open("E:\zDropbox\Dropbox\PhDMedPhysi\scripts-medphys\multipage_pdf.pdf")  # open the PDF
        rect = fitz.Rect(0, 0, 100, 32)  # where to put image: use upper left corner

        for page in doc:
            page.insertImage(rect, filename="E:\zDropbox\Dropbox\PhDMedPhysi\scripts-medphys\\ahs-logo.png")

        doc.saveIncr()


    # plt.show(block=True)  #this shows all the figures!! (most important plt.show)





    # Normal mode:
    print()
    print("Directory folder.........:", dirname)
    print("Storage type.....:", dataset.SOPClassUID)
    print()

    pat_name = dataset.PatientName
    display_name = pat_name.family_name + ", " + pat_name.given_name
    print("Patient's name...:", display_name)
    print("Patient id.......:", dataset.PatientID)
    print("Modality.........:", dataset.Modality)
    print("Study Date.......:", dataset.StudyDate)
    print("Gantry angle......", dataset.GantryAngle)
    #
    # if 'PixelData' in dataset:
    #     rows = int(dataset.Rows)
    #     cols = int(dataset.Columns)
    #     print("Image size.......: {rows:d} x {cols:d}, {size:d} bytes".format(
    #         rows=rows, cols=cols, size=len(dataset.PixelData)))
    #     if 'PixelSpacing' in dataset:
    #         print("Pixel spacing....:", dataset.PixelSpacing)
    #
    # # use .get() if not sure the item exists, and want a default value if missing
    # print("Slice location...:", dataset.get('SliceLocation', "(missing)"))



    # print(type(dataset))
    # print(dataset.pixel_array)
    # # exit(0)
    # # plot the image using matplotlib
    # plt.imshow(dataset.pixel_array, cmap=plt.cm.bone)
    # plt.xlabel('X Distance from reference point')
    # plt.ylabel('Y Distance from reference point')
    # plt.show(block=True)



try:
    dirname = str(sys.argv[1])
    print(dirname)
except:
    print('Please enter a valid filename')
    print("Use the following command to run this script")
    print("python test_pydicom3D.py \"[dirname]\"")


while True:  # example of infinite loops using try and except to catch only numbers
    line = input('Are the files compressed [yes(y)/no(n)]> ')
    try:
        ##        if line == 'done':
        ##            break
        poption = str(line.lower())
        if poption.startswith(('y', 'yeah', 'yes','n','no','nope')):
            break

    except:
        print('Please enter a valid option:')

read_dicom3D(dirname,poption)

